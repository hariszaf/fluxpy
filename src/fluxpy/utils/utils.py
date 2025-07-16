""" Functionalities to handle metabolic modeling analysis objects"""
import json
import numpy as np
import pandas as pd
import cobra
from typing import Union, Optional, Dict, List
from difflib import get_close_matches
from Levenshtein import distance as levenshtein_distance
from ..constants import *

NestedDataFrameType = pd.DataFrame


def extract_model_namespace(model):
    """
    Returns a model's namespace considering only ModelSEED and BiGG

    Args:
        model (cobra.Model, mandatory) -- the model to extract its namespace from

    Returns:
        namespace (str) -- a Literal["modelseed", "bigg"]
    """
    return "modelseed" if model.metabolites[0].id.startswith("cpd") else "bigg"

# %% Util functions: parsing

def nonzero_fluxes(fluxes:  Union[cobra.Solution, pd.Series]) -> None:
    """
    Returns nonzero fluxes from a pandas Series.

    This function takes a pandas Series containing flux values and returns a new Series with only the nonzero fluxes.

    Args:
        sol (pd.Series, mandatory): A pandas Series representing flux values.

    Returns:
        pd.Series: A pandas Series containing only the nonzero flux values from the input Series.
    """
    if not isinstance(fluxes, cobra.Solution) and not isinstance(fluxes, pd.Series):
        return ValueError("Provide a cobra.Solution or a pandas.Series element as fluxes.")
    if isinstance(fluxes, cobra.Solution):
        fluxes = fluxes.fluxes
    mask = fluxes.to_numpy().nonzero()
    return fluxes.iloc[mask]



def get_rxns_producing_consuming_met(met_id, model: cobra.Model = None, flux_vector: pd.Series = None):
    """
    Returns a DataFrame with the reaction IDs of those consuming and producing a given metabolite in a model.

    This function optimizes the given metabolic model to identify reactions involving a specified metabolite.
    It categorizes reactions into those producing and those consuming the metabolite and returns this information
    in a DataFrame.

    Args:
        model (cobra.Model, mandatory): The metabolic model to be analyzed.
        met_id (str, optional): The ID of the metabolite of interest.
        flux_vector (pd.Series, optional):

    Returns:
        pd.DataFrame: A DataFrame with two columns:
                    - `consuming_{met_id}`: Reactions consuming the specified metabolite in flux vector provided
                    - `producing_{met_id}`: Reactions producing the specified metabolite in flux vector provided

    Examples:

        >>> import cobra
        >>> model = cobra.io.read_sbml_model('e_coli_core.xml')
        >>> df = get_reactions_producing_met(model, 'met_id')
        >>> print(df)
            consuming_met_id producing_met_id
        0  reaction1          reaction2

    Notes:
        - The function assumes the model is properly optimized and that each reaction's flux can be accessed via the solution object.
        - Reactions with a flux less than 0 are considered to be consuming reactants and producing products in the reverse direction.
        - Reactions with a flux greater than 0 are considered to be producing reactants and consuming products in the forward direction.
        - This applies for the specific medium provided with the model
    """
    if model is None or not isinstance(model, cobra.Model):
        return TypeError("model needs to be a cobra.Model object.")

    if flux_vector is not None and not isinstance(flux_vector, pd.Series):
        return TypeError("flux_vector needs to be a pd.Series object.")

    if flux_vector is None:
        flux_vector = model.optimize()

    reactions_of_interest = model.metabolites.get_by_id(met_id).reactions

    reactions_producing_met = []
    reactions_consuming_met = []

    for reaction in reactions_of_interest:

        if flux_vector[reaction.id] < 0:

            for reactant in reaction.reactants:
                if met_id == reactant.id:
                    reactions_producing_met.append(reaction)
                    break

            for product in reaction.products:
                if met_id == product.id:
                    reactions_consuming_met.append(reaction)
                    break

        if flux_vector[reaction.id] > 0:

            for reactant in reaction.reactants:
                if met_id == reactant.id:
                    reactions_consuming_met.append(reaction)
                    break
            for product in reaction.products:
                if met_id == product.id:
                    reactions_producing_met.append(reaction)
                    break

    cons_df = pd.DataFrame({
        f'consuming {met_id}': reactions_consuming_met,
    })
    prod_df = pd.DataFrame({
        f'producing {met_id}': reactions_producing_met
    })
    return (cons_df, prod_df)


def track_metabolite_pathways(
    metabolite,
    model,
    skip_objective=True,
    _namespace=None,  # Literal["modelseed", "bigg"]
    _cofactors=None,
    _visited_metabolites=None,
    _inner_reactions=None,
    ):
    """
    Returns pathways that are initiated from a nutrient.

    """
    if _namespace is None:
        _namespace = extract_model_namespace(model)
    if _cofactors is None:
        _cofactors = MODELSEED_COFACTORS if _namespace == "modelseed" else BIGG_COFACTORS if _namespace == "bigg" else \
                    (_ for _ in ()).throw(ValueError("Invalid namespace. Only 'modelseed' or 'bigg' are allowed."))

    if skip_objective:
        for r in model.reactions:
            if r.objective_coefficient != 0 :
                model.remove_reactions([r])
        skip_objective = False

    _inner_reactions = _inner_reactions if _inner_reactions is not None else set()
    _visited_metabolites = _visited_metabolites if _visited_metabolites is not None else set()

    # Mark metabolite as visited
    _visited_metabolites.add(metabolite.id)

    reactions_with_the_met = model.metabolites.get_by_id(metabolite.id).reactions
    has_irreversible_reaction = any(not reaction.reversibility for reaction in reactions_with_the_met)

    print("\n\n============\n\n", "Metabolite under study:", metabolite.name, metabolite.id)
    print("Metabolite's total related reactions to be parsed:", reactions_with_the_met)

    # Check reactions involving this metabolite
    for reaction in reactions_with_the_met:

        # Add reaction to _inner_reactions if it's not already there
        if reaction not in _inner_reactions:
            _inner_reactions.add(reaction)

            print("\nReaction", reaction.id, "added in inner_reactions and is now under proc.")

            # ----
            # Traverse through metabolites in the reaction
            # ----

            # In case of a irreversible reaction..
            if not reaction.reversibility:
                if metabolite in reaction.reactants:
                    for in_metabolite in reaction.products:
                        # [NOTE] If it's an exchange reaction, then products are empty
                        print("1.Run recursive for IRREVERSIBLE met", in_metabolite.id, "from reaction in proc:", reaction.id)
                        _perform_recursive(in_metabolite, reaction, model, _visited_metabolites, _cofactors, _inner_reactions)
                else:
                    print("2.Metabolite:", metabolite.id, " not as reactant in:", reaction.build_reaction_string(),
                          "from reaction in proc:", reaction.id, ". Not new recursive.")

            # In case of a reversible reaction..
            else:
                # where the metabolite under study either is found only in this reaction or there is at least 1 reaction in those it participates that is reversible
                if len(reactions_with_the_met) == 1 or has_irreversible_reaction is False:

                    # then, parse through the reaction's products if the met under study is among its reactants
                    if metabolite in reaction.reactants:
                        for in_metabolite in reaction.products:
                            print("3.New recursive for REVERSIBLE reaction in proc", reaction.id, "in products, for met:", in_metabolite.id)
                            _perform_recursive(in_metabolite, reaction, model, _visited_metabolites, _cofactors, _inner_reactions)
                    # or the other way around
                    else:
                        for in_metabolite in reaction.reactants:
                            print("4.New recursive for REVERSIBLE reaction in proc", reaction.id, "in reactants, for met:", in_metabolite.id)
                            _perform_recursive(in_metabolite, reaction, model, _visited_metabolites, _cofactors, _inner_reactions)
                # where the met under study
                else:
                    print("5.Reaction", reaction.id, "is irreversible and there are more than 1 reactions with the met under study (", metabolite.id,
                          ") among which there is at least one reversible reaction. Thus, we skip further parsing for reac.")

    return(_inner_reactions, _visited_metabolites)




"""
Internal routines
"""
def _convert_list_to_binary(lst):
    """
    Converts an ndarray to a binary decision based on the presence of at least one nonzero/True value.

    This function checks if the input `lst` is a NumPy ndarray. If so, it evaluates whether there is
    at least one nonzero (or True) element in the array. If such an element exists, the function returns 1;
    otherwise, it returns 0. If `lst` is not a NumPy ndarray, it returns the original `lst` unchanged.

    Parameters:
    lst (any): The input that might be a NumPy ndarray.

    Returns:
    int or any: Returns 1 if `lst` is a NumPy ndarray with at least one nonzero/True value;
                returns 0 if `lst` is a NumPy ndarray with all elements zero/False;
                returns the original `lst` if it is not a NumPy ndarray.
    """
    if isinstance(lst, np.ndarray):
        if any(lst):
            return 1
        else:
            return 0
    else:
        return lst


def _convert_single_element_set(cell):
    """
    Converts a set with a single element to that element.

    This function checks if the input `cell` is a set containing exactly one element.
    If so, it converts the set to that single element. If `cell` is not a set with
    a single element, it returns the original `cell` unchanged.

    Parameters:
    cell (any): The input that might be a set with a single element.

    Returns:
    any: The single element if `cell` is a set with one element; otherwise, the original `cell`.
    """
    if isinstance(cell, set) and len(cell) == 1:
        return next(iter(cell))  # Convert set to string
    return cell


def _perform_recursive(in_met, reaction, model, visited_metabolites, cofactors, inner_reactions):
    # Check if metabolite is not in BIGG_COFACTORS and not visited
    if (
        in_met.id not in cofactors and
        in_met.id not in visited_metabolites
    ):
        print(">>", in_met.id, "from", reaction.build_reaction_string())
        # Recursively find inner reactions for this metabolite
        track_metabolite_pathways(
            metabolite=in_met,
            model=model,
            _cofactors = cofactors,
            _visited_metabolites=visited_metabolites,
            _inner_reactions=inner_reactions
        )




# %% Mapping metabolite and reaction terms

def search_names_to_ontology(terms: List[str], all_terms: Union[Dict, str], threshold=0.90, classic=False):
    """
    Map a list of terms to a reference ontology using string similarity.

    Parameters:
        terms (List[str]): List of terms to map.
        all_terms (Union[Dict, str]): Reference ontology as a dictionary or a path to a JSON file.
        threshold (float): Similarity threshold for fuzzy matching (used in classic mode).
        classic (bool): Whether to use `difflib` for fuzzy matching (default is False).

    Returns:
        Dict[str, str]: Mapping of input terms to the best matching reference term IDs.
    """
    # Load reference terms if provided as a file
    if isinstance(all_terms, str):
        try:
            with open(all_terms, "r") as f:
                all_terms = json.load(f)
        except Exception as e:
            raise ValueError("Provide a valid JSON file.") from e
    elif not isinstance(all_terms, Dict):
        raise TypeError("Provide reference terms as a dictionary or a valid JSON file.")

    # Extract aliases from reference terms
    alias = {k:v["alias"] for k,v in all_terms.items()}

    mmap = {}

    for term in terms:

        if classic:
            # Use difflib to find the best match
            best_hit_terms = []
            best_hit_ids = []
            for ref in alias:
                matches = get_close_matches(term, alias[ref], n=1, cutoff=threshold)
                if len(matches) > 0:
                    best_hit_terms.append(matches[0])
                    best_hit_ids.append(ref)

            if len(best_hit_terms) > 1:
                flattened_list = best_hit_terms
                best_hit_term = get_close_matches(term, flattened_list, n=1)[0]  # best_match[0]
                best_hit_id = best_hit_ids[flattened_list.index(best_hit_term)]

            elif len(best_hit_terms) == 1:
                best_hit_term, best_hit_id = best_hit_terms[0], best_hit_ids[0]

            else:
                print("no match for term:", term)
                best_hit_term, best_hit_id = None, None

        else:
            # Use Levenshtein distance to find the closest match
            best_hit_id = None
            best_hit_score = 100
            for ref in alias:
                for case in alias[ref]:
                    lev_score = levenshtein_distance(term, case)
                    if lev_score < best_hit_score:
                        best_hit_score = lev_score
                        best_hit_id = ref

        mmap[term] = best_hit_id

    return mmap


# %% Integrating growth data

def decide_compounds_being_produced_consumed(growth_data: pd.DataFrame):
    """
    Based on the first and last entry of a metabolite's concentration, decides if the metabolite is being
    produced or consumed.

    Usage:
        hplc_data = pd.read_excel(growth_data_file, sheet_name="HPLC")
        decide_compounds_being_produced_consumed(hplc_data)

    """
    decision = {}
    for _, met in growth_data.iterrows():
        t0, t1 = met.iloc[1], met.iloc[-1]
        diff = t1 - t0
        threshold = 0.1*t0

        if diff > threshold:
            decision[met["metabolite"]] = "produced"
        elif diff < threshold:
            decision[met["metabolite"]] = "consumed"
        else:
            decision[met["metabolite"]] = np.nan

    decision_df = pd.DataFrame.from_dict(decision, orient='index', columns=['Status']).reset_index()
    decision_df.columns = ['metabolite', 'sign']

    return decision_df


def add_exchange_reactions_for(model: cobra.Model, signs: pd.DataFrame):
    """
    signs:

    metabolite	sign	modelseed_id
	Glucose	consumed	cpd00027
    """

    with open(MSEED_COMPOUNDS, "r") as f:
        all_mseed_metabolites = json.load(f)

    for _, met in signs.iterrows():
        try:
            cpd = met["modelseed_id"]
        except KeyError:
            raise("Your signs data frame needs to include a `modelseed_id` column with the corresponding id.")

        if met["sign"] == np.nan:
            continue
        compound_name = met["metabolite"]
        try:
            c0 = cpd + "_c0"
            model.metabolites.get_by_id(c0)
        except:
            print(f"Metabolite {compound_name} ({c0}) is not present in the cytosol. it will be skipped")
            continue

        try:
            e0 = cpd + "_e0"
            model.metabolites.get_by_id(e0)
        except KeyError:
            print(f"We need to create the external metabolite for {e0}.")
            met_annotation = all_mseed_metabolites[cpd]
            e0_met = cobra.Metabolite(
                e0,
                formula = met_annotation["formula"],
                charge = met_annotation["charge"],
                compartment="e0",
                name=met_annotation["name"]
            )
            model.add_metabolites(e0_met)

        try:
            ex = "EX_" + cpd + "_e0"
            model.reactions.get_by_id(ex)
        except:
            print("we need to create an exchange reaction")
            print(e0)
            model.add_boundary(model.metabolites.get_by_id(e0), type="exchange")

    return model







