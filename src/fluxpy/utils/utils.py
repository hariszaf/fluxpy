""" Functionalities to handle metabolic modeling analysis objects"""

import numpy as np
import pandas as pd
import cobra
from typing import Union, Optional
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


def find_connected_components(graph):
    # Define a function to find connected components using Depth-First Search (DFS)

    # A set to keep track of visited nodes
    visited = set()
    # A list to store the connected components
    components = []

    # Helper function to perform a recursive DFS
    def dfs(node, component):
        # Mark the current node as visited
        visited.add(node)
        # Add the node to the current component
        component.append(node)
        # Visit all unvisited neighbors of the current node
        for neighbor in graph[node]:
            if neighbor not in visited:
                dfs(neighbor, component)

    # Iterate through all nodes in the graph
    for node in graph:
        # If the node has not been visited, it's the start of a new connected component
        if node not in visited:
            # Create a new component (list) to hold connected nodes
            component = []
            # Perform DFS starting from this node
            dfs(node, component)
            # Add the completed component to the list of components
            components.append(component)

    # Return the list of connected components
    return components
        

def order_cofactors_by_abundance(sorted_counts_data, specific_reaction_cofactors):           
    """
    Function that takes the cofactors of a single reaction and sorts them 
    based on their overall abundance in the model

    Example input:
    
    sorted_counts_data = [('h_c', 35), ('h2o_c', 18), ('h_e', 17), ('atp_c', 13), 
                        ('adp_c', 12), ('nad_c', 12), ('nadh_c', 12), ('pi_c', 12)]
                        
    specific_reaction_cofactors = ['h_e', 'h_c', 'adp_c', 'nadh_c']
    
    Example output:
    
    [35, 17, 12, 12]
    """

    # Convert the list of tuples to a dictionary for fast lookup
    abundance_dict = dict(sorted_counts_data)
    
    # Sort the items based on their abundance in the dictionary
    sorted_cofactors = sorted(specific_reaction_cofactors, key=lambda x: abundance_dict.get(x, 0), reverse=True)
    sorted_abundances = [abundance_dict.get(item, 0) for item in sorted_cofactors]
    
    return sorted_abundances

        
def find_first_winning_sublist_with_losers(sublists):
    """
    Find the "winning" sublist and indices of non-winning sublists.
    
    Parameters:
        sublists (list of lists): Input list of sublists to compare.

    Returns:
        tuple: (winning_sublist, losers_indices)

    Example input:
    
    sublists = [
        [9, 3, 5, 7],
        [9, 3, 4, 6],
        [0, 5, 3, 3],
        [4, 2, 1, 6],
        [7, 3, 2, 0]
    ]
    
    Example output:
    
    [9, 3, 5, 7]
    
    """
    
    winning_sublist = max(sublists, key=lambda x: x)
    losers_indices = [i for i, sublist in enumerate(sublists) if sublist != winning_sublist]
    
    return winning_sublist, losers_indices