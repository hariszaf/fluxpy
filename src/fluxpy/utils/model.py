""" Utilities to parse, convert and edit cobra models """

import cobra
import numpy as np
import pandas as pd
from typing import List, Dict, Union
from enum import Enum
from mergem import merge
from ..constants import *
from .utils import _convert_list_to_binary, _convert_single_element_set

# %% Inner functions
all_namespaces = ["chebi", "metacyc", "kegg", "reactome", "metanetx", "hmdb", "biocyc", "bigg", "seed", "sabiork", "rhea"]

class _Conversions(Enum):
    REACTION_NAMES = "reactionNames"
    REACTION_IDS = "reactionIds"
    METABOLITE_NAMES = "metaboliteNames"
    METABOLITE_IDS = "metaboliteIds"

# %% Util classes

class Models:
    """
    Load a list of sbml models with cobra as attributes of a class.
    """
    def __init__(self, model_list: List[str], model_names: List[str] = None):
        if model_names is None:
            model_names =  ["".join(model.split("/")[-1].split(".")[:-1]) for model in model_list]
        for model_name, model_file in zip(model_names, model_list):
            setattr(self, model_name, cobra.io.read_sbml_model(model_file))


class NameIdConverter:
    """
    A class to switch indexes of a dataframe from ids to names and vice versa for metabolites and reactions of a model.

    """
    def __init__(self, model, to_convert: Union[pd.DataFrame, Dict, List, cobra.Solution], convert: str, to: str):
        """
        Args:
            model (cobra.Model, mandatory):
            to_convert (pd.DataFrame, mandatory): a pandas dataframe to which the conversion will be applied to
            convert (str, mandatory): current namespace of the df; values from ['reactionNames', reactionIds, 'metaboliteNames', 'metaboliteIds']
            to (str, mandatory): namespace to switch to; values as in `convert`

        """

        try:
            convert, to = _Conversions(convert), _Conversions(to)
        except:
            raise ValueError(f"Both convert and to need to be among {[x.value for x in _Conversions.mro()]}")

        met_map_df = pd.DataFrame([(met.id, met.name) for met in model.metabolites], columns=["id", "name"])
        react_map_df = pd.DataFrame([(react.id, react.name) for react in model.reactions], columns=["id", "name"])
        self.met_map_df = met_map_df
        self.react_map_df = react_map_df
        self.convert = convert
        self.to = to
        self.model = model
        if isinstance(to_convert, pd.DataFrame):
            self.to_convert = to_convert.copy()
        elif isinstance(to_convert, cobra.core.solution.Solution):
            self.to_convert = to_convert.to_frame().copy()
        elif isinstance(to_convert, Dict):
            df = pd.DataFrame(list(to_convert.items()), columns=[to, 'Flux'])  # actually here is convert but it will be replaced
            df.set_index(to, inplace=True)
            self.to_convert = df.copy()
        elif isinstance(to_convert, List):
            df = pd.DataFrame(list(model.medium), columns=[to])
            df.set_index(to, inplace=True)
            self.to_convert = df.copy()

        self.isIndex = isinstance(to_convert, pd.Index)

    def run_convert(self):

        convert_cases = {
            (self.convert.REACTION_NAMES.value,   self.to.REACTION_IDS.value): self._react_names_to_react_ids,
            (self.convert.REACTION_IDS.value,     self.to.REACTION_NAMES.value): self._react_ids_to_react_names,
            (self.convert.METABOLITE_NAMES.value, self.to.METABOLITE_IDS.value): self._met_names_to_met_ids,
            (self.convert.METABOLITE_IDS.value,   self.to.METABOLITE_NAMES.value): self._met_ids_to_met_names
        }
        # Retrieve the method associated with the specific combination of convert and to values
        convert_res = convert_cases.get((self.convert.value, self.to.value))
        if convert_res is not None:
            # Execute the corresponding function
            s = convert_res()
        else:
            # Handle the case when the combination doesn't exist
            raise ValueError("Invalid combination of convert and to values")

        self.to_convert = s
        return s


    def _met_ids_to_met_names(self):

        self.met_map_df.set_index('id', inplace=True)
        id_to_name_mapping = self.met_map_df['name'].to_dict()
        self.to_convert.index = self.to_convert.index.map(id_to_name_mapping)
        if self.isIndex:
            self.to_convert.columns = ["id"]
            self.to_convert = pd.DataFrame(self.to_convert.index, index = self.to_convert["id"])
            self.to_convert.columns = ["name"]
        return self.to_convert

    def _met_names_to_met_ids(self):

        self.met_map_df.set_index('name', inplace=True)
        id_to_name_mapping = self.met_map_df['id'].to_dict()
        self.to_convert.index = self.to_convert.index.map(id_to_name_mapping)
        if self.isIndex:
            self.to_convert.columns = ["id"]
        return self.to_convert


    def _react_ids_to_react_names(self):

        self.react_map_df.set_index('id', inplace=True)
        id_to_name_mapping = self.react_map_df['name'].to_dict()
        self.to_convert.index = self.to_convert.index.map(id_to_name_mapping)
        if self.isIndex:
            self.to_convert.columns = ["name"]
        return self.to_convert


    def _react_names_to_react_ids(self):

        self.react_map_df.set_index('name', inplace=True)
        id_to_name_mapping = self.react_map_df['id'].to_dict()
        self.to_convert.index = self.to_convert.index.map(id_to_name_mapping)
        if self.isIndex:
            self.to_convert.columns = ["id"]
        return self.to_convert


class CompareModels:
    """
    Using mergem <https://mergem.readthedocs.io/> to compare a pair of models and returns a dataframe with
    reactions unique in model1 and model2 and those they share.

    Key arguments:\n
    model_list -- list of paths to model files to compare \n
    trans_to_db -- mergem allows to map metabolite and reaction ids to a series of namespaces.
    """

    def __init__(self, model_list: List[str], trans_to_db=None):

        self.number_of_models = len(model_list)
        self.model_names = ["".join(model.split("/")[-1].split(".")[:-1]) for model in model_list]
        self.models = Models(model_list)
        if trans_to_db is not None and trans_to_db not in all_namespaces:
            raise ValueError(f"namespace needs to be among {all_namespaces}")
        self.namespace = trans_to_db

        # Run mergem
        self.res = merge(model_list, set_objective='merge', trans_to_db=self.namespace)

        # Create presence-absence dataframes
        res = self.res.copy()
        for lst in res["met_sources"]:
            for i in range(self.number_of_models):
                if i not in res["met_sources"][lst]:
                    res["met_sources"][lst][i] = np.nan
        for lst in res["reac_sources"]:
            for i in range(self.number_of_models):
                if i not in res["reac_sources"][lst]:
                    res["reac_sources"][lst][i] = np.nan

        metabolites_df = pd.DataFrame.from_dict(res["met_sources"], orient="index")
        metabolites_df = metabolites_df.applymap(lambda x: np.where(pd.notnull(x), 1, 0))
        self.metabolites_df = metabolites_df.applymap(lambda x: _convert_list_to_binary(x))
        self.metabolites_df.columns = self.model_names

        reactions_df = pd.DataFrame.from_dict(self.res["reac_sources"], orient="index")
        reactions_df = reactions_df.applymap(lambda x: np.where(pd.notnull(x), 1, 0))
        self.reactions_df = reactions_df.applymap(lambda x: _convert_list_to_binary(x))
        self.reactions_df.columns = self.model_names


    def unique_metabolites(self, model_name):
        """
        Returns:
        unique_mets -- a pd.Index with the metabolites that are only present in model_name and not in any other of those in model_list
        """
        filtered_df = self.metabolites_df[self.metabolites_df[model_name] == 1]
        unique_mets = filtered_df.index[filtered_df.drop(columns=[model_name]).sum(axis=1) == 0]
        return unique_mets

    def unique_reactions(self, model_name):
        """
        Key arguments:
        model_name -- the name (column name) of the model for which unique reactions will be found

        Returns:
        model_unique_reactions -- a pd.Index with the reactions that are only present in the model_name and not in any other of those in model_list
        """
        filtered_df = self.reactions_df[self.reactions_df[model_name] == 1]
        unique_reacts = filtered_df.index[filtered_df.drop(columns=[model_name]).sum(axis=1) == 0]
        return unique_reacts


    def compare_model_pair(self, base_model=None, model_to_compare=None):
        """
        Compare pairwise models

        Key arguments:
        base_model -- name of the model as in dataframe, i.e. filename without the extension
        model_to_compare -- name of the model as in dataframe, i.e. filename without the extension

        Returns:
        r1 -- reactions only present in base model
        r2 -- reactions only present in model_to_compare
        m1 -- metabolites only present in base_model
        m2 -- metabolites only present in model_to_compare
        """

        if base_model is None and model_to_compare is None and self.number_of_models == 2:
            base_model, model_to_compare = self.model_names[0], self.model_names[1]

        reactions_m1 = self.reactions_df[self.reactions_df[base_model] == 1].index
        reactions_m2 = self.reactions_df[self.reactions_df[model_to_compare] == 1].index
        r1 = reactions_m1.difference(reactions_m2)
        r2 = reactions_m2.difference(reactions_m1)

        metabolites_m1 = self.reactions_df[self.metabolites_df[base_model] == 1].index
        metabolites_m2 = self.reactions_df[self.metabolites_df[model_to_compare] == 1].index
        m1 = metabolites_m1.difference(metabolites_m2)
        m2 = metabolites_m2.difference(metabolites_m1)
        return r1, r2, m1, m2


    def shared_metabolites(self, models_subset=None):
        """
        Key arguments:
        models_subset -- list of models to get shared metabolites and reactions

        Returns:
        shared_reactions -- a pd.Index with rearctions present among all models of models_subset
        """
        if models_subset is None:
            models_subset = self.model_names
        df = self.metabolites_df[models_subset]
        shared_metabolites = df[df.eq(1).all(axis=1)]

        return shared_metabolites.index


    def shared_reactions(self, models_subset=None):
        """
        Key arguments:
        models_subset -- list of models to get shared metabolites and reactions

        Returns:
        shared_reactions -- a pd.Index with rearctions present among all models of models_subset
        """
        if models_subset is None:
            models_subset = self.model_names
        df = self.reactions_df[models_subset]
        shared_reactions = df[df.eq(1).all(axis=1)]

        return shared_reactions.index


    # def compare_sublist(self, sublist):

    #     return 1



        # self.shared_metabolites =
        # self.shared_reactions =
        # self.unique_reactions_model1 =
        # self.unique_reactions_model2 =
        # self.unique_metabolites_model1 =
        # self.unique_metabolites_model2 =


# %% Util functions
def _check_if_modelseed_model(model: cobra.Model):
    return model.metabolites[0].id.startswith("cpd")


def map_namespace_to_ModelSEEDDatabase(seed_compounds: List[str],
                                       path_to_modelSEEDDatabase: str,
                                       annotations_to_keep=['BiGG', 'BiGG1']):
    """
    This function makes use of the ModelSEEDpy and the ModelSEEDDatabase libraries to map a list of ModelSEED compounds
    to their annotations returning a pandas dataframe with those

    seed_compouds -- a list with seed compound ids
    path_to_modelSEEDDatabase -- path to the modelSEEDDatabase. To get it: git clone https://github.com/ModelSEED/ModelSEEDDatabase.git
    """
    from modelseedpy.biochem import from_local
    import gc  ## gc (garbage collection) module to force garbage collection.

    if path_to_modelSEEDDatabase is None:
        raise ValueError("path_to_modelSEEDdatabase is required.")
    try:
        modelseed_local = from_local(path_to_modelSEEDDatabase)
    except:
        raise ValueError("The modelSEEDDatabase provided was not correct. Please give correct path. For example: /Users/workspace/data/ModelSEEDDatabase/")
    dic = {}
    for medium_seed_met in seed_compounds:
        for met in modelseed_local.compounds:
            if met.id == medium_seed_met:
                found_annotations = {bigg: (bigg in met.annotation) for bigg in annotations_to_keep}
                if any(found_annotations.values()):
                    tmp = {annotation: met.annotation[annotation] if found else None for annotation, found in found_annotations.items()}
                    dic[medium_seed_met] = tmp
                else:
                    dic[medium_seed_met] = None
                break
    for key, value in dic.items():
        if value is None:
            dic[key] = {}
        for subkey in annotations_to_keep:
            dic[key].setdefault(subkey, None)
    t = pd.DataFrame.from_dict(dic, orient="index")
    df = t.applymap(_convert_single_element_set)
    del modelseed_local
    gc.collect()
    return df


def map_to_namespace(compounds, from_namespace="seed", to_namespace="bigg"):
    """
    compounds -- a list
    from_namesapce --
    to_namespace --
    """
    seed2metanex = pd.read_json(SEED2MNX)
    bigg2metanex = pd.read_json(BIGG2MNX)
    if from_namespace == "seed":
        seed_metanex = seed2metanex.loc[compounds].to_dict()["metanex_id"]
        metanex_ids = list(set(seed_metanex.values()))
        bigg_metanex = bigg2metanex[bigg2metanex["metanex_id"].isin(metanex_ids)]
        bigg_metanex = bigg_metanex["metanex_id"].to_dict()
        mapped_dict = {
            seed_key: next((bigg_key for bigg_key, bigg_value in bigg_metanex.items() if bigg_value == seed_value), None)
            for seed_key, seed_value in seed_metanex.items()
        }
        data_for_df = {'seed_metanex': list(mapped_dict.keys()), 'bigg_metanex': list(mapped_dict.values())}
    elif from_namespace == "bigg":
        bigg_metanex = bigg2metanex.loc[compounds].to_dict()["metanex_id"]
        metanex_ids = list(set(seed_metanex.values()))
        seed_metanex = seed2metanex[seed2metanex["metanex_id"].isin(metanex_ids)]
        seed_metanex = seed_metanex["metanex_id"].to_dict()
        mapped_dict = {
            bigg_key: next((seed_key for seed_key, seed_value in seed_metanex.items() if seed_value == bigg_value), None)
            for bigg_key, bigg_value in bigg_metanex.items()
        }
        data_for_df = {'bigg_metanex': list(mapped_dict.keys()), 'seed_metanex': list(mapped_dict.values())}
    df = pd.DataFrame(data_for_df)
    return df


def sync_with_medium(model: cobra.Model, medium: Dict):
    """
    Adds metabolites and reactions required to allow a medium to be assigned to a cobra.Model

    Args:
        model (cobra.Model, mandatory): model to which the new medium will be assigned to and to which metabolites and reactions will be added to
        medium (Dict, mandatory): a dictionary with the medium to be used with the metabolites as keys and their boundaries as values

    Returns:
        A cobra.Model with added metabolites and reactions to enable medium to be used

    """
    if not isinstance(model, cobra.Model) or not isinstance(medium, Dict):
        return TypeError("")
    try:
        model.medium = medium
        return "The cobra model is already capable of supporting this medium."
    except:
        pass
    # Check if model is using ModelSEED ontology
    if not _check_if_modelseed_model(model):
        return "Currently, cobra.Model needs to use ModelSEED ontology."
    complete_model = cobra.io.read_sbml_model(COMPLETE_MODEL)
    rxns_to_add = set()
    user_model_exchange_ids = [x.id for x in model.exchanges]
    suggested_medium = {}
    for ex_react in medium:
        edit_ex_react = ex_react
        if edit_ex_react.startswith("cpd"):
            edit_ex_react = "EX_" + edit_ex_react
        if not edit_ex_react.endswith("e0"):
            edit_ex_react = edit_ex_react + "_e0"
        if edit_ex_react not in user_model_exchange_ids:
            try:
                rxn_to_add = complete_model.reactions.get_by_id(edit_ex_react)
                rxns_to_add.add(rxn_to_add)
                suggested_medium[edit_ex_react] = medium[ex_react]
                # logging.info('Reaction added: ', edit_ex_react)
                print('Reaction added: ', edit_ex_react)
            except:
                print("Reaction", ex_react,
                      "is not part of the complete model and most likely has not related reactions or it is obsolete. \
                      We suggest to remove it from your medium."
                )
                pass
        else:
            suggested_medium[edit_ex_react] = medium[ex_react]

    rxns_to_add = list(rxns_to_add)
    model.add_reactions(rxns_to_add)

    return (model, suggested_medium)


def objective_reaction_name(model: cobra.Model):
    """
    Returns the name and the id of the objective function
    """
    obj_index = [r.objective_coefficient for r in model.reactions].index(True)
    return model.reactions[obj_index].id, model.reactions[obj_index].name



def fix_constraints_based_on_growth(model: cobra.Model, signs: pd.DataFrame):
    """
    Gets a model and a series of data frames or dictionaries as input to constrain the model based on growth data.


    [NOTE] check how different the solutions go if you fix an exchange reaction to something compared to if you just specify the boundary sign

    sigs:
        The `metabolite` column is optional. Also, modelseed_id column may refer to specific exchange reactions, e.g. EX_cpd00027_e0 or include the `_e0` suffix in case of metabolites.
        metabolite	status	modelseed_id
        Glucose	consumed	cpd00027

    """

    for _, row in signs.iterrows():

        try:
            cpd = row["modelseed_id"]
        except KeyError:
            raise("Your signs data frame needs to include a `modelseed_id` column with the corresponding id.")

        if cpd.endswith("_e"):
            cpd = cpd.replace("_e", "_e0")
        if "_e0" not in cpd:
            cpd += "_e0"
        if "EX" not in cpd:
            rxn = "EX_" + cpd
        else:
            rxn = cpd

        try:
            model.metabolites.get_by_id(cpd)
        except:
            print(f"Metabolite {cpd} is not present in the model. It will be skipped.")
            continue

        model_rxn = model.reactions.get_by_id(rxn)

        if row["sign"] == "consumes":
            model_rxn.bounds = model_rxn.lower_bound, 0

        elif row["sign"] == "produces":
            model_rxn.bounds = 0, model_rxn.upper_bound

        else:
            continue

        return model




"""
loopless_S1 = cobra.flux_analysis.loopless_solution(cobra_model_banimalis, banimalis_10opt_s1)
# Calculate the absolute differences
abs_diff = np.abs(dingo_sample_loopless.fluxes - s1)

# Find the indices where the difference is greater than 0.001
indices = np.where(abs_diff > 0.001)

# Get the corresponding values from both arrays
differences = {
    "dingo_sample_loopless.fluxes": dingo_sample_loopless.fluxes[indices],
    "s1": s1[indices],
    "absolute_difference": abs_diff[indices]
}
"""