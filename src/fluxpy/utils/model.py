""" Utilities to parse, convert and edit cobra models """

import cobra
import numpy as np
import pandas as pd
from typing import List, Dict, Union, Literal
from enum import Enum
from mergem import merge
from ..constants import *
from .utils import _convert_list_to_binary, _convert_single_element_set, order_cofactors_by_abundance, find_connected_components, find_first_winning_sublist_with_losers
from cobra.io import load_json_model
from collections import defaultdict
from collections import Counter


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
    if not _check_if_modelseed_model:
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


"""
Inner routines
"""
def _check_if_modelseed_model(model: cobra.Model):
    return model.metabolites[0].id.startswith("cpd")



class CofactorSpecificity:
    
    """
    Class is used to identify cases where the same biochemical conversion is carried out by different cofactors in the model, 
    while in reality the organism only uses one of the cofactors. If the cofactor specificity information is available,
    the reaction with non-specific cofactor can be removed or turned off.
    
    This is how to initialize the class:
    cofactor_specificity = CofactorSpecificity(cobra_model)
    """
    
    def __init__(self, model):
        
        self.model = model


    def find_reactants_products_cofactors(self):
        
        self.reactions_ids =  [ reaction.id for reaction in self.model.reactions ]
        
        self.reactants_list_all_reactions = []
        self.products_list_all_reactions = []
        self.cofactors_list_all_reactions = []
        self.reversibility_list_all_reactions = []
        
        for reaction in self.reactions_ids:
            reactants_list_single_reaction = []
            products_list_single_reaction = []
            cofactors_list_single_reaction = []

            reaction_information = self.model.reactions.get_by_id(reaction)
            
            reactants = reaction_information.reactants
            products = reaction_information.products
            
            reversibility = reaction_information.reversibility
            self.reversibility_list_all_reactions.append(reversibility)
            
            for reactant in reactants:
                reactant = str(reactant)
                if reactant in BIGG_COFACTORS or reactant in BIGG_BUILDING_BLOCLS:
                    cofactors_list_single_reaction.append(reactant)
                else:
                    reactants_list_single_reaction.append(reactant)
                    
            for product in products:
                product = str(product)
                if product in BIGG_COFACTORS or product in BIGG_BUILDING_BLOCLS:
                    cofactors_list_single_reaction.append(product)
                else:
                    products_list_single_reaction.append(product)         
                    
            #print(reaction_information, reactants_list_single_reaction, products_list_single_reaction, cofactors_list_single_reaction, end ="\n")
            self.reactants_list_all_reactions.append(reactants_list_single_reaction)
            self.products_list_all_reactions.append(products_list_single_reaction)
            self.cofactors_list_all_reactions.append(cofactors_list_single_reaction)


    def find_reactions_combinations(self):
        
        self.find_reactants_products_cofactors()
        
        combinations = []
        
        for i in range(len(self.reactions_ids)):
            for j in range(len(self.reactions_ids)):
                
                # avoid comparison of the same reaction
                if self.reactions_ids[i] != self.reactions_ids[j]:

                    # avoid comparison with empty lists
                    if len(self.reactants_list_all_reactions[i]) > 0 and len(self.reactants_list_all_reactions[j]) > 0 and \
                    len(self.products_list_all_reactions[i]) > 0 and len(self.products_list_all_reactions[j]) > 0:
                    
                        # boolean to check if pairwise reactants/products are the same
                        identical_reactants = (set(self.reactants_list_all_reactions[i]) == set(self.reactants_list_all_reactions[j]))
                        identical_products = (set(self.products_list_all_reactions[i]) == set(self.products_list_all_reactions[j]))
                        identical_cofactors = (set(self.cofactors_list_all_reactions[i]) == set(self.cofactors_list_all_reactions[j]))

                        if identical_reactants == True and identical_products == True and identical_cofactors == False:
                            # finish here and keep combination
                            combinations.append((self.reactions_ids[i], self.reactions_ids[j]))
                        else:
                            if self.reversibility_list_all_reactions[i] == True or self.reversibility_list_all_reactions[j] == True:
                                # switch the reversibility of reaction i and compare updated reactants-products
                                identical_reactants = (set(self.reactants_list_all_reactions[i]) == set(self.products_list_all_reactions[j]))
                                identical_products = (set(self.products_list_all_reactions[i]) == set(self.reactants_list_all_reactions[j]))

                                if identical_reactants == True and identical_products == True and identical_cofactors == False:
                                    # finish here and keep combination
                                    combinations.append((self.reactions_ids[i], self.reactions_ids[j]))

        self.unique_combinations = {tuple(sorted(pair)) for pair in combinations}


    def merge_connected_combinations(self):
        
        """
        This function takes as input the pairwise combinations of reactions from above and
        creates groups of reactions that are candidates for cofactor specificity. It uses a graph-based approach
        to find connected components in the graph of reactions.
        """
        
        self.find_reactions_combinations()
        
        # We'll use a defaultdict where each key is a node, and its value is a set of neighboring nodes
        graph = defaultdict(set)

        # Build the adjacency list by iterating through each pair in the input
        for x, y in self.unique_combinations:
            # Add y as a neighbor of x
            graph[x].add(y)
            # Add x as a neighbor of y (because the graph is undirected)
            graph[y].add(x)

        # At this point, for an input of this format: [('A', 'B'), ('A', 'C'), ('B', 'C'), ('D', 'K')]
        # the `graph` looks like:
        # {'A': {'B', 'C'}, 'B': {'A', 'C'}, 'C': {'A', 'B'}, 'D': {'K'}, 'K': {'D'}}

        self.groups = find_connected_components(graph)

        # At this point, `self.groups` will contain:
        # [['A', 'B', 'C'], ['D', 'K']]
        

    def cofactors_abundance(self):
        
        self.merge_connected_combinations()
                
        # Count occurrences of each cofactor
        flat_list = [item for sublist in self.cofactors_list_all_reactions for item in sublist]
        element_counts = Counter(flat_list)  
        self.sorted_counts = element_counts.most_common()
        
        
    def remove_reactions(self):
        
        for knockout_index in self.knockout_indices_all:
            knockout_reaction = self.reactions_ids[knockout_index]
            self.model.reactions.get_by_id(knockout_reaction).lower_bound = 0
            self.model.reactions.get_by_id(knockout_reaction).upper_bound = 0
                
        
    def filter_reactions(self):
        
        self.cofactors_abundance()
        self.knockout_indices_all = []
        
        for group in self.groups:
            keep_indices_group = []
            cofactors_count_min = float('+inf')
                        
            # find minimum number of cofactors across reactions of a single group
            for reaction in group:
                reaction_index = self.reactions_ids.index(reaction)   
                reaction_cofactors = self.cofactors_list_all_reactions[reaction_index]
                cofactors_count = len(reaction_cofactors)
                
                if cofactors_count < cofactors_count_min:
                    cofactors_count_min = cofactors_count

                #print(reaction, reaction_cofactors, cofactors_count, cofactors_count_min)
                
            # find which reactions match the minimum number of cofactors
            for reaction in group:
                # find the index from the reaction list (general)
                reaction_index = self.reactions_ids.index(reaction)                   
                reaction_cofactors = self.cofactors_list_all_reactions[reaction_index]
                cofactors_count = len(reaction_cofactors)
   
                # equal to minimum value
                if cofactors_count == cofactors_count_min:
                    keep_indices_group.append(reaction_index)
                    
                # greater than minimum value ==> knockout
                else:
                    self.knockout_indices_all.append(reaction_index)

            
            # if 2 or more reactions match the minimum number of cofactors
            if len(keep_indices_group) > 1:
                # find which reaction has the most abundant cofactor (general abundance from the model) 
                # if the most abundant cofactor appears in all reactions, then check the next most abundant cofactor

                # apply the 'order_cofactors_by_abundance' function to reactions in the 'keep_indices' list
                # of the current group. Add all sorted abundance of cofactors for each reactions in a list
                # then call the 'find_first_winning_sublist_with_losers' and add indices of the loser sublists into a list
                                    
                sorted_abundance_all_reactions = [
                order_cofactors_by_abundance(self.sorted_counts, self.cofactors_list_all_reactions[idx])
                for idx in keep_indices_group ]
                
                _, losers_indices = find_first_winning_sublist_with_losers(sorted_abundance_all_reactions)
                
                for loser_index in losers_indices:
                    reaction_general_index = self.reactions_ids.index(group[loser_index])
                    self.knockout_indices_all.append(reaction_general_index)                            
            
        
        self.remove_reactions()
