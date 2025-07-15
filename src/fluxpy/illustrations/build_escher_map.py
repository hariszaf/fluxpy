"""
Using cobramod we build Escher maps based on:

- a group (of the metabolic model)
- a list of reactions
- a KEGG pathway 

"""
import json
from pathlib import Path
from typing import Literal
from contextlib import suppress
# contextlib.suppress is a context manager in Python that allows you to ignore specific exceptions inside a with block — cleanly and explicitly.
from collections import defaultdict

import cobra
from cobramod.core import pathway as pt
from cobramod.visualization.converter import JsonDictionary
from ..constants import EXCLUDED_COMPOUNDS

# %%

def remove_cycles(hierarchy):
    """
    Takes a {parent: (child,…)} hierarchy and returns the same structure
    with cycles removed (based on DFS). All original keys are preserved,
    even if they end up with no children.
    """
    graph = {k: set(v) for k, v in hierarchy.items()}

    # Ensure all nodes are in the graph (even leaf-only ones)
    for children in graph.values():
        for c in children:
            graph.setdefault(c, set())

    visited, on_stack = set(), set()

    def dfs(u):
        visited.add(u)
        on_stack.add(u)

        for v in list(graph[u]):
            if v not in visited:
                dfs(v)
            elif v in on_stack:
                # Back-edge → would close a cycle
                graph[u].remove(v)

        on_stack.remove(u)

    for node in graph:
        if node not in visited:
            dfs(node)

    # ✅ Keep all nodes, even those with no children
    return {
        k: (None if not v else (sorted(v)[0] if len(v) == 1 else tuple(sorted(v))))
        for k, v in graph.items()
    }

# %%
from collections import defaultdict

def build_reaction_hierarchy(reaction_dict):
    def extract_id(m): return m.id if hasattr(m, 'id') else str(m)
    def get_ids(metabolite_list): return {extract_id(m) for m in metabolite_list}

    def get_effective_products(data):
        if data.get('reversibility'):
            return get_ids(data.get('f_products', [])) | get_ids(data.get('f_reactants', []))
        return get_ids(data.get('f_products', []))

    def get_effective_reactants(data):
        if data.get('reversibility'):
            return get_ids(data.get('f_reactants', [])) | get_ids(data.get('f_products', []))
        return get_ids(data.get('f_reactants', []))

    hierarchy = defaultdict(set)
    used_fallback = set()

    # Step 1: filtered fields (f_products/f_reactants), respecting reversibility
    for parent, data1 in reaction_dict.items():
        parent_outputs = get_effective_products(data1)
        for child, data2 in reaction_dict.items():
            if parent != child:
                child_inputs = get_effective_reactants(data2)
                if parent_outputs & child_inputs:
                    hierarchy[parent].add(child)

    # Step 2: add missing reactions using unfiltered products/reactants
    all_children = set(child for children in hierarchy.values() for child in children)
    missing_parents = all_children - set(hierarchy)

    for parent in missing_parents:
        data1 = reaction_dict[parent]
        parent_outputs = (
            get_ids(data1.get('products', [])) | get_ids(data1.get('reactants', []))
            if data1.get('reversibility')
            else get_ids(data1.get('products', []))
        )
        for child, data2 in reaction_dict.items():
            if parent != child:
                child_inputs = (
                    get_ids(data2.get('reactants', [])) | get_ids(data2.get('products', []))
                    if data2.get('reversibility')
                    else get_ids(data2.get('reactants', []))
                )
                if parent_outputs & child_inputs:
                    hierarchy[parent].add(child)
                    used_fallback.add(parent)

    return {k: tuple(sorted(v)) for k, v in hierarchy.items()}


# %%
def build_escher_map(
    model, 
    type: Literal["pathway", "graph"] = "pathway", 
    reaction_list=[], 
    pathway=None, KEGG_pathway_id=None, # TODO (Haris Zafeiropoulos, 2025-07-15): use pathway name or id as input, instead of list of reactions
    map_name="test_map",
    vertical=False,
    prev_gc = {}  # dev
):
    """
    Build Escher map from cobra model.

    Example:
        glycolysis = ['ALCD2x', 'ENO', 'FBA', 'FBP', 'GAPD', 'PFK', 'PGK', 'PGM', 'PPCK', 'PPS', 'PYK', 'TPI']
        ppp = ['FBA', 'FBP', 'GND', 'PFK', 'PGL', 'RPE', 'RPI', 'TKT1']
    """
    test_path = Path.cwd().joinpath("test_map.html")
    if len(reaction_list) > 0:

        try:
            members = {model.reactions.get_by_id(x) for x in reaction_list}
        except KeyError as e:
            print("Reaction(s) in list not part of the model.")

    members_parts = {}
    for reaction in members:

        members_parts[reaction.id] = {}

        members_parts[reaction.id]["reversibility"] = reaction.reversibility

        products  = [x for x in reaction.products]
        reactants = [x for x in reaction.reactants]

        members_parts[reaction.id]["products"]  = products
        members_parts[reaction.id]["reactants"] = reactants

        filtered_products  = [x for x in reaction.products if x.id not in EXCLUDED_COMPOUNDS]
        filtered_reactants = [x for x in reaction.reactants if x.id not in EXCLUDED_COMPOUNDS]

        members_parts[reaction.id]["f_products"]  = filtered_products
        members_parts[reaction.id]["f_reactants"] = filtered_reactants

    if len(prev_gc) > 0:
        child_graph = prev_gc
    else:
        child_graph = remove_cycles(build_reaction_hierarchy(members_parts))

    if type == "pathway":
        cobramod_obj = pt.Pathway(id="test_group", members=members)

    elif type == "graph":
        cobramod_obj = JsonDictionary()
        cobramod_obj.reaction_strings = {x.id: x.build_reaction_string() for x in members}

    with suppress(FileNotFoundError):
        test_path.unlink()

    cobramod_obj.graph = child_graph

    if type == "pathway":
        builder = cobramod_obj.visualize(vis = "escher-custom")
    elif type == "graph":
        builder = cobramod_obj.visualize(
            filepath           = test_path,
            vertical           = vertical,
            custom_integration = True
        )

    data = json.loads(builder.map_json)

    for _, r in data[1]["reactions"].items():
        for rxn in members:
            if r["name"] == rxn.id:
                r["reversibility"] = rxn.reversibility

    outfile = map_name + ".json"
    with open(outfile, "w") as f:
        json.dump(data, f)

    return child_graph, data, members


# %%
if __name__ == "__main__":

    # this is a test 
    model = cobra.io.read_sbml_model()
    glycolysis = ['ALCD2x', 'ENO', 'FBA', 'FBP', 'GAPD', 'PFK', 'PGK', 'PGM', 'PPCK', 'PPS', 'PYK', 'TPI']
    gc, data, members = build_escher_map(model, reaction_list=glycolysis, map_name="test", type="graph")

