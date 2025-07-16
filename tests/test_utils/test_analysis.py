import os
import unittest
from fluxpy.utils.analysis import producing_or_consuming_a_met, get_reactions_producing_a_met, \
    trace_path, find_shortest_path_in_reaction_list
import cobra


tests_dir = os.path.dirname(__file__)
root_dir = os.path.dirname(os.path.dirname(tests_dir))
models_dir = os.path.join(root_dir, "ext_data/models")
ecoli = os.path.join(models_dir, "e_coli_core.xml")


class TestUtils(unittest.TestCase):

    def test_producing_or_consuming_a_met(self):
        self.assertTrue(producing_or_consuming_a_met(model=ecoli, reaction_id="CS", metabolite_id="accoa_c") == "consuming")

    def test_get_reactions_producing_a_met(self):
        self.assertTrue(get_reactions_producing_a_met(model=ecoli, metabolite_id="accoa_c")[0].id == "PDH" )


    def test_trace_path(self):
        cofactors = ["h2o_c", "h_c", "pi_c", "nad_c", "atp_c", "adp_c", "nadph_c", "nadp_c"]
        trace_path(model=ecoli, start_reaction_id="", target_reaction_id="", ignore_mets=cofactors)

        self.assertTrue()

    def test_find_shortest_path_in_reaction_list(self):

        find_shortest_path_in_reaction_list()

        self.assertTrue()

