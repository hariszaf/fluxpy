import os
import json
import pandas as pd
import unittest
import cobra
from fluxpy.utils.model import *

tests_dir = os.path.dirname(__file__)
root_dir = os.path.dirname(os.path.dirname(tests_dir))
models_dir = os.path.join(root_dir, "ext_data/models")
ecoli = os.path.join(models_dir, "e_coli_core.xml")
ba = os.path.join(models_dir, "B_animalis_ModelSEED.xml")
ec_model = cobra.io.read_sbml_model(ecoli)
ba_model = cobra.io.read_sbml_model(ba)

class TestUtils(unittest.TestCase):

    def test_NameIdConverter(self):
        sol = ba_model.optimize()
        sol_fluxes = NameIdConverter(
            model=ba_model,
            to_convert=sol.fluxes,
            convert="reactionIds",
            to="reactionNames"
        )
        f = sol_fluxes.run_convert()
        self.assertTrue("pyruvate" in list(f.index)[0])

    def test_sync_with_medium(self):
        edit_model = ba_model.copy()
        with open(os.path.join(root_dir, "ext_data/media/test_medium.json"), "r") as f:
            new_medium = json.load(f)
        edit_model, suggested_medium = sync_with_medium(model=edit_model, medium=new_medium)
        self.assertTrue("EX_cpd00122_e0" in edit_model.exchanges)
        self.assertTrue("EX_cpd37276_e0" not in suggested_medium)
        
    def test_cofactor_specificity(self):
        
        self.assertTrue(ec_model.reactions.FORt2.bounds == (0.0, 1000.0))
         
        cofactor_specificity = CofactorSpecificity(ec_model)
        cofactor_specificity.filter_reactions()
        self.assertTrue(ec_model.reactions.FORt2.bounds == (0, 0))

if __name__ == "__main__":
    unittest.main()
