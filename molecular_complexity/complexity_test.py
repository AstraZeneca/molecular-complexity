import unittest

from .complexity import molecular_complexity

TEST_CASES = [
    {
        "smiles": "CCOCC",
        "cm": 19.67970000576925,
        "cm_star": 6.272887502165801,
        "cse": 1.5219280948873621,
    },
    {
        "smiles": "CCCCC",
        "cm": 22.399498959955935,
        "cm_star": 6.8526094185572735,
        "cse": 1.5219280948873621,
    },
    {
        "smiles": "CCCCCC",
        "cm": 27.399498959955935,
        "cm_star": 7.197303740040691,
        "cse": 1.584962500721156,
    },
    {
        "smiles": "CCCCCCC",
        "cm": 31.899498959955935,
        "cm_star": 7.404173141982172,
        "cse": 1.950212064914747,
    },
    {
        "smiles": "CCCCC(C)C",
        "cm": 32.597515671483464,
        "cm_star": 7.513437964818001,
        "cse": 2.521640636343318,
    },
    {
        "smiles": "CCCC(C)CC",
        "cm": 33.06197681726271,
        "cm_star": 7.655873311439077,
        "cse": 2.807354922057604,
    },
    {"smiles": "C1CC1", "cm": 13.5, "cm_star": 6.084962500721156, "cse": -0.0},
    {"smiles": "C1CCC1", "cm": 18.0, "cm_star": 6.5, "cse": -0.0},
    {"smiles": "C1CCCC1", "cm": 22.5, "cm_star": 6.821928094887363, "cse": -0.0},
    {"smiles": "C1CCCCC1", "cm": 27.0, "cm_star": 7.084962500721156, "cse": -0.0},
    {
        "smiles": "c1ccccc1",
        "cm": 23.06313713864835,
        "cm_star": 6.428818690495881,
        "cse": -0.0,
    },
    {
        "smiles": "C1CCOCC1",
        "cm": 26.50977500432694,
        "cm_star": 7.056813643492505,
        "cse": 1.9182958340544893,
    },
    {
        "smiles": "c1ccncc1",
        "cm": 23.331568569324176,
        "cm_star": 6.527502710047034,
        "cse": 1.9182958340544893,
    },
    {
        "smiles": "C1CCNCC1",
        "cm": 28.27270096091705,
        "cm_star": 7.327304070946181,
        "cse": 1.9182958340544893,
    },
    {
        "smiles": "C1=CCNCC1",
        "cm": 29.068190250958274,
        "cm_star": 7.441773908420807,
        "cse": 2.584962500721156,
    },
    {
        "smiles": "C1NCC2CC12",
        "cm": 31.131841911895073,
        "cm_star": 7.849225522677325,
        "cse": 1.9182958340544893,
    },
    {
        "smiles": "CC1CCCN1",
        "cm": 30.58902752925798,
        "cm_star": 7.743244843716409,
        "cse": 2.584962500721156,
    },
    {
        "smiles": "c1ccccc1C(=O)N",
        "cm": 35.466395473423994,
        "cm_star": 7.272421335447814,
        "cse": 2.725480556997868,
    },
    {
        "smiles": "COc1ccc(C(C)=O)cc1",
        "cm": 42.76368520454086,
        "cm_star": 7.549097942245354,
        "cse": 3.095795255000934,
    },
    {
        "smiles": "CN(C)CCOC(c1ccccc1)c1ccccc1",
        "cm": 82.70587659328551,
        "cm_star": 8.658352930496191,
        "cse": 3.0900327766014803,
    },
    {
        "smiles": "CCCCCc1cc(c2c(c1)OC([C@H]3[C@H]2C=C(CC3)C)(C)C)O",
        "cm": 109.68643136962028,
        "cm_star": 9.43178673022628,
        "cse": 4.436605434317882,
    },
    {
        "smiles": "COc1ccc2nc([nH]c2c1)[S@](=O)Cc1ncc(C)c(OC)c1C",
        "cm": 99.09550093215142,
        "cm_star": 8.904435585569017,
        "cse": 4.501629167387822,
    },
    {
        "smiles": "CC#CC(=O)N1CCC[C@H]1c2nc(c3n2ccnc3N)c4ccc(cc4)C(=O)Nc5ccccn5",
        "cm": 148.40551799733018,
        "cm_star": 9.59548257488669,
        "cse": 4.900711588373538,
    },
    {
        "smiles": "Cc1ccccc1S(=O)(=O)NC(=O)c2cc(OC)c(cc2)Cc3cn(C)c4ccc(cc43)NC(=O)OC5CCCC5",
        "cm": 179.71828246252102,
        "cm_star": 9.922332688822255,
        "cse": 5.1624300533985705,
    },
    {
        "smiles": "O=C(N1[C@H](C#N)C[C@@H]2C[C@H]12)[C@@H](N)C35CC4CC(C3)CC(O)(C4)C5",
        "cm": 106.80503550362938,
        "cm_star": 9.517133268814485,
        "cse": 4.2626923908396215,
    },
    {
        "smiles": "O=C(NCc1ccc(C(=N\\O)\\N)cc1)[C@H]3N(C(=O)[C@H](NCC(=O)OCC)C2CCCCC2)CC3",
        "cm": 151.2083891818507,
        "cm_star": 9.84540386395404,
        "cse": 4.85216872360328,
    },
]


class ComplexityTest(unittest.TestCase):
    def test_complexity_calculation(self):
        """Test that we can calculate molecular complexity"""
        for test_case in TEST_CASES:
            with self.subTest(smiles=test_case["smiles"]):
                cm, cm_star, cse = molecular_complexity(test_case["smiles"])
                self.assertAlmostEqual(cm, test_case["cm"])
                self.assertAlmostEqual(cm_star, test_case["cm_star"])
                self.assertAlmostEqual(cse, test_case["cse"])


if __name__ == "__main__":
    res = unittest.main(verbosity=3, exit=False)
