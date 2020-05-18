import unittest
import pandas
import descriptastorus.descriptors

class TestPandasDataFrame(unittest.TestCase):
    def test_pd_descriptors(self):
        smiles = ["C"*i for i in range(1,10)]
        df = pandas.DataFrame(dict(smiles=smiles))
        df2 = descriptastorus.descriptors.create_descriptors(
            df, mols_column_name="smiles", generator_names=["Morgan3Counts","RDKit2D"])
        g = descriptastorus.descriptors.MakeGenerator(["Morgan3Counts","RDKit2D"])
        self.assertEqual(len(df2), 9)
        for a,b in zip(df2.columns, g.GetColumns()):
            self.assertEqual(a,b)

        
if __name__ == "__main__":
    print("running")
    unittest.main()
