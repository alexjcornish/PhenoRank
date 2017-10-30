#!/usr/bin/env python -B

#==============================================================================
#title          :test_prince.py
#description    :unit tests of functions in prince.py
#author         :Alex J Cornish
#date_created   :20 November 2015
#version        :0.1.0
#usage          :
#python_version :2.7.9
#==============================================================================

import pandas as pd
import unittest
import phenorank.prince



class test_run_prince(unittest.TestCase):
  """
  Tests:
    1: outputs pandas dataframe
    2: pandas dataframe contains the correct number of dimensions
    3: output contains correct genes
    5: computes correct phenotype similarity scores (Y)
    6: computes correct propagated scores (SCORE)
    7: correctly masks associations
  """

  def test_run_prince_toy(self):
    # run function
    output = phenorank.prince.run_prince("D2", a=0.5, n=10, c=-15, factor=1, dir_data="test/data_prince/")

    # run tests
    self.assertTrue(isinstance(output, pd.DataFrame)) # 1
    self.assertItemsEqual(output.shape, (6, 3)) # 2

    # 3
    self.assertEquals(output["GENE"]["ENSG1"], "ENSG1")
    self.assertEquals(output["GENE"]["ENSG2"], "ENSG2")
    self.assertEquals(output["GENE"]["ENSG3"], "ENSG3")
    self.assertEquals(output["GENE"]["ENSG4"], "ENSG4")
    self.assertEquals(output["GENE"]["ENSG5"], "ENSG5")
    self.assertEquals(output["GENE"]["ENSG6"], "ENSG6")

    # 5
    self.assertAlmostEqual(output["Y"]["ENSG1"], 0.9969506, places=5)
    self.assertAlmostEqual(output["Y"]["ENSG2"], 0.0, places=5)
    self.assertAlmostEqual(output["Y"]["ENSG3"], 0.0, places=5)
    self.assertAlmostEqual(output["Y"]["ENSG4"], 0.0, places=5)
    self.assertAlmostEqual(output["Y"]["ENSG5"], 0.0, places=5)
    self.assertAlmostEqual(output["Y"]["ENSG6"], 0.1531325, places=5)

    # 6
    self.assertTrue(output["SCORE"]["ENSG1"] > output["SCORE"]["ENSG2"])
    self.assertEquals(output["SCORE"]["ENSG2"], output["SCORE"]["ENSG3"])
    self.assertEquals(output["SCORE"]["ENSG2"], output["SCORE"]["ENSG4"])
    self.assertAlmostEqual(output["SCORE"]["ENSG5"], 0.0, places=5)
    self.assertAlmostEqual(output["SCORE"]["ENSG6"], 0.1531325, places=5)

class test_compute_condition_scores(unittest.TestCase):
  """
  Tests:
    1: returns dictionary
    2: dictionary contains correct disease IDs
    3: if doid in pheno_sim then scores returned
    4: if doid not in pheno_sim then error produced
  """

  @classmethod
  def setUpClass(self):
    # setup
    self.pheno_sim = pd.DataFrame({"D1":[1, 0.1, 0.2], "D2":[0.1, 1, 0.3], "D3":[0.2, 0.3, 1]}, index=["D1", "D2", "D3"])

  def test_compute_condition_scores_toy(self):
    # run function
    output = phenorank.prince.compute_condition_scores("D2", self.pheno_sim, c=-15, factor=1)

    # run tests
    self.assertEquals(type(output), type({})) # 1
    self.assertEquals(len(output), 3) # 2
    self.assertAlmostEqual(output["D1"], 0.0004480129, places=5) # 2, 3
    self.assertAlmostEqual(output["D2"], 0.9969506, places=5) # 2, 3
    self.assertAlmostEqual(output["D3"], 0.008922289, places=5) # 2, 3

  def test_compute_condition_scores_error(self):

    with self.assertRaises(KeyError):
      phenorank.prince.compute_condition_scores("D0", self.pheno_sim, c=-15, factor=1)



class test_compute_gene_pheno_scores(unittest.TestCase):
  """
  Tests:
    1: returns pandas series
    2: pandas series contains correct genes
    3: returns correct score when one disease ID associated and it has a phenotypic score (ENSG1)
    4: returns correct score when one disease ID associated and it does not have a phenotypic score (ENSG2)
    5: returns correct score when multiple disease IDs associated and one one has a phenotypic score (ENSG3)
    6: returns correct score when multiple disease IDs assocaited and multiple have phenotypic scores (ENSG4)
  """

  def test_compute_gene_pheno_scores_toy(self):
    # setup
    gc = {"ENSG1":["D2"], "ENSG2":["D3"], "ENSG3":["D4", "D5"], "ENSG4":["D6", "D7"]}
    condition_score = {"D1":0.1, "D2":0.2, "D4":0.3, "D6":0.4, "D7":0.5}

    # run function
    output = phenorank.prince.compute_gene_pheno_scores(gc, condition_score)
    # run tests
    self.assertTrue(isinstance(output, pd.Series)) # 1
    self.assertItemsEqual(output.index.tolist(), gc.keys()) # 2
    self.assertEquals(output["ENSG1"], 0.2)
    self.assertEquals(output["ENSG2"], 0.0)
    self.assertEquals(output["ENSG3"], 0.3)
    self.assertEquals(output["ENSG4"], 0.5)



if __name__ == "__main__":
  unittest.main()
