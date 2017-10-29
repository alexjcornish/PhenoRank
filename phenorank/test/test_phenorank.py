#!/usr/bin/env python -B

#==============================================================================
#title          :test_phenorank.py
#description    :unit test functions in phenorank.py
#author         :Alex J Cornish
#date_created   :21 May 2015
#version        :0.1.0
#usage          :
#python_version :2.7.9
#==============================================================================

import logging
import pandas as pd
import sys
import unittest
import phenorank



class test_phenorank(unittest.TestCase):
  """
  Tests:
    1: outputs pandas dataframe
    2: pandas dataframe contains the correct number of dimensions
    3: output contains correct genes
    5: correct human conditions identified as associated with each gene, missing conditions are NaN
    6: correct mouse conditions identified as associated with each gene, missing conditions are NaN
    7: correctly computes simGIC
    8: correctly propagates simGIC (just for none)
    9: correctly computes p-value (just for none)
    10: masking works
  """

  def test_phenorank_no_masking(self):

    logger = logging.getLogger(__name__)
    format = '%(message)s'
    logging.basicConfig(stream=sys.stdout, format=format, level=logging.DEBUG)

    # run function
    output = phenorank.phenorank.run_phenorank("DOID:1", nperm=10, r=0.5, ni=1, dir_data="test/data_phenorank")

    # run tests
    self.assertTrue(isinstance(output, pd.DataFrame)) # 1
    self.assertItemsEqual(output.shape, (5, 6)) # 2

    # 3
    self.assertEquals(output["GENE"]["ENSG1"], "ENSG1")
    self.assertEquals(output["GENE"]["ENSG2"], "ENSG2")
    self.assertEquals(output["GENE"]["ENSG3"], "ENSG3")
    self.assertEquals(output["GENE"]["ENSG4"], "ENSG4")
    self.assertEquals(output["GENE"]["ENSG7"], "ENSG7")

    # 5
    self.assertEquals(output["CONDITIONS_ALL"]["ENSG1"], "DOID:1")
    self.assertEquals(output["CONDITIONS_ALL"]["ENSG2"], "DOID:2")
    self.assertEquals(output["CONDITIONS_ALL"]["ENSG3"], "")
    self.assertEquals(output["CONDITIONS_ALL"]["ENSG4"], "")
    self.assertEquals(output["CONDITIONS_ALL"]["ENSG7"], "")

    #self.assertEquals(output["MOUSE_ALL"]["ENSG1"], "MGI:1")
    #self.assertEquals(output["MOUSE_ALL"]["ENSG2"], "")
    #self.assertEquals(output["MOUSE_ALL"]["ENSG3"], "MGI:4")
    #self.assertEquals(output["MOUSE_ALL"]["ENSG4"], "MGI:5")
    #self.assertEquals(output["MOUSE_ALL"]["ENSG7"], "MGI:2")

    # 6
    self.assertAlmostEqual(output["SCORE_RANKED_PROP"]["ENSG1"], 5.0, places=3)
    self.assertAlmostEqual(output["SCORE_RANKED_PROP"]["ENSG2"], 4.0, places=3)
    self.assertAlmostEqual(output["SCORE_RANKED_PROP"]["ENSG3"], 2.5, places=3)
    self.assertAlmostEqual(output["SCORE_RANKED_PROP"]["ENSG4"], 2.5, places=3)
    self.assertAlmostEqual(output["SCORE_RANKED_PROP"]["ENSG7"], 1.0, places=3)

  def test_phenorank_masking(self):

    # run function
    output = phenorank.phenorank.run_phenorank("DOID:1", nperm=10, r=0.9, ni=1, gene_mask="ENSG1", dir_data="test/data_phenorank")

    # run tests
    self.assertTrue(isinstance(output, pd.DataFrame)) # 1
    self.assertItemsEqual(output.shape, (5, 6)) # 2

    # 5, 6
    self.assertEquals(output["CONDITIONS_ALL"]["ENSG1"], "")
    self.assertEquals(output["CONDITIONS_ALL"]["ENSG2"], "DOID:2")
    self.assertEquals(output["CONDITIONS_ALL"]["ENSG3"], "")
    self.assertEquals(output["CONDITIONS_ALL"]["ENSG4"], "")
    self.assertEquals(output["CONDITIONS_ALL"]["ENSG7"], "")

    #self.assertEquals(output["MOUSE_ALL"]["ENSG1"], "MGI:1")
    #self.assertEquals(output["MOUSE_ALL"]["ENSG2"], "")
    #self.assertEquals(output["MOUSE_ALL"]["ENSG3"], "MGI:4")
    #self.assertEquals(output["MOUSE_ALL"]["ENSG4"], "MGI:5")
    #self.assertEquals(output["MOUSE_ALL"]["ENSG7"], "MGI:2")

    # 8, 11
    # these exact results have been calculated by hand
    self.assertAlmostEqual(output["SCORE_RANKED_PROP"]["ENSG1"], 4.0, places=3)
    self.assertAlmostEqual(output["SCORE_RANKED_PROP"]["ENSG2"], 5.0, places=3)
    self.assertAlmostEqual(output["SCORE_RANKED_PROP"]["ENSG3"], 2.0, places=3)
    self.assertAlmostEqual(output["SCORE_RANKED_PROP"]["ENSG4"], 2.0, places=3)
    self.assertAlmostEqual(output["SCORE_RANKED_PROP"]["ENSG7"], 2.0, places=3)



if __name__ == "__main__":
  unittest.main()
