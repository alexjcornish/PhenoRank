#!/usr/bin/env python -B

#==============================================================================
#title          :test_inout.py
#description    :unit test functions in inout.py
#author         :Alex J. Cornish
#date_created   :02 February 2015
#version        :0.1.0
#usage          :
#python_version :2.7.9
#==============================================================================

import unittest
import phenorank.inout


      
class test_import_dictionary(unittest.TestCase):
  """
  Tests:
    1: returns dictionary
    2: values are split by the correct character
    3: if values are not split, they are returned as integers by default
    4: values are correctly converted to floats
  """

  def test_import_dictionary_split(self):
    # run function
    con = open("data_phenorank/disease_pheno_toy.tsv", "r")
    output = phenorank.inout.import_dictionary(con, split_by="|")
    con.close()

    # run tests
    self.assertEquals(type(output), type({})) # 1
    self.assertEquals(len(output), 3)
    self.assertItemsEqual(output["DISEASE_ID:1"], []) # 2
    self.assertItemsEqual(output["DISEASE_ID:2"], ["PHENO_ID:1"]) # 2
    self.assertItemsEqual(output["DISEASE_ID:3"], ["PHENO_ID:2", "PHENO_ID:3"]) # 2

  def test_import_dictionary_not_split_integer(self):
    # run function
    con = open("data_phenorank/ic_toy.tsv", "r")
    output = phenorank.inout.import_dictionary(con, value_int=True)
    con.close()

    # run tests
    self.assertEquals(type(output), type({})) # 1
    self.assertEquals(len(output), 3) # 2
    self.assertEquals(output["ID:1"], int(1.0)) # 3
    self.assertEquals(output["ID:2"], int(2.0)) # 3
    self.assertEquals(output["ID:3"], int(0.0)) # 3

  def test_import_dictionary_not_split_float(self):
    # run function
    con = open("data_phenorank/ic_toy.tsv", "r")
    output = phenorank.inout.import_dictionary(con, value_float=True)
    con.close()

    # run tests
    self.assertEquals(type(output), type({})) # 1
    self.assertEquals(len(output), 3) # 2
    self.assertEquals(output["ID:1"], 1.0) # 3
    self.assertEquals(output["ID:2"], 2.0) # 3
    self.assertEquals(output["ID:3"], 0.0) # 3



if __name__ == "__main__":
  unittest.main()
