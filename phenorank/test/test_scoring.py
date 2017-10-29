#!/usr/bin/env python -B

#==============================================================================
#title          :test_scoring.py
#description    :unit test functions in scoring.py
#author         :Alex J Cornish
#date_created   :06 November 2015
#version        :0.1.0
#usage          :
#python_version :2.7.9
#==============================================================================

import numpy as np
import pandas as pd
import scipy.sparse as sp
import unittest
import phenorank.scoring



class test_score_genes(unittest.TestCase):
  """
  Tests:
    1: returns list
    2: produces correct ranks given input
    3: handles genes that are unconnected in the network correctly
    4: error if pheno is of length 0
  """

  def test_score_genes_basic(self):
    # setup
    disease_pheno = [0, 1]
    pheno_ancestors = {0: [1]}
    pheno_ic = np.array([0.5, 0.5, 0.1, 0.2, 0.3])
    omim_ic = np.array([1.0, 0.1, 0.2, 0.5, 0.3])
    pheno_omim_ic_matrix = sp.csc_matrix(np.array([[0.5, 0.0, 0.0, 0.0, 0.0], [0.5, 0.0, 0.0, 0.5, 0.0], [0.0, 0.1, 0.0, 0.0, 0.0], [0.0, 0.0, 0.2, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.3]]))
    gc_h = sp.csr_matrix(np.array([[1, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]]))
    gc_m = sp.csr_matrix(np.array([[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 1, 1]]))
    W = sp.csc_matrix(np.array([[0.0, 0.5, 0.0, 0.0], [1.0, 0.0, 1.0, 0.0], [0.0, 0.5, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0]]))

    correct = {
      "r05_ni0": [4.0, 1.5, 1.5, 3.0],
      "r05_ni2": [4.0, 2.0, 1.0, 3.0],
      "r00_ni1": [1.5, 4.0, 1.5, 3.0]
    }

    # run function
    output = {}
    foo, bar, output["r05_ni0"] = phenorank.scoring.score_genes(disease_pheno, pheno_ancestors, pheno_ic, omim_ic, pheno_omim_ic_matrix, gc_h, gc_m, W, r=0.5, ni=0)
    foo, bar, output["r05_ni2"] = phenorank.scoring.score_genes(disease_pheno, pheno_ancestors, pheno_ic, omim_ic, pheno_omim_ic_matrix, gc_h, gc_m, W, r=0.5, ni=2)
    foo, bar, output["r00_ni1"] = phenorank.scoring.score_genes(disease_pheno, pheno_ancestors, pheno_ic, omim_ic, pheno_omim_ic_matrix, gc_h, gc_m, W, r=0.0, ni=1)

    # complete tests
    for key in output:
      self.assertTrue(isinstance(output[key], np.ndarray)) # 1
      self.assertEquals(list(output[key]), list(correct[key])) # 2, 3

    with self.assertRaises(ValueError):
      phenorank.scoring.score_genes([], pheno_ancestors, pheno_ic, omim_ic, pheno_omim_ic_matrix, gc_h, gc_m, W, r=0.0, ni=1) # 4

  #def test_score_genes_normalised(self):
    # setup
    disease_pheno = [0, 1]
    pheno_ancestors = {0: [1]}
    pheno_ic = np.array([0.1, 0.4, 1.1, 1.1])
    omim_ic = np.array([0.5, 1.5, 1.5, 0.5, 1.5, 1.5])
    pheno_omim_ic_matrix = sp.csr_matrix(np.array([[0.1, 0.4, 0.0, 0.0], [0.0, 0.4, 1.1, 0.0], [0.0, 0.4, 0.0, 1.1], [0.1, 0.4, 0.0, 0.0], [0.0, 0.4, 1.1, 0.0], [0.0, 0.4, 0.0, 1.1]]).T)
    gc_h = sp.csr_matrix(np.array([[1, 0, 0, 0, 0, 0], [1, 1, 0, 0, 0, 0], [1, 1, 1, 0, 0, 0], [1, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]))
    gc_m = sp.csr_matrix(np.array([[0, 0, 0, 1, 1, 1], [0, 0, 0, 1, 1, 0], [0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 1, 0]]))
    W = sp.csc_matrix(np.array([[0.0, 0.5, 0.0, 0.0, 0.0], [1.0, 0.0, 1.0, 0.0, 0.0], [0.0, 0.5, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 0.0, 1.0]]))

    correct = {"r05_ni0": [4.0, 1.5, 1.5, 3.0], "r05_ni2": [4.0, 2.0, 1.0, 3.0], "r00_ni1": [1.5, 4.0, 1.5, 3.0]}

    # run function
    output = {}
    foo, bar, output["r05_ni0"] = phenorank.scoring.score_genes(disease_pheno, pheno_ancestors, pheno_ic, omim_ic, pheno_omim_ic_matrix, gc_h, gc_m, W, r=0.5, ni=0)
    foo, bar, output["r05_ni2"] = phenorank.scoring.score_genes(disease_pheno, pheno_ancestors, pheno_ic, omim_ic, pheno_omim_ic_matrix, gc_h, gc_m, W, r=0.5, ni=2)
    foo, bar, output["r00_ni1"] = phenorank.scoring.score_genes(disease_pheno, pheno_ancestors, pheno_ic, omim_ic, pheno_omim_ic_matrix, gc_h, gc_m, W, r=0.0, ni=1)



class test_add_ancestors(unittest.TestCase):
  """
  Tests:
    1: returns list
    2: ancestors added if found
    3: ancestors not added if not found
    4: input phenotypes (leaves) included
    5: returned values unique
    6: if input is an empty list, output is an empty list
  """

  def test_add_ancestors_toy(self):
    # setup
    terms = [0, 1, 2]
    ancestors = {0: [1, 3], 1: [4, 5], 3: [6, 7]}

    # run function
    output = phenorank.scoring.add_ancestors(terms, ancestors)

    # run tests
    # 1, 2, 3, 4, 5
    self.assertEquals(type(output), type([]))
    self.assertEquals(len(output), 6)
    self.assertItemsEqual(output, [0, 1, 2, 3, 4, 5])

  def test_add_ancestors_empty(self):
    # setup
    terms = []
    ancestors = {0: [1, 2]}

    # run function
    output = phenorank.scoring.add_ancestors(terms, ancestors)

    # run tests
    # 6
    self.assertEquals(type(output), type([]))
    self.assertEquals(len(output), 0)



class test_compute_simgic(unittest.TestCase):
  """
  Tests:
    1: returns list
    2: computes correct score (1) when intersection complete
    3: computes correct score when intersection parial
    4: computes correct score (0) when intersection empty
  """

  def test_compute_simgic_toy(self):
    # setup
    disease_pheno = [0, 1]
    disease_ic = 3.0
    ic_omim = np.array([3.0, 5.0, 4.0])
    pheno_omim_ic_matrix = np.array([[1.0, 0, 0], [2.0, 2.0, 0], [0.0, 3.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 5.0]])
    pheno_omim_ic_matrix = sp.csr_matrix(pheno_omim_ic_matrix)

    # run function
    output = phenorank.scoring.compute_simgic(disease_pheno, disease_ic, ic_omim, pheno_omim_ic_matrix)

    # run tests
    self.assertTrue(isinstance(output, np.ndarray)) # 1
    self.assertEquals(output[0], 1.0) # 2
    self.assertEquals(output[1], 1.0 / 3) # 3
    self.assertEquals(output[2], 0.0) # 4



class test_simulate_disease(unittest.TestCase):
  """
  Tests:
    1: returns list
    2: returns correct number of terms
    3: produces correct results when size of the set of the most frequently co-occuring terms is greater than the number of terms to select from
    4: produces correct results when size of the set of the most frequently co-occuring terms is equal to the number of terms to select from
    5: produces correct results when size of the set of the most frequently co-occuring terms is smaller than the number of terms to select from
    6: if no phenotype term co-occurs with enough other phenotype terms, then error produced
  """

  def test_simulate_disease(self):
    # complete all of the tests with the uniform pheno_assoc

    # setup
    pheno_cooccur = {0: {0.4: [3, 2], 0.15: [1], 0.05: [0]}, 1: {}, 2:{}, 3:{}, 4:{}}

    # run function
    output_n1 = phenorank.scoring.simulate_disease(1, pheno_cooccur)
    output_n2 = phenorank.scoring.simulate_disease(2, pheno_cooccur)
    output_n3 = phenorank.scoring.simulate_disease(3, pheno_cooccur)
    output_n4 = phenorank.scoring.simulate_disease(4, pheno_cooccur)

    # run tests
    self.assertTrue(isinstance(output_n1, list)) # 1
    self.assertTrue(isinstance(output_n2, list)) # 1
    self.assertEquals(len(output_n1), 1) # 2
    self.assertEquals(len(output_n2), 2) # 2
    self.assertEquals(len(output_n3), 3) # 2
    self.assertEquals(len(output_n4), 4) # 2
    self.assertTrue(output_n1[0] in [2, 3])
    self.assertItemsEqual(output_n2, [2, 3]) # 3
    self.assertItemsEqual(output_n3, [1, 2, 3]) # 4
    self.assertItemsEqual(output_n4, [0, 1, 2, 3]) # 5

    with self.assertRaises(Exception):
      phenorank.scoring.simulate_disease(5, pheno_cooccur) # 6



class test_propagate_scores(unittest.TestCase):
  """
  Tests:
    1: returns numpy array of scores
    2: if restart probability is 0, same as rw
    3: correctly uses restart probability of 0.5
    4: if restart probability is 1, input scores same as output
    5: if number of iterations is 0, input scores same as output
    6: correctly uses number of iterations of 2
    7: lone vertices correctly scored
    8: if p is not a numpy array, then an error is produced
    9: if W is not a sprase matrix, then an error is produced
    10: if the number of scores is not equal to the number of rows and columns in W, then an error is produced
  """

  def test_propagate_scores_toy(self):
    """
    NETWORK
    0 - 1 - 2 - 3    5
    """

    # setup
    p = np.array([1.0, 0.0, 0.0, 0.0, 2.0])
    W_dense = np.transpose(np.array([[0.0, 1.0, 0.0, 0.0, 0.0], [0.5, 0.0, 0.5, 0.0, 0.0], [0.0, 0.5, 0.0, 0.5, 0.0], [0.0, 0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 1.0]]))
    W = sp.csc_matrix(W_dense)

    correct = {}
    correct["r0"] = np.array([0.0, 1.0, 0.0, 0.0, 2.0])
    correct["r05"] = np.array([0.5, 0.5, 0.0, 0.0, 2.0])
    correct["r1"] = p
    correct["ni0"] = p
    correct["ni2"] = np.array([0.625, 0.25, 0.125, 0.0, 2.0])

    # run functions
    outputs = {}
    outputs["r0"] = phenorank.scoring.propagate_scores(p, W, r=0, ni=1)
    outputs["r05"] = phenorank.scoring.propagate_scores(p, W, r=0.5, ni=1)
    outputs["r1"] = phenorank.scoring.propagate_scores(p, W, r=1, ni=1)
    outputs["ni0"] = phenorank.scoring.propagate_scores(p, W, r=0.5, ni=0)
    outputs["ni2"] = phenorank.scoring.propagate_scores(p, W, r=0.5, ni=2)

    # run tests
    for output_name in correct:
      self.assertTrue(isinstance(outputs[output_name], np.ndarray)) # 1
      self.assertEquals(len(outputs[output_name]), len(correct[output_name]))
      for i in range(len(correct[output_name])):
        self.assertAlmostEqual(outputs[output_name][i], correct[output_name][i], places=3) # 2, 3, 4, 5, 6, 7

    with self.assertRaises(TypeError):
      phenorank.scoring.propagate_scores(list(p), W, r=0.5, ni=1) # 8

    with self.assertRaises(TypeError):
      phenorank.scoring.propagate_scores(p, W_dense, r=0.5, ni=1) # 9

    with self.assertRaises(ValueError):
      phenorank.scoring.propagate_scores(p[:3], W, r=0.5, ni=1) # 10



if __name__ == "__main__":
  unittest.main()
