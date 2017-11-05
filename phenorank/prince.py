#!/usr/bin/env python -B

#==============================================================================
#title          :prince.py
#description    :functions for running PRINCE
#author         :Alex J Cornish
#date_created   :20 November 2015
#version        :0.1.0
#usage          :
#python_version :2.7.9
#==============================================================================

import cPickle
import getopt
import numpy as np
import pandas as pd
import pkg_resources
import scipy.sparse as sp
import scipy.stats as ss
import sys
import inout
import scoring

def compute_condition_scores(id, pheno_sim, c, factor):
  """
  score diseases by their phenotypic similarity to the disease ID

  Args:
    id: the disease ID (the same ID format as pheno_sim)
    pheno_sim: a pandas DataFrame of phenotype similarity scores
    c: the value of parameter c
    factor: pheno_sim values have been scaled to values between 0 and this value (to save space)

  Returns:
    a dictionary of disease ID-phenotype similarity score key-value pairs
  """

  d = np.log(9999) # defined in Vanunu et al. 2010 (PLOS Computational Biology)
  keys = pheno_sim[id].index.tolist()
  values = pheno_sim[id].tolist()
  values = 1 / (1 + (np.exp([c * i / factor for i in values] + d)))
  return {keys[i]: values[i] for i in range(len(pheno_sim[id]))}

def compute_gene_pheno_scores(gc, condition_score):
  """
  score genes in gc by the maximum phenotype-similarity score of conditions in condition_score

  Args:
    gc: a dictionary of gene-OMIM list key-value pairs
    condition_score: a dictionary of OMIM-phenotype similarity scores for the disease of interest
``
  Returns:
    a pandas Series of scores for each gene in gc
  """

  genes = gc.keys()
  scores = {gene: 0 for gene in genes}
  for gene in genes:
    for condition in gc[gene]:
      try:
        if condition_score[condition] > scores[gene]:
          scores[gene] = condition_score[condition]
      except KeyError:
        pass
  return pd.Series(scores.values(), index=scores.keys())

def run_prince(condition, a, n, c, factor=100, gene_mask=None, dir_data="data_prince"):
  """
  run the main PRINCE function

  Args:
    condition: the ID for the condition of interest
    a: the value of parameter alpha used
    n: the number of iterations to completed
    c: the value of parameter c used
    factor: pheno_sim values have been scaled to values between 0 and this value (to save space)
    gene_mask: a gene to mask
    dir_data: path to the data, used for testing

  Returns:
    a pandas DataFrame of scores for the input GWAS loci
  """

  # import data
  con = open(pkg_resources.resource_filename('phenorank', dir_data + '/gc_h_omim.tsv'), "r")
  gc = inout.import_dictionary(con, split_by="|")
  con.close()

  con = open(pkg_resources.resource_filename('phenorank', dir_data + '/phenosim.tsv'), "r")
  pheno_sim = pd.read_csv(con, sep="\t", header=0, index_col=0)
  con.close()

  con = open(pkg_resources.resource_filename('phenorank', dir_data + '/prince_genes.tsv'), "r")
  genes = list(pd.read_csv(con, header=None)[0])
  con.close()

  con = open(pkg_resources.resource_filename('phenorank', dir_data + '/W_generic.pickle'), "r")
  W = cPickle.load(con)
  con.close()

  # mask association between gene and DO ID of interest, if there is an association
  if gene_mask:
     try:
       gc[gene_mask].remove(condition)
     except (KeyError, ValueError):
       pass

  # score genes
  condition_score = compute_condition_scores(condition, pheno_sim, c, factor)
  gene_pheno_scores = compute_gene_pheno_scores(gc, condition_score)

  # propagate scores
  scores = np.array([gene_pheno_scores[gene] if gene in gene_pheno_scores else 0.0 for gene in genes])
  score_propagated_obs = scoring.propagate_scores(scores, W, r=a, ni=n)

  # return collated scores
  return pd.DataFrame({"GENE": pd.Series(genes, index=genes), "SCORE": pd.Series(score_propagated_obs, index=genes), "Y": pd.Series(scores, index=genes)})
