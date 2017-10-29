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



def main(argv):
  """
  run PRINCE from the command line

  Returns:
    nothing
  """

  # default arguments
  gene_mask = None

  # assign arguments to variables
  opts, args = getopt.getopt(argv, "o:i:g:a:n:c:", ["file_output=", "condition=", "gene_mask=", "a=", "n=", "c="])
  for opt, arg in opts:
    if opt in ("-o", "--file_output"): filename_output = arg
    elif opt in ("-i", "--condition"): condition = arg
    elif opt in ("-g", "--gene_mask"): gene_mask = arg
    elif opt in ("-a", "--a"): a = arg
    elif opt in ("-n", "--n"): n = arg
    elif opt in ("-c", "--c"): c = arg

  # run prince
  gene_scores = run_prince(condition, a=a, n=n, c=c, gene_mask=gene_mask)

  # save results
  gene_scores.to_csv(filename_output, header=True, sep="\t", index=False)



def run_prince(condition, a, n, c, gene_mask=None, dir_data="data_prince"):
  """
  run the main PRINCE function

  Args:
    condition: the ID for the condition of interest
    a: the value of alpha parameter used
    n: the number of iterations to completed
    c: the value of c parameter used
    gene_mask: a gene to mask
    dir_data: path to the data, used for testing

  Returns:
    a pandas DataFrame of scores for the input GWAS loci
  """

  # import data
  con = open(pkg_resources.resource_filename('phenorank', dir_data + '/gc_h_omim.tsv'), "r")
  gc = inout.import_dictionary(con, split_by="|")
  con.close()

  con = open(pkg_resources.resource_filename('phenorank', dir_data + '/pheno_sim.tsv'), "r")
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
  condition_score = compute_condition_scores(condition, pheno_sim)
  gene_pheno_scores = compute_gene_pheno_scores(gc, condition_score)

  # propagate scores
  scores = np.array([gene_pheno_scores[gene] if gene in gene_pheno_scores else 0.0 for gene in genes])
  score_propagated_obs = scoring.propagate_scores(scores, W, r=a, ni=n)

  # collate scores
  res = pd.DataFrame({"GENE": pd.Series(genes, index=genes), "SCORE": pd.Series(score_propagated_obs, index=genes)})

  return res



def compute_condition_scores(id, pheno_sim):
  """
  score diseases by their phenotypic similarity to the disease ID

  Args:
    id: the disease ID (the same ID format as pheno_sim)
    pheno_sim: a pandas DataFrame of phenotype similarity scores

  Returns:
    a dictionary of disease ID-phenotype similarity score key-value pairs
  """

  keys = pheno_sim[id].index.tolist()
  values = pheno_sim[id].tolist()
  output = {keys[i]: values[i] for i in range(len(pheno_sim[id]))}

  return output



def compute_gene_pheno_scores(gc, condition_score):
  """
  score genes in gc by the maximum phenotype score of phenotypes in pheno_sim

  Args:
    gc: a dictionary of gene-OMIM list key-value pairs
    condition_score: a dictionary of OMIM-phenotype similarity scores for the disease of interest

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



if __name__ == "__main__":
  main(sys.argv[1:])
