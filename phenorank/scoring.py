#!/usr/bin/env python -B

#==============================================================================
#title          :scoring.py
#description    :functions related to scoring
#author         :Alex J Cornish
#date_created   :06 November 2015
#version        :0.1.0
#usage          :
#python_version :2.7.9
#==============================================================================

import copy
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.stats as ss


def score_genes(disease_pheno, pheno_ancestors, pheno_ic, omim_ic, pheno_omim_ic_matrix, gc_h, gc_m, W, r, ni, include_h=True, include_m=True):
  """
  score genes in W_genes by the phenotypic relevance to the phenotype terms in pheno, and then propagate these across the network W
  this is the core part of the PhenoRank algorithm

  Args:
    disease_pheno: indices of phenotypes associated with disease (without ancestors)
    pheno_ancestors: a dictionary of term indices (keys) with their ancestor terms indices (values)
    pheno_ic: an numpy array of ICs for each indexed phenotype
    omim_ic:  a numpy array of ICs for each indexed OMIM ID
    pheno_omim_ic_matrix: a csr sparse matrix with rows representing OMIM IDs and columns represent phenotypes and the phenotype IC if the phenotype is associated with the OMIM ID and nothing otherwise
    gc_h, gc_m: sparse matrices containing 1 if the gene is associated with the OMIM ID and 0 otherwise for humans and mice
    W: a sparse matrix column-normalised (all columns sum to 1) adjacency matrix
    r: the restart probability for the score propagation
    ni: the number of interations to complete for the score propagation
    include_h, include_m: if True, include human data and mouse data respectively


  Returns:
    numpy array of score rankings for the genes in W_genes
  """

  if len(disease_pheno) == 0: raise ValueError("1 or more phenotype terms must be specified")

  disease_pheno = add_ancestors(disease_pheno, pheno_ancestors)

  disease_ic = sum(pheno_ic[disease_pheno])

  omim_scores = compute_simgic(disease_pheno, disease_ic, omim_ic, pheno_omim_ic_matrix)

  gene_scores_h = gc_h * omim_scores
  gene_scores_m = gc_m * omim_scores

  n_omims_h = np.array(gc_h.sum(1))[:,0]
  n_mutants_m = np.array(gc_m.sum(1))[:,0]
  n_omims_h[n_omims_h == 0] = 1
  n_mutants_m[n_mutants_m == 0] = 1

  if include_h and include_m:
    gene_scores = gene_scores_h / n_omims_h + gene_scores_m / n_mutants_m
  if include_h and not include_m:
    gene_scores = gene_scores_h / n_omims_h
  if not include_h and include_m:
    gene_scores = gene_scores_m / n_mutants_m

  gene_scores_prop = propagate_scores(gene_scores, W, r, ni)

  return gene_scores, gene_scores_prop, ss.rankdata(gene_scores_prop)



def add_ancestors(terms, ancestors):
  """
  expand a list of terms to include their ancestors

  Args:
    terms: list of term indices
    ancestors: a dictionary of term indices (keys) with their ancestor terms indices (values)

  Returns:
    list of term indices
  """

  terms_expanded = copy.copy(terms)
  for term in terms:
    try:
      terms_expanded += ancestors[term]
    except KeyError:
      pass
  return list(set(terms_expanded))



def compute_simgic(disease_pheno, disease_ic, omim_ic, pheno_omim_ic_matrix):
  """
  compute the simGIC scores for a set of phenotypes associated with the disease and the set of phenotypes associated with each OMIM ID

  Args:
    disease_pheno: indices of phenotypes associated with disease
    disease_ic: the IC of the disease
    omim_ic: an numpy array of ICs for each indexed OMIM ID
    pheno_omim_ic_matrix: a csr sparse matrix with rows representing OMIM IDs and columns represent phenotypes and the phenotype IC if the phenotype is associated with the OMIM ID and nothing otherwise

  Returns:
    a numpy array of simGIC scores
  """

  # compute the IC of the intersection
  intersection_ic = np.array(pheno_omim_ic_matrix[disease_pheno].sum(0))[0]

  # ensure that we don't divide by 0
  np.seterr(divide='raise')
  try:
    simgic = intersection_ic / (disease_ic + omim_ic - intersection_ic)
  except FloatingPointError:
    union_ic = disease_ic + omim_ic - intersection_ic
    union_ic[union_ic == 0] = 0.1
    simgic = intersection_ic / union_ic

  # return simGIC
  return simgic



def simulate_disease(n, pheno_cooccur):
  """
  simulate disease of n phenotype terms using sets of co-occuring phenotypes

  Args:
    n: the number of phenotype terms associated with the original disease (without expansion with ancestor terms)
    pheno_cooccur: a dictionary of dictionaries, containing the probabilities of selecting each phenotype when selecting a term co-occuring with the phenotype
  Returns:
    a list of phenotype term indices for the simulated disease (without expansion with ancestor terms)
  """

  # setup
  if n < 1: raise Exception("n should be greater than or equal to 1")

  # loop until a phenotype is chosen with enough co-occuring phenotypes (or we just give up and settle with what we have)
  i = 0
  loop_limit = 999
  while i < loop_limit:

    pheno_set = list()
    probs_pheno = pheno_cooccur[np.random.choice(pheno_cooccur.keys(), 1)[0]]
    probs = probs_pheno.keys()
    probs.sort(reverse=True)

    # loop until we have added enough phenotypes
    j = 0
    while j < len(probs_pheno):

      if len(probs_pheno[probs[j]]) <= (n - len(pheno_set)):
        pheno_set += probs_pheno[probs[j]] # add the complete set
      else:
        pheno_set += list(np.random.choice(probs_pheno[probs[j]], n - len(pheno_set))) # sample from the set

      # if the pheno_set is of the correct size, continue
      if len(pheno_set) == n:
        break

      j += 1

    if len(pheno_set) == n:
      break

    i += 1
    if i == loop_limit: raise Error("failed to find phenotype that co-occurs with enough phenotypes")

  return pheno_set



def propagate_scores(p, W, r=0.5, ni=10):
  """
  propagate score using the random walk with restart (rwr) method

  Args:
    p: a numpy array of gene scores, matching the order of genes in W
    W: a sparse matrix column-normalised (all columns sum to 1) adjacency matrix (the probability of moving from vertex i to j should be W[j][i])
    r: the restart probability
    ni: the number of interations to complete

  Returns:
    a numpy array of scores
  """

  if not isinstance(p, np.ndarray): raise TypeError("p should be a numpy array")
  if not sp.issparse(W): raise TypeError("W should be a sparse matrix")
  if len(p) != W.shape[0]: raise ValueError("number of scores not equal to the number of rows in W")
  if len(p) != W.shape[1]: raise ValueError("number of scores not equal to the columns of rows in W")

  # propagate scores using the rwr method
  p0 = np.copy(p)
  for _ in range(int(ni)):
    p = (1 - r) * (W * p) + r * p0

  return p
