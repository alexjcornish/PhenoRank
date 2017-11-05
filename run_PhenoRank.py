#!/usr/bin/env python -B

#==============================================================================
#title          :run_PhenoRank.py
#description    :run PhenoRank from the command line
#author         :Alex J Cornish
#date_created   :12 May 2015
#version        :0.1.2
#usage          :
#python_version :2.7.9
#==============================================================================

import datetime
import getopt
import logging
import numpy as np
import pandas as pd
import sys
import phenorank

def main(argv):
  """
  run PhenoRank from the command line

  Returns:
    nothing
  """

  # default arguments
  omim_obs = None # OMIM ID
  phenotypes_obs = None
  nperm = 1000
  r = 0.1
  ni = 20
  gene_mask = None
  human_only = False
  mouse_only = False

  # assign arguments to variables
  opts, args = getopt.getopt(argv, "o:d:p:n:r:i:g:hm", ["file_output=", "omim_id=", "phenotype_ids=", "nperm=", "r=", "ni=", "gene_mask=", "human_only", "mouse_only"])
  for opt, arg in opts:
    if opt in ("-o", "--file_output"): filename_output = arg
    elif opt in ("-d", "--omim_id"): omim_obs = arg
    elif opt in ("-p", "--phenotype_ids"): phenotypes_obs = arg.split(";")
    elif opt in ("-n", "--nperm"): nperm = int(arg)
    elif opt in ("-r", "--r"): r = float(arg)
    elif opt in ("-i", "--ni"): ni = int(arg)
    elif opt in ("-g", "--gene_mask"): gene_mask = arg
    elif opt in ("-h", "--human_only"): human_only = True
    elif opt in ("-m", "--mouse_only"): mouse_only = True

  # start log
  name_ascii = "  ___ _                 ___           _   \n | _ \ |_  ___ _ _  ___| _ \__ _ _ _ | |__\n |  _/ ' \/ -_) ' \/ _ \   / _` | ' \| / /\n |_| |_||_\___|_||_\___/_|_\__,_|_||_|_\_\ \n                                          "
  version = "0.1.1"
  logger.info(name_ascii)
  logger.info("Reducing study bias in gene prioritisation through simulation")
  logger.info("Version: {}".format(version))
  logger.info("Start time: {}".format(datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S')))

  # log input parameters
  logger.info("")
  logger.info("# Input parameters")
  logger.info("Input query disease OMIM term: {}".format(omim_obs))
  logger.info("Input query disease phenotype terms: {}".format(phenotypes_obs))
  logger.info("Number of permuations to complete: {}".format(nperm))
  logger.info("RWR restart probability (r): {}".format(r))
  logger.info("RWR number of iterations (ni): {}".format(ni))
  logger.info("Masked gene: {}".format(gene_mask))
  logger.info("Output filename: {}".format(filename_output))

  # check input
  if omim_obs is None and phenotypes_obs is None:
    raise Exception("Either an OMIM ID or a set of phenotype terms should be input. Neither have been input.")
  if omim_obs is not None and phenotypes_obs is not None:
    raise Exception("Either an OMIM ID or a set of phenotype terms should be input. Both have been input.")
  if mouse_only and human_only:
    raise Exception("Can't specify both only human and only mouse data")

  # run PhenoRank
  gene_scores = phenorank.phenorank.run_phenorank(omim_obs, phenotypes_obs, nperm, r, ni, gene_mask, include_h=not mouse_only, include_m=not human_only, filename_output=filename_output)
  gene_scores.to_csv(filename_output, header=True, sep="\t", index=False, columns=["GENE", "SCORE_UNRANKED_UNPROP", "SCORE_UNRANKED_PROP", "SCORE_RANKED_PROP", "PVALUE"])

  # complete log
  logger.info("")
  logger.info("All done, thanks for using PhenoRank!")



if __name__ == "__main__":
  logger = logging.getLogger(__name__)
  format = '%(message)s'
  logging.basicConfig(stream=sys.stdout, format=format, level=logging.DEBUG)
  main(sys.argv[1:])
