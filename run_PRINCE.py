#!/usr/bin/env python -B

#==============================================================================
#title          :run_PRINCE.py
#description    :run PRINCE from the command line
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
import phenorank
import sys

def main(argv):
  """
  run PRINCE from the command line

  Returns:
    nothing
  """

  # default arguments
  filename_output = None
  omim_obs = None
  a = 0.5
  n = 20
  c = -15
  gene_mask = None

  # assign arguments to variables
  opts, args = getopt.getopt(argv, "o:d:a:n:c:g:", ["file_output=", "omim_id=", "a=", "n=", "c=", "gene_mask="])
  for opt, arg in opts:
    if opt in ("-o", "--file_output"): filename_output = arg
    elif opt in ("-d", "--omim_id"): omim_obs = arg
    elif opt in ("-a", "--a"): a = float(arg)
    elif opt in ("-n", "--n"): n = int(arg)
    elif opt in ("-c", "--c"): c = float(arg)
    elif opt in ("-g", "--gene_mask"): gene_mask = arg

  # check input
  if filename_output is None:
    raise Exception("An output filename must be specified")
  if omim_obs is None:
    raise Exception("An OMIM ID must be specified")

  # start log
  logger.info("PRINCE")
  logger.info("Start time: {}".format(datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S')))

  # log input parameters
  logger.info("")
  logger.info("# Input parameters")
  logger.info("Input query disease OMIM term: {}".format(omim_obs))
  logger.info("Alpha parameter (a): {}".format(a))
  logger.info("RWR number of iterations (n): {}".format(n))
  logger.info("C parameter (c): {}".format(c))
  logger.info("Masked gene: {}".format(gene_mask))
  logger.info("Output filename: {}".format(filename_output))

  # run prince
  gene_scores = phenorank.prince.run_prince(omim_obs, a=a, n=n, c=c, gene_mask=gene_mask)
  gene_scores.to_csv(filename_output, header=True, sep="\t", index=False)

  # complete log
  logger.info("")
  logger.info("All done")

if __name__ == "__main__":
  logger = logging.getLogger(__name__)
  format = '%(message)s'
  logging.basicConfig(stream=sys.stdout, format=format, level=logging.DEBUG)
  main(sys.argv[1:])
