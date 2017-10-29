#!/usr/bin/env python -B

#==============================================================================
#title          :inout.py
#description    :functions for inputting and outputting data
#author         :Alex J Cornish
#date_created   :08 May 2015
#version        :0.1.0
#usage          :
#python_version :2.7.9
#==============================================================================

import cPickle



def import_dictionary(con, split_by=None, key_int=False, value_int=False, value_float=False):
  """
  import dictionary of values that are not concatenated
  file should contain 2 columns: 1 containing the keys and 2 values

  Args:
    con: the connection
    split_by: if not None, the keys are split by these values
    value_float: if True then the values are converted to floats

  Returns:
    dictionary of key-value pairs, value is a single value
  """

  # import data
  lines = con.read().splitlines()

  # create dictionary
  d = {}

  for line in lines[1:]:
    key, value = line.split("\t")
    if key_int: key = int(key)
    if split_by is not None:
      values = value.split(split_by)
      if values == [""]: values= [] # CONTINUE
      if value_int: values = [int(value) for value in values]
      if value_float: values = [float(value) for value in values]
      d[key] = values
    else:
      if value_int: value = int(value)
      if value_float: value = float(value)
      d[key] = value

  return d
