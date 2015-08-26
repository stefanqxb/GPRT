#!/usr/bin/python

import numpy as np
from matplotlib import pyplot as pp
pp.ion()

import data_tools
import feature_extraction
import ml_tools

if __name__=="__main__":
	# Loading the peptides from the file. We should give the filename to this function later.
	peptides = data_tools.read_data()

	benchmark = ml_tools.rt_benchmark( peptides, 'bow' )
	benchmark.cross_validation()
