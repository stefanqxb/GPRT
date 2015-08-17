#!/usr/bin/python

import numpy as np
from matplotlib import pyplot as pp
import feature_extraction as extraction
pp.ion()

def read_data():
	data_path = "/Users/heydar/Stuff/tmp/rt/20120126_EXQ5_KiSh_SA_LabelFree_HeLa_Phospho_Control_rep1_Fr3.rtimes.tsv"

	# Loading Data
	ff = open( data_path, 'r' )
	lines = ff.readlines()
	ff.close()

	peptides = [];

	for l in lines :
		values = l[:-2].split('\t')
		peptides.append( extraction.peptide( values[0], float(values[1]) ) );

	return peptides
