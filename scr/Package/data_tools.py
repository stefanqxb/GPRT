#train_inds train_inds !/usr/bin/python

import numpy as np
from matplotlib import pyplot as pp
import feature_extraction as extraction
import sys
import data_manager as dm
import random
pp.ion()

def read_data():
	data_path = "./data/20121212_S25_3ug_300min_LCMS_PM3_Tryp_GRADIENT_15sec_MZ.rtimes.tsv"
    #data_path = "./data/20120126_EXQ5_KiSh_SA_LabelFree_HeLa_Phospho_Control_rep1_Fr3.rtimes.tsv"
	ff = open(data_path, 'r')
	lines = ff.readlines()
	ff.close()
	peptides = []
	for l in lines :
		values = l[:-2].split('\t')
		peptides.append(extraction.peptide(values[0], float(values[1])))
	return np.array(peptides)

	#if str == 'bow':
	#
	#    for l in lines:
	#        values = l[:-2].split('\t')
	#        peptides.append(extraction.peptide(values[0], float(values[1]), " "))
	#elif str == 'elude':
	#    for l in lines:
	#        values = l[:-2].split('\t')
	#        peptides.append(extraction.peptide(values[0],values[1], data_path))
