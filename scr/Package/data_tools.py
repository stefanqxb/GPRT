#train_inds train_inds !/usr/bin/python

import numpy as np
from matplotlib import pyplot as pp
import feature_extraction as extraction
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

def checked_duplicated(peptides):
    candiates = []
    rt = []
    for pp in peptides:
        if pp.sequence not in candiates:
            candiates.append(pp.sequence)
            rt.append(pp.rt)

    if len(candiates) == len(peptides):
        return "No duplicated peptide in this data set"
    else:
        for i in range(len(candiates)):
            for pp in peptides:
                if candiates[i] == pp.sequence and rt[i] != pp.rt:
                    return "Error: Same peptide has different RT"

    return "All same peptide has same RT"