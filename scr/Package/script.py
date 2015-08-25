#!/usr/bin/python

import numpy as np
from matplotlib import pyplot as pp
pp.ion()

import data_tools
import feature_extraction


data_path = "./data/20121212_S25_3ug_300min_LCMS_PM3_Tryp_GRADIENT_15sec_MZ.rtimes.tsv"

peptides = data_tools.read_data(data_path)
fe = feature_extraction.model_generator(peptides)
voc = fe.bow_voc(25)
aaAlphabet = fe.amino_list

str = "elude"

# forming feature:
# if "bow" then forming bow feature for each peptide and store in bow_descriptor,
# if "elude": forming elude feature for each peptide and store in elude_descriptor
# then use for training and testing.

if str == "bow":
    for pp in peptides:
        feature_extraction.peptide.bow_descriptor(pp, voc)
elif str == "elude":
    for pp in peptides:
        feature_extraction.peptide.elude_descriptor(pp, aaAlphabet, peptides) # adding a new item elude_descriptor in class peptide
        temp = pp.elude_descriptor # test to see if the elude_descriptor is working
        print temp



	#fe.gen_features_parallel( voc )

