#!/usr/bin/python

import numpy as np
from matplotlib import pyplot as pp
pp.ion()

import data_tools
import feature_extraction

if __name__=="__main__":
	# Loading the peptides from the file. We should give the filename to this function later.
	peptides = data_tools.read_data()
	mg = feature_extraction.model_generator( peptides )
	# The peptides that are supposed to be used for building feature extraction model are given to
	# feature extractor. Later these peptides will be only the training set.    
	# Elude model is built only based on the peptides given to feature_extractor
	elude = mg.get_elude_model()
	# Bow Vocabulary is built only based on the peptides that are given to feature_extractor
	voc = mg.get_bow_voc( 25 );
	# These descritors are not normalized. A seperate module will be responsible for the normalization.
	fe = feature_extraction.feature_extractor()
	fe.ext_elude( peptides[0:500], elude )
