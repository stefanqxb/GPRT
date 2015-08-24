#!/usr/bin/python

import numpy as np
from matplotlib import pyplot as pp
pp.ion()

import data_tools
import feature_extraction

if __name__=="__main__":
	peptides = data_tools.read_data() 


	# The peptides that are supposed to be used for building feature extraction model are given to
	# feature extractor. Later these peptides will be only the training set.
	fe = feature_extraction.feature_extractor(peptides)

	# Elude model is built only based on the peptides given to feature_extractor
	em = fe.get_elude_model()
	# Bow Vocabulary is built only based on the peptides that are given to feature_extractor
	voc = fe.get_bow_voc( 25 );

	# These descritors are not normalized. A seperate module will be responsible for the normalization.
	elude_desc = peptides[0].elude_descriptor( em )
	bow_desc = peptides[0].bow_descriptor( voc )
