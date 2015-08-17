#!/usr/bin/python

import numpy as np
from matplotlib import pyplot as pp
pp.ion()

import data_tools
import feature_extraction

if __name__=="__main__":
	peptides = data_tools.read_data()
	fe = feature_extraction.feature_extractor(peptides)
	voc = fe.bow_voc(15)

	fe.gen_features_parallel( voc )

