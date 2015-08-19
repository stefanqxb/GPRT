#!/usr/bin/python

import numpy as np
from matplotlib import pyplot as pp
pp.ion()

import data_tools
import feature_extraction

if __name__=="__main__":
	peptides = data_tools.read_data('elude') # elude or bow
	#fe = feature_extraction.feature_extractor(peptides)
	#voc = fe.bow_voc(20)

temp = peptides[0]
temp1 = peptides[1]
temp2 = peptides[2]
#jully = feature_extraction.peptide.bow_descriptor(peptides[0], voc)




	#fe.gen_features_parallel( voc )

