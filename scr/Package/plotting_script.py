#!/usr/bin/python

import numpy as np
import data_tools
import pickle as pk
import ml_tools
import gp_tools

def load_test_vectors( n = 100 ):
    data_root = "/Users/heydar/Stuff/tmp/gprt"
    model_tmp = "%s/models_ntrain_%d.pk"
    path = model_tmp % ( data_root, n )

    peptides = data_tools.read_data()
    benchmark = ml_tools.rt_benchmark(peptides, 'elude', 'gp', n, 5)
    models, kernels = ml_tools.load_rt_models( path )

    test_vectors = benchmark.get_test_vectors( 0, models[0] )
    return models[0], kernels[0], test_vectors

def calc_values( kernel, test_vectors ):

    dist = [];
    var = [];

    for t in test_vectors : 
        d = kernel.mean_dist( t[0] )
        dist.append( d )
        var.append(t[3])

    return np.array( dist ), np.array( var )
