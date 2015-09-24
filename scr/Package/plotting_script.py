#!/usr/bin/python

import platform
import numpy as np
import data_tools
import pickle as pk
import ml_tools
import gp_tools
from matplotlib import pyplot as pp
pp.ion()

class data_plotter:
    def __init__( self, n=100 ):
        self.n = 100
    def load_data( self ):
        if platform.system() == 'Darwin' :
            data_root = "/Users/heydar/Stuff/tmp/gprt"
        else :
            data_root = "/media/hdd/heydar/data/gprt"
        model_tmp = "%s/models_ntrain_%d.pk"
        path = model_tmp % ( data_root, self.n )
        self.peptides = data_tools.read_data()
        self.benchmark = ml_tools.rt_benchmark(peptides, 'elude', 'gp', self.n, 5)
        self.models, self.kernels = ml_tools.load_rt_models( path )

        test_vectors = benchmark.get_test_vectors( 0, models[0] )

        return models[0], kernels[0]

    def dist_vs_var( pind=0 ):
        test_vectors = self.benchmark.get_test_vectors(pind,self.models[pind])
        dist = [];
        var = [];

        for t in test_vectors : 
            d = self.kernels[pind].mean_dist( t[0] )
            dist.append( d )
            var.append(t[3])

        dist = np.array( dist )
        var = np.array( var )

        pp.plot( dist, var, 'r.')
        pp.grid()
        pp.xlabel('Distance From Training Data')
        pp.ylabel('Predicted Variance')
        pp.title('Distance vs. Variance')
        pp.savefig('dist_vs_var.pdf')
