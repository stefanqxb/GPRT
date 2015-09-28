#!/usr/bin/python

import matplotlib as mpl
mpl.use('Agg')
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
        self.n = n
    def load_data( self ):
        if platform.system() == 'Darwin' :
            data_root = "/Users/heydar/Stuff/tmp/gprt"
        else :
            data_root = "/media/hdd/heydar/data/gprt"
        model_tmp = "%s/models_ntrain_%d.pk"
        path = model_tmp % ( data_root, self.n )
        self.peptides = data_tools.read_data()
        self.benchmark = ml_tools.rt_benchmark(self.peptides, 'elude', 'gp', self.n, 10)
        self.models, self.kernels = ml_tools.load_rt_models( path )

    def dist_vs_var( self, pind=0 ):
        test_vectors = self.benchmark.get_test_vectors(pind,self.models[pind])
        dist = [];
        var = [];

        for t in test_vectors : 
            d = self.kernels[pind].min_dist( t[0] )
            dist.append( d )
            var.append(t[3])

        dist = np.array( dist )
        var = np.array( var )

        pp.figure()
        pp.plot( dist, var, 'r.')
        pp.grid()
        pp.xlabel('Distance From Training Data')
        pp.ylabel('Predicted Variance')
        pp.title('Distance vs. Variance')
        pp.savefig('dist_vs_var.pdf')

    def section_error_independent( self ):
        f_m, m_m, e_m, e_v = ml_tools.section_error_independent( self.benchmark, self.models )

        pp.figure()
        pp.plot( m_m, e_m,'.-' )
        pp.fill_between( m_m, e_m-e_v , e_m+e_v, alpha=0.8, facecolor='0.75' )

        [ x1,x2,y1,y2 ] = pp.axis();
        pp.axis( [ np.min(m_m), np.max(m_m), y1,y2 ])

        pp.xlabel('Predicted Standard Deviation')
        pp.ylabel('Root Mean Square Error')
        pp.title('Sections Error')
        pp.grid()
        pp.savefig('ind_section.pdf')

    def section_error_overall( self ):
        f_m, m_m, e_m, e_v = ml_tools.section_error_overall( self.benchmark, self.models )

        pp.figure()
        pp.plot( f_m, e_m,'.-' )
        pp.fill_between( f_m, e_m-e_v , e_m+e_v, alpha=0.8, facecolor='0.75' )

        [ x1,x2,y1,y2 ] = pp.axis();
        pp.axis( [ np.min(f_m), np.max(f_m), y1,y2 ])

        pp.xlabel('Data Fraction')
        pp.ylabel('Root Mean Square Error')
        pp.title('Cumulative Overall Error')
        pp.grid()
        pp.savefig('cum_overall.pdf')
