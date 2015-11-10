#!/usr/bin/python

import platform
if platform.system() != 'Darwin':
    import matplotlib as mpl
    mpl.use('Agg')
import numpy as np
import data_tools
import pickle as pk
import ml_tools
import gp_tools
from matplotlib import pyplot as pp
pp.ion()
from sklearn.linear_model import LinearRegression
from common import parameters

class data_plotter:
    def __init__( self,n=100 ):
        self.params = parameters()
        self.n = n
    def load_data( self ):
        path = self.params.save_tmp % ( self.params.data_root, self.params.models_tag, self.n )
        self.peptides = data_tools.read_data( self.params.data_path )
        self.benchmark = ml_tools.rt_benchmark(self.peptides, 'elude', 'gp', self.n, self.params.nparts, self.params.train_ratio )
        self.models, self.kernels = ml_tools.load_rt_models( path )
    def rt_range( self ):
        rt = [];
        for p in self.peptides : 
            rt.append( p.rt )
        return np.min( rt ), np.max( rt )
    def benchmark_gp( self ):
        drt_values = []
        rmse_values = []

        for pind in range( self.benchmark.parts.nfolds ):
            [ drt, rmse ] = self.benchmark.eval_model(pind,self.models[pind])
            drt_values.append( drt )
            rmse_values.append( rmse )
        
        #drt_mean = np.mean( drt_values )
        #drt_std = np.std( drt_values )
        #rmse_mean = np.mean( rmse_values )
        #rmse_std = np.std( rmse_values )

        return [ self.n, drt_values, rmse_values ];

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

    def actual_vs_predicitive_variance( self ):
        rt_min,rt_max = self.rt_range()

        base = [];

        actual = [];
        predicitive = [];
    
        for pind in range( self.benchmark.parts.nfolds ):
            print pind
            a,p,s = self.benchmark.predict( pind, self.models[pind] )
            v = np.sqrt(s)
            b, av, pv  = ml_tools.actual_vs_predictive_variance( a,p,v,rt_min,rt_max,50 )

            base = b;
            actual.append( av )
            predicitive.append( pv )
            

        actual = np.matrix( actual )
        predicitive = np.matrix( predicitive )

        av_m = np.mean( actual, axis=0 )
        av_s = np.std( actual, axis=0 )
        pv_m = np.mean( predicitive, axis=0 )
        pv_s = np.std( predicitive, axis=0 )

        av_m = np.array( np.squeeze( np.asarray( av_m ) ) )
        av_s = np.array( np.squeeze( np.asarray( av_s ) ) )
        pv_m = np.array( np.squeeze( np.asarray( pv_m ) ) )
        pv_s = np.array( np.squeeze( np.asarray( pv_s ) ) )

        pp.figure()

        pp.plot( base, av_m, 'r' )
        pp.fill_between( base, av_m-av_s , av_m+av_s, alpha=0.8, facecolor='0.75' )

        pp.plot( base, pv_m, 'b' )
        pp.fill_between( base, pv_m-pv_s , pv_m+pv_s, alpha=0.8, facecolor='0.75' )

        [ x1,x2,y1,y2 ] = pp.axis();
        pp.axis( [ np.min(base), np.max(base), y1,y2 ])
        pp.grid()

        pp.xlabel('Actual RT')
        pp.ylabel('Variance')

        pp.legend(['Actual Variance','Predicitive Variance'],loc=4)

        pp.savefig('actual_vs_predicted_variance.pdf')



    def dummy():
        a,p,s = self.benchmark.predict( pind, self.models[pind] ) 
        var = np.sqrt(s)

        min_a = np.min(a)
        max_a = np.max(a)

        n = 50

        print min_a, max_a

        s = ( max_a - min_a )/(n-1)

        base = [];
        means = [];
        vars_gp = [];
        vars_dist = []
        low = [];
        high = [];

        
        base_pol = [];

        for i in range(n-2):
            
            ii = i+1

            b0 = (ii-1)*s+min_a
            b1 = (ii+1)*s+min_a
            c = (ii)*s+min_a 
            inds = np.where( (a >= b0) & (a<=b1) )[0]

            p_sub = p[inds]
            v_sub = var[inds] 

            m = np.mean( p_sub )
            v = np.std( p_sub )

            pol = [];

            for j in range(4):
                pol.append( c ** j );

            base.append( c )
            base_pol.append( pol )
            means.append(m)
            vars_gp.append( np.mean( v_sub ))
            vars_dist.append( v ) 
            low.append( m-2*v )
            high.append( m+2*v )

        base_pol = np.matrix( base_pol )

        base_mat = np.matrix(base).T
        LR_mean = LinearRegression()
        LR_mean.fit( base_pol,means )
        means_p = LR_mean.predict( base_pol )

        LR_low = LinearRegression()
        LR_low.fit( base_pol, low )
        low_p = LR_low.predict( base_pol )

        LR_high = LinearRegression()
        LR_high.fit( base_pol, high )
        high_p = LR_high.predict( base_pol ) 

        pp.figure()

        pp.plot( a,p,'r.')

        plt1 = pp.plot( base, means_p,linewidth=2,label='Local Mean' )
        plt2 = pp.plot( base, low_p,'k',linewidth=2,label='Local Mean + 2*variance' )
        plt3 = pp.plot( base, high_p,'g',linewidth=2, label='Local Up - 2*variance' )

        pp.grid();
        pp.xlabel('Actual RT')
        pp.ylabel('Predicted RT')
        pp.legend(loc=4)

        pp.savefig('actual_vs_predicted.pdf')

        pp.figure()

        pp.plot( base, (high_p - low_p)/4 )
        pp.plot( base, vars_gp )
        
        pp.grid()
        pp.xlabel('Actual RT')
        pp.ylabel('Variance')

        pp.legend(['Overtime Variance','GP Variance'],loc=4)

        pp.savefig('actual_vs_predicted_variance.pdf')


