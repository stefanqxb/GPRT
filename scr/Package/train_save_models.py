#!/usr/bin/python

import platform
import numpy as np
import data_tools
import ml_tools
import pickle as pk

if __name__ == "__main__":
    if platform.system() == 'Darwin' :
        data_root = "/Users/heydar/Stuff/tmp/gprt"
    else :
        data_root = "/media/hdd/heydar/data/gprt"
    peptides = data_tools.read_data()
    for n in [100, 200, 300, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 10000 ]:
        benchmark = ml_tools.rt_benchmark(peptides, 'elude', 'gp', n, 10)
        models = ml_tools.parallel_train( benchmark )
        save_path = "%s/models_ntrain_%d.pk" % ( data_root, n )
        print save_path
        with open( save_path, 'w' ) as ff :
            pk.dump( [ models ], ff )
            ff.close()
    #    res = ml_tools.parallel_cross_validataion(benchmark)

    #    means1 = np.mean( res, axis=0 ).tolist()[0]
    #    stds1 = np.std( res, axis=0 ).tolist()[0]

    #    res = ml_tools.parallel_cross_validataion_multi(benchmark)

    #    means2 = np.mean( res, axis=0 ).tolist()[0]
    #    stds2 = np.std( res, axis=0 ).tolist()[0]
    #    print n,means1[0],stds1[0],means2[0],stds2[0]


