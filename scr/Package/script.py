#!/usr/bin/python

import numpy as np
from matplotlib import pyplot as pp

pp.ion()

import data_tools
import ml_tools


if __name__ == "__main__":
    peptides = data_tools.read_data()
    # duplicated_message = data_tools.checked_duplicated(peptides)
    # print duplicated_message

    for n in [100, 200, 300, 500, 1000, 1500 ]:
        benchmark = ml_tools.rt_benchmark(peptides, 'elude', 'gp', n, 10)
        res = ml_tools.parallel_cross_validataion(benchmark)

        means = np.mean( res, axis=0 ).tolist()[0]
        stds = np.std( res, axis=0 ).tolist()[0]

        print n,means,stds


