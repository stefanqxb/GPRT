#!/usr/bin/python

import numpy as np
import data_tools
import ml_tools


if __name__ == "__main__":
    peptides = data_tools.read_data()
    # duplicated_message = data_tools.checked_duplicated(peptides)
    # print duplicated_message

    bench = ml_tools.rt_benchmark(peptides, 'elude', 'gp', 100, 5)

    fmat = [];
    mmat = [];
    dmat = [];

    for i in range( bench.parts.nfolds ):
        print i
        model = bench.train_model(i)
        f,m,d = bench.test_sorted(i,model)

        fmat.append(f)
        mmat.append(m)
        dmat.append(d)

    fmat = np.matrix(fmat)
    mmat = np.matrix(mmat)
    dmat = np.matrix(dmat)
