#!/usr/bin/python

import numpy as np
from matplotlib import pyplot as pp
pp.ion()

import data_tools
import ml_tools


if __name__== "__main__":
    peptides = data_tools.read_data()
    duplicated_message = data_tools.checked_duplicated(peptides)
    print duplicated_message
    benchmark = ml_tools.rt_benchmark(peptides, 'elude', 200)
    benchmark.cross_validation()