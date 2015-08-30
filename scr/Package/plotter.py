#!/usr/bin/python

from matplotlib import pyplot as pp
pp.ion()

def gp_vs_svm():
	n = [100,200,300,500,1000,2000,3000]
	gp_m = [0.224127,0.202901,0.16615,0.163617,0.155707,0.158368,0.158403]
	gp_s = [0.0514926,0.0564519,0.0190633,0.0158663,0.00379348,0.00552266,0.00474206]

	svr_m = [0.306034,0.226197,0.173972,0.173279,0.153411,0.155354,0.157249]
	svr_s = [0.0202813,0.0635182,0.0384707,0.0397791,0.0062629,0.00728697,0.00528372]

	pp.figure()

	pp.plot(n,gp_m,'r')
	pp.plot(n,svr_m,'b')
	pp.xlabel('Size of training data')
	pp.ylabel('Delta t')
	pp.legend(['GP','SVR'])
	pp.title('GP vs SVR')
	pp.grid()

	pp.savefig('gp_vs_svr.pdf')
