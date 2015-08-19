#!/usr/bin/python
from scipy.stats.stats import pearsonr

import matplotlib.pyplot as plt
import pylab

def mini_time_win(hist,diff,step,max_total_t):
    max_t = max(diff)
    min_t = min(diff)
    total = sum(hist)
    threshold = round(0.95 * total)
    count = 0
    counter = 0
    for i in range(len(hist)):
        count = hist[i] + count
        if count >= threshold:
            counter = i
            break

    time_interval = 2*counter * (max_t - min_t)/(step*max_total_t)
    return  time_interval

class evaluation:
    def evaluation(pv_gp,pv_svr,test_tag,test_tag_svr,t_gp,t_svr):
        step = 10
        diff_svr = abs(pv_svr - test_tag_svr)
        diff_gp = abs(pv_gp -test_tag)
        histo_svr = plt.hist(diff_svr,step)
        histo_gp = plt.hist(diff_gp,step)

        max_total_gp = max(test_tag) - min(test_tag)
        max_total_svr = max(test_tag_svr) - min(test_tag_svr)

        mtw_gp =mini_time_win(histo_gp[0],diff_gp,step,max_total_gp)
        corrcoef_gp = pearsonr(pv_gp,test_tag)
        mtw_svr = mini_time_win(histo_svr[0],diff_svr,step,max_total_svr)
        corrcoef_svr = pearsonr(pv_svr,test_tag_svr)

        print 'corrcoef of GP = ',corrcoef_gp[0],'corrcoef of SVR = ',corrcoef_svr[0]
        print 'mtw of GP = ',mtw_gp, 'mtw of SVR  =', mtw_svr
        print 'running time of GP = ',t_gp , 'sec', 'running time of SVR = ', t_svr
        return step ,max_total_gp,max_total_svr

    def ploting(pv_gp,pv_svr, test_tag,test_tag_svr,max_total_gp,max_total_svr,step):
        diff_svr = pv_svr - test_tag_svr
        diff_gp = pv_gp - test_tag
        pylab.figure(1)
        plt.subplot(121)
        plt.plot(pv_svr,test_tag_svr, 'g*')
        plt.plot([1,max_total_svr],[1,max_total_svr],'y-',linewidth = 2)
        plt.title('SVR')
        plt.subplot(122)
        plt.plot(pv_gp,test_tag, 'r*')
        plt.plot([1,max_total_gp],[1,max_total_gp],'b-',linewidth = 2)
        plt.title('GP')

        pylab.figure(2)
        plt.hold('on')
        plt.hist(diff_gp,bins = step,label = 'GP histogram',color = 'red')
        plt.legend()
        plt.hist(diff_svr,bins = step,label = 'SVR histogram',color = 'green')
        plt.legend()
        pylab.show()

