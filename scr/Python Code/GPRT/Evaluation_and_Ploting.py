__author__ = 'Administrator'
from scipy.stats.stats import pearsonr
import calc_mtw
import matplotlib.pyplot as plt
import pylab

def evaluation(pv_gp,test_tag,t_gp):
        step = 10
        #diff_svr = abs(pv_svr - test_tag_svr)
        diff_gp = abs(pv_gp -test_tag)
        #histo_svr = plt.hist(diff_svr,step)
        histo_gp = plt.hist(diff_gp,step)

        max_total_gp = max(test_tag)-min(test_tag)
        #max_total_svr = max(test_tag)

        mtw_gp =calc_mtw.mini_time_win(histo_gp[0],diff_gp,step,max_total_gp)
        corrcoef_gp = pearsonr(pv_gp,test_tag)
        #mtw_svr = calc_mtw.mini_time_win(histo_svr[0],diff_svr,step,max_total_svr)
        #corrcoef_svr = pearsonr(pv_svr,test_tag_svr)

        print 'corrcoef of GP = ',corrcoef_gp[0]
        print 'mtw of GP = ',mtw_gp
        print 'running time of GP = ',t_gp , 'sec'
        return diff_gp,step ,max_total_gp

def ploting(pv_gp,test_tag,max_total_gp,diff_gp,step):
        pylab.figure(1)
        plt.subplot(121)
        #plt.plot(pv_svr,test_tag, 'g*')
        #plt.plot([1,max_total_svr],[1,max_total_svr],'y-',linewidth = 2)
        #plt.title('SVR')
        #plt.subplot(122)
        plt.plot(pv_gp,test_tag, 'g*')
        plt.plot([1,max_total_gp],[1,max_total_gp],'b-',linewidth = 2)
        plt.title('GP')

        #pylab.figure(2)
        #plt.hold('on')
        plt.subplot(122)
        plt.hist(diff_gp,bins = step,label = 'GP histogram',color = 'red')
        plt.legend()
        #plt.subplot(122)
        #line2 = plt.hist(diff_svr,bins = step,label = 'SVR histogram',color = 'green')
        #plt.legend()
        pylab.show()



