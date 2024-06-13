import statistics
import msprime
from IPython.display import SVG, display
import tsinfer
import scipy
import math
import numpy
import tskit
import io
import builtins
import sys
from tqdm.notebook import tqdm
from tskit import MISSING_DATA
import pickle
import random
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import cvxopt as opt
from cvxopt import matrix, spmatrix, sparse
from cvxopt.solvers import qp, options
from cvxopt import blas
from cvxopt import spdiag
import statsmodels
from statsmodels.nonparametric.kernel_regression import KernelReg
from scipy.signal import savgol_filter
import bisect
import time


seq_len = 10**7

########################################################################
################################################################################


muhat = []
muup = []
mulow = []
muhatplt = []






for sam_size in [5,10,20,40]:
    muhatlist = []    
    for byrun in range(1,26):
        
    
        name1 = "fig6_musimall_seq_len="+str(seq_len)+"run="+str(byrun) +"_sam_size="+str(sam_size)+"_noeps=" +"_noconver_constne_efibd.dat"
        name3 = "fig6_sta1simall_seq_len="+str(seq_len)+"run="+str(byrun) +"_sam_size="+str(sam_size)+"_noeps(est)=" +"_noconver_constne_efibd.dat"
        name5 = "fig6_sta1inf_seq_len="+str(seq_len)+"run="+str(byrun) +"_sam_size="+str(sam_size)+"_noeps(est)=" +"_noconver_constne_efibd.dat"
        name7 = "fig6_abdlen_seq_len="+str(seq_len)+"run="+str(byrun) +"_sam_size="+str(sam_size)+"_noeps(est)=" +"_noconver_constne_efibd.dat"



        file = open(name1,"rb") 
        muacceptlist = pickle.load(file)
        file.close()



        file = open(name3,"rb") 
        summary_sta1_simlist = pickle.load(file)
        file.close()


        file = open(name5,"rb") 
        summary_sta1_inf = pickle.load(file)
        file.close()


        file = open(name7,"rb") 
        lllabd = pickle.load(file)
        file.close()
        
        muacceptlist = numpy.array(muacceptlist)*10**8


#############   

### linear regression

        yreg = muacceptlist
        xreg = numpy.zeros( (len(summary_sta1_simlist),1) )
        for i in range(len(summary_sta1_simlist)):
            xreg[i,:] = [summary_sta1_simlist[i] - summary_sta1_inf ]
        reg = LinearRegression().fit(xreg, yreg)

        beta =  reg.coef_
 
        for i in range(len(muacceptlist)):
            muacceptlist[i] = muacceptlist[i] - (xreg[i][0]*beta[0]  )

####

        
        judgevalue = []
        xreg = numpy.zeros( (len(summary_sta1_simlist),1) )
        for i in range(len(summary_sta1_simlist)):
            xreg[i,:] = [summary_sta1_simlist[i] - summary_sta1_inf]
        smat = numpy.cov(numpy.transpose( xreg) )
        for i in range(len(muacceptlist)):
            ss = numpy.array( [summary_sta1_simlist[i] - summary_sta1_inf ])
            aaa = ss/smat
            bbb = aaa*ss
            judgevalue.append( numpy.sqrt(bbb[0]) )
#################
                    
        
        
        muuse = []
        epsuse =[]
        
        
        sind = numpy.argsort(judgevalue)[0:125]
        muuse = numpy.array(muacceptlist)[sind]
         
 

        muhatlist.append( statistics.mean(muuse))
    
    muhat = statistics.mean(muhatlist)
    muhatciup = muhat + 1.96*statistics.stdev(muhatlist) 
    muhatcilow = muhat - 1.96*statistics.stdev(muhatlist) 

    muhatplt.append(muhat)
    muup.append(muhatciup)    
    mulow.append(muhatcilow)
    

x1 = [2-0.2, 4-0.2,6-0.2,8-0.2]
plt.scatter(x1,muhatplt,color='r',label="Model C")




for i in range(len(muhatplt)):
    plt.plot([x1[i],x1[i]], [ mulow[i],muup[i] ], 'r')
    plt.plot( [x1[i]-0.1, x1[i]+0.1] , [mulow[i],mulow[i] ],'r' )
    plt.plot( [x1[i]-0.1, x1[i]+0.1] , [muup[i],muup[i] ],'r' )
#######################################################
muhat = []
muup = []
mulow = []
muhatplt = []

for sam_size in [5,10,20,40]:
    muhatlist = []
    for byrun in range(1,26):
        name1 = "fig6_musimall_seq_len="+str(seq_len)+"run="+str(byrun) +"_sam_size="+str(sam_size)+"_noeps=" +"_noconver_expne_efibd.dat"
        name3 = "fig6_sta1simall_seq_len="+str(seq_len)+"run="+str(byrun) +"_sam_size="+str(sam_size)+"_noeps(est)=" +"_noconver_expne_efibd.dat"
        name5 = "fig6_sta1inf_seq_len="+str(seq_len)+"run="+str(byrun) +"_sam_size="+str(sam_size)+"_noeps(est)=" +"_noconver_expne_efibd.dat"
        name7 = "fig6_abdlen_seq_len="+str(seq_len)+"run="+str(byrun) +"_sam_size="+str(sam_size)+"_noeps(est)=" +"_noconver_expne_efibd.dat"



        file = open(name1,"rb") 
        muacceptlist = pickle.load(file)
        file.close()



        file = open(name3,"rb") 
        summary_sta1_simlist = pickle.load(file)
        file.close()


        file = open(name5,"rb") 
        summary_sta1_inf = pickle.load(file)
        file.close()


        file = open(name7,"rb") 
        lllabd = pickle.load(file)
        file.close()
        
        muacceptlist = numpy.array(muacceptlist)*10**8

#############   

### linear regression

        yreg = muacceptlist
        xreg = numpy.zeros( (len(summary_sta1_simlist),1) )
        for i in range(len(summary_sta1_simlist)):
            xreg[i,:] = [summary_sta1_simlist[i] - summary_sta1_inf ]
        reg = LinearRegression().fit(xreg, yreg)

        beta =  reg.coef_
 
        for i in range(len(muacceptlist)):
            muacceptlist[i] = muacceptlist[i] - (xreg[i][0]*beta[0]  )

####

        
        judgevalue = []
        xreg = numpy.zeros( (len(summary_sta1_simlist),1) )
        for i in range(len(summary_sta1_simlist)):
            xreg[i,:] = [summary_sta1_simlist[i] - summary_sta1_inf]
        smat = numpy.cov(numpy.transpose( xreg) )
        for i in range(len(muacceptlist)):
            ss = numpy.array( [summary_sta1_simlist[i] - summary_sta1_inf ])
            aaa = ss/smat
            bbb = aaa*ss
            judgevalue.append( numpy.sqrt(bbb[0]) )
#################
                    
        
        
        muuse = []
        epsuse =[]
        
        
        sind = numpy.argsort(judgevalue)[0:125]
        muuse = numpy.array(muacceptlist)[sind]
         
 

        muhatlist.append( statistics.mean(muuse))
    
    muhat = statistics.mean(muhatlist)
    muhatciup = muhat + 1.96*statistics.stdev(muhatlist) 
    muhatcilow = muhat - 1.96*statistics.stdev(muhatlist) 

    muhatplt.append(muhat)
    muup.append(muhatciup)    
    mulow.append(muhatcilow)

x1 = [ 2+0.2,4+0.2,6+0.2,8+0.2]
plt.scatter(x1,muhatplt,color='b',marker="x",label="Model Gb")

 

for i in range(len(muhatplt)):
    plt.plot([x1[i],x1[i]], [ mulow[i],muup[i] ], 'b')
    plt.plot( [x1[i]-0.1, x1[i]+0.1] , [mulow[i],mulow[i] ],'b' )
    plt.plot( [x1[i]-0.1, x1[i]+0.1] , [muup[i],muup[i] ],'b' )
    




plt.plot([-1,14],[1.3,1.3],'k-.',label="True")
    
    
plt.legend(fontsize=14,loc = 'upper left')    

 
plt.xlabel('sample size',fontsize=16)
plt.ylabel('estimate of ' + r'$\mu$' + ' ' +r'($10^{-8}$)',fontsize=16)
plt.title('Estimates using inferred TS '+r'($\ell$'+'=' +r'$10^7$)',fontsize=16)
plt.xticks([2,4,6,8],['10','20','40','80'],fontsize=12)
plt.yticks(fontsize=12)
plt.ylim(1.1,1.6)
plt.xlim(0.4,9)
