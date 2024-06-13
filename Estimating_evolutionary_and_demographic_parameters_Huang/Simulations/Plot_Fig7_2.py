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

 

sam_size = 100



xx = numpy.linspace(0,100000,1000)
yy = 20000*numpy.exp(-0*xx)
plt.plot(xx,yy,'k-.')
##############################################################
 

inipopaccept = []
increrateaccept = []
for run in range(1,31):   
    name = "aaacceptlist"+str(31-run)+"pop10000_neconst.dat"  
    file = open(name,"rb")
    acceptlist = pickle.load(file)
    file.close()
    
    inipopacceptcan = []
    increrateacceptcan = []
    for j in range(len(acceptlist)):
        inipopacceptcan.append(acceptlist[j][0])
        increrateacceptcan.append(acceptlist[j][1])
    inipopaccept.append(statistics.mean(inipopacceptcan))
    increrateaccept.append(statistics.mean(increrateacceptcan))
    
    
    
    

    xrange = numpy.linspace(0,200000,2000)
    yrange = 2*inipopaccept[-1]*numpy.exp(-increrateaccept[-1]*xrange)
#     yrangeup = 2*(pophatup)*numpy.exp(-increrhatlow*xrange)
#     yrangelow = 2*(pophatlow)*numpy.exp(-increrhatup*xrange)


    plt.plot(xrange,yrange)
#     plt.fill_between(xrange, yrangelow, yrangeup, color='r', alpha=.3)
    
# pophat = statistics.mean(inipopaccept)
# increrhat = statistics.mean(increrateaccept)
# pophatup = pophat+1.96*statistics.stdev(inipopaccept)
# pophatlow = pophat-1.96*statistics.stdev(inipopaccept) 
# increrhatup = increrhat+1.96*statistics.stdev(increrateaccept) 
# increrhatlow = increrhat-1.96*statistics.stdev(increrateaccept) 


# # pophatup = pophat+1.96*statistics.stdev(inipopaccept) 
# # pophatlow = pophat-1.96*statistics.stdev(inipopaccept) 
# # increrhatup = increrhat+1.96*statistics.stdev(increrateaccept) 
# # increrhatlow = increrhat-1.96*statistics.stdev(increrateaccept) 


# xrange = numpy.linspace(0,200000,2000)
# yrange = 2*pophat*numpy.exp(-increrhat*xrange)
# yrangeup = 2*(pophatup)*numpy.exp(-increrhatlow*xrange)
# yrangelow = 2*(pophatlow)*numpy.exp(-increrhatup*xrange)


# plt.plot(xrange,yrange,'r')
# plt.fill_between(xrange, yrangelow, yrangeup, color='r', alpha=.3)


############################################# 
# name = "xx_inferred_list="+str(sam_size)+"_constne_seq6.dat"
# file = open(name,"rb")
# xx= pickle.load(file)
# file.close()

# xx=xx[0]
 
# name = "nt_inferred_list="+str(sam_size)+"_constne_seq6.dat"
# file = open(name,"rb")
# yylist = pickle.load(file)
# file.close()

# #delete outlier 28
# del yylist[28]

# yylist = numpy.array(yylist) 
# for i in range(len(yylist)):
#     for j in range(len(yylist[0])):
#         if abs(yylist[i,j]) > 10**10:
#             yylist[i,j] = 10**10
        

# yy = []
# yyup = []
# yylow = []
# for i in range(len(yylist[0])):
#     yy.append(statistics.mean(yylist[:,i]))
#     yyup.append( yy[-1] + 1.96*statistics.stdev(yylist[:,i]) /numpy.sqrt(30) ) 
#     yylow.append( yy[-1] - 1.96*statistics.stdev(yylist[:,i]) /numpy.sqrt(30) )

# plt.plot(xx,yy,'b')
# plt.fill_between(xx, yylow, yyup, color='b', alpha=.3)

############################################### 




xx = numpy.linspace(0,100000,1000)
yy = 20000*numpy.exp(-0*xx)
plt.plot(xx,yy,'k-.',label="True", linewidth = 2.5)

plt.legend(fontsize=12,loc = 'upper left')



 

plt.ylim(0,50000)
plt.xlim(0,70000)

plt.xlabel('time (generations)',fontsize=16)
plt.ylabel('population size',fontsize=16)
plt.title('Estimates using \'TSABC\' ' ,fontsize=16)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)


print(statistics.mean(inipopaccept)*2)
print(statistics.stdev(inipopaccept)*numpy.sqrt(2) )

print(statistics.mean(increrateaccept))
print(statistics.stdev(increrateaccept))

