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

muhat = []
muup = []
mulow = []
muhatplt = []

epsilon = 0.001
for sam_size in [5,10,20,40,80 ]:
    name = "epsilonhat_true_" + "eps=" + str(epsilon) + "_sam_size=" + str(sam_size) + "constne_len108.dat"
    file = open(name,"rb")
    mulist = pickle.load(file)
    file.close()
    
    mulist = numpy.array(mulist)*10**3
    muhat.append(statistics.mean(mulist))
    muup.append( muhat[-1] + 1.96*statistics.stdev(mulist) )
    mulow.append( muhat[-1] - 1.96*statistics.stdev(mulist)    )
    muhatplt.append(muhat[-1])


x1 = [2-0.2, 4-0.2,6-0.2,8-0.2,10-0.2 ]
plt.scatter(x1,muhatplt,color='r',label="Model C")




for i in range(len(muhat)):
    plt.plot([x1[i],x1[i]], [ mulow[i],muup[i] ], 'r')
    plt.plot( [x1[i]-0.1, x1[i]+0.1] , [mulow[i],mulow[i] ],'r' )
    plt.plot( [x1[i]-0.1, x1[i]+0.1] , [muup[i],muup[i] ],'r' )
#######################################################
muhat = []
muup = []
mulow = []
muhatplt = []

epsilon = 0.001
for sam_size in [5,10,20,40,80 ]:
    name = "epsilonhat_true_" + "eps=" + str(epsilon) + "_sam_size=" + str(sam_size) + "expne0001_len108.dat"
    file = open(name,"rb")
    mulist = pickle.load(file)
    file.close()
    
    mulist = numpy.array(mulist)*10**3
    muhat.append(statistics.mean(mulist))
    muup.append( muhat[-1] + 1.96*statistics.stdev(mulist) )
    mulow.append( muhat[-1] - 1.96*statistics.stdev(mulist)    )
    muhatplt.append(muhat[-1])
    
x1 = [2, 4,6,8,10 ]
plt.scatter(x1,muhatplt,color='b',marker="x",label="Model Ga")

 

for i in range(len(muhat)):
    plt.plot([x1[i],x1[i]], [ mulow[i],muup[i] ], 'b')
    plt.plot( [x1[i]-0.1, x1[i]+0.1] , [mulow[i],mulow[i] ],'b' )
    plt.plot( [x1[i]-0.1, x1[i]+0.1] , [muup[i],muup[i] ],'b' )
        
#######################################################
muhat = []
muup = []
mulow = []
muhatplt = []

epsilon = 0.001
for sam_size in [5,10,20,40,80 ]:
    name = "epsilonhat_true_" + "eps=" + str(epsilon) + "_sam_size=" + str(sam_size) + "bottleneck_len108.dat"
    file = open(name,"rb")
    mulist = pickle.load(file)
    file.close()
    
    mulist = numpy.array(mulist)*10**3
    muhat.append(statistics.mean(mulist))
    muup.append( muhat[-1] + 1.96*statistics.stdev(mulist) )
    mulow.append( muhat[-1] - 1.96*statistics.stdev(mulist)    )
    muhatplt.append(muhat[-1])

x1 = [ 2+0.2, 4+0.2,6+0.2,8+0.2,10+0.2]
plt.scatter(x1,muhatplt,color='g',marker="s",label="Model S")

 

for i in range(len(muhat)):
    plt.plot([x1[i],x1[i]], [ mulow[i],muup[i] ], 'g')
    plt.plot( [x1[i]-0.1, x1[i]+0.1] , [mulow[i],mulow[i] ],'g' )
    plt.plot( [x1[i]-0.1, x1[i]+0.1] , [muup[i],muup[i] ],'g' )
    




plt.plot([-1,14],[1,1],'k-.',label="True")
    
    
plt.legend(fontsize=14,loc = 'upper left')    
plt.xlabel('sample size',fontsize=16)
plt.ylabel('estimate of ' + r'$\epsilon$' + ' ' +r'($10^{-3}$)',fontsize=16)
plt.title(r'$\ell$'+' = ' +r'$10^8$',fontsize=16)
plt.xticks([2,4,6,8,10],['10','20','40','80','160'],fontsize=12)
plt.yticks(fontsize=12)
plt.ylim(0.95,1.1)
plt.xlim(1,11)