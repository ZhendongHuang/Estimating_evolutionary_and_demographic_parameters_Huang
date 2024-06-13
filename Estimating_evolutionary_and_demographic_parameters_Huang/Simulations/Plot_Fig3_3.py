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

plt.figure(figsize=(6.4, 4.8))

muhat = []
muup = []
mulow = []
muhatplt = []

sam_size = 5
for epsilon in [0,0.0002,0.0004,0.0006,0.0008,0.001]:
    name = "muhat_true_" + "eps=" + str(epsilon) + "_sam_size=" + str(sam_size) + "constne.dat"
    file = open(name,"rb")
    mulist = pickle.load(file)
    file.close()
    
    mulist = numpy.array(mulist)*10**8
    muhat.append(statistics.mean(mulist))
    muup.append( muhat[-1] + 1.96*statistics.stdev(mulist) )
    mulow.append( muhat[-1] - 1.96*statistics.stdev(mulist)    )
    muhatplt.append(muhat[-1])


x1 = [0-0.2, 2-0.2, 4-0.2,6-0.2,8-0.2,10-0.2]
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

sam_size = 5
for epsilon in [0,0.0002,0.0004,0.0006,0.0008,0.001]:
    name = "muhat_true_" + "eps=" + str(epsilon) + "_sam_size=" + str(sam_size) + "expne0001.dat"
    file = open(name,"rb")
    mulist = pickle.load(file)
    file.close()
    
    mulist = numpy.array(mulist)*10**8
    muhat.append(statistics.mean(mulist))
    muup.append( muhat[-1] + 1.96*statistics.stdev(mulist)  )
    mulow.append( muhat[-1] - 1.96*statistics.stdev(mulist)    )
    muhatplt.append(muhat[-1])


x1 = [0, 2, 4,6,8,10]
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

sam_size = 5
for epsilon in [0,0.0002,0.0004,0.0006,0.0008,0.001]:
    name = "muhat_true_" + "eps=" + str(epsilon) + "_sam_size=" + str(sam_size) + "bottleneck.dat"
    file = open(name,"rb")
    mulist = pickle.load(file)
    file.close()
    
    mulist = numpy.array(mulist)*10**8
    muhat.append(statistics.mean(mulist))
    muup.append( muhat[-1] + 1.96*statistics.stdev(mulist)  )
    mulow.append( muhat[-1] - 1.96*statistics.stdev(mulist)    )
    muhatplt.append(muhat[-1])


x1 = [0+0.2, 2+0.2, 4+0.2,6+0.2,8+0.2,10+0.2]
plt.scatter(x1,muhatplt,color='g',marker="s",label="Model S")

 

for i in range(len(muhat)):
    plt.plot([x1[i],x1[i]], [ mulow[i],muup[i] ], 'g')
    plt.plot( [x1[i]-0.1, x1[i]+0.1] , [mulow[i],mulow[i] ],'g' )
    plt.plot( [x1[i]-0.1, x1[i]+0.1] , [muup[i],muup[i] ],'g' )
    




plt.plot([-1,14],[1.3,1.3],'k-.',label="True")
    
    
plt.legend(fontsize=14,loc = 'upper left')    
plt.xlabel('sequencing error rate ' + r'$\epsilon$' + ' ' + r'($10^{-4}$)',fontsize=16)
plt.ylabel('estimate of ' + r'$\mu$' + ' ' +r'($10^{-8}$)',fontsize=16)
plt.title(r'$\ell$'+' = ' +r'$10^7$',fontsize=16)
plt.xticks([0,2,4,6,8,10],fontsize=12)
plt.yticks(fontsize=12)
plt.ylim(1.1,1.6)
plt.xlim(-0.4,10.4)
