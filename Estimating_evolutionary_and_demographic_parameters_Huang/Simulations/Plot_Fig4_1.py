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


plt.plot([0,200000],[20000,20000],'k-.',linewidth = 2)


sam_size = 40

name = "xx_true_list="+str(sam_size)+"_constne_len108.dat"
file = open(name,"rb")
xx= pickle.load(file)
file.close()

xx=xx[0]

name = "nt_true_list="+str(sam_size)+"_constne_len108.dat"
file = open(name,"rb")
yylist = pickle.load(file)
file.close()

yylist = yylist[0:25]

for i in range(len(yylist)):
    for j in range(len(yylist[0])):
        if abs(yylist[i][j]) > 10**10:
            yylist[i][j] = 10**10
        


for i in range(len(yylist)):
    plt.plot(xx,yylist[i],alpha =1)
 
 

################################

# xx = numpy.linspace(0,100000,100)
# yy = 20000*numpy.exp(-0.00001*xx)
# plt.plot(xx,yy,'k-.')

plt.plot([0,200000],[20000,20000],'k-.',linewidth = 2)

plt.legend(['True'],fontsize=16,loc = 'upper left')



 

plt.ylim(0,40000)
plt.xlim(0,70000)

plt.xlabel('time (generations)',fontsize=18)
plt.ylabel('population size',fontsize=18)
plt.title('Estimates of Model C',fontsize=18)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
