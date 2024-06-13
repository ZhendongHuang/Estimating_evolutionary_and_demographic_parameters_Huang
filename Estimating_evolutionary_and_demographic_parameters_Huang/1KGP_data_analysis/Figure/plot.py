import argparse
import os
import re
import numpy as np
import stdpopsim

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
import time
import bisect




curpoplist = []
popsimlist = []
ratesimlist = []
s1list = []
s2list = []
s3list = []
s4list = []
for byrun in range(1,26):
    
    name1 = "1M2_ne_realpopsimall_"+"run="+str(byrun) +"BEB_efibd.dat"

    name3 = "1M2_ne_realsta1simall_"+"run="+str(byrun) +"BEB_efibd.dat"
    name4 = "1M2_ne_realsta2simall_"+"run="+str(byrun) +"BEB_efibd.dat"
    name5 = "1M2_ne_realsta1inf_"+"run="+str(byrun) +"BEB_efibd.dat"
    name6 = "1M2_ne_realsta2inf_"+"run="+str(byrun) +"BEB_efibd.dat"
    name7 = "1M2_ne_realsta3simall_"+"run="+str(byrun) +"BEB_efibd.dat"
    name8 = "1M2_ne_realsta4simall_"+"run="+str(byrun) +"BEB_efibd.dat"
    name9 = "1M2_ne_realsta3inf_"+"run="+str(byrun) +"BEB_efibd.dat"
    name10 = "1M2_ne_realsta4inf_"+"run="+str(byrun) +"BEB_efibd.dat"



    file = open(name1,"rb") 
    inipopsimlist= pickle.load(file)
    file.close()



    file = open(name3,"rb") 
    summary_sta1_simlist = pickle.load(file)
    file.close()

    file = open(name4,"rb") 
    summary_sta2_simlist = pickle.load(file)
    file.close()

    file = open(name5,"rb") 
    summary_sta1_inf = pickle.load(file)
    file.close()

    file = open(name6,"rb") 
    summary_sta2_inf = pickle.load(file)
    file.close()

    file = open(name7,"rb") 
    summary_sta3_simlist = pickle.load(file)
    file.close()

    file = open(name8,"rb") 
    summary_sta4_simlist = pickle.load(file)
    file.close()

    file = open(name9,"rb") 
    summary_sta3_inf = pickle.load(file)
    file.close()

    file = open(name10,"rb") 
    summary_sta4_inf = pickle.load(file)
    file.close()
    

    popsimlist = popsimlist + inipopsimlist

    s1list = s1list + summary_sta1_simlist
    s2list = s2list + summary_sta2_simlist
    s1inf = summary_sta1_inf
    s2inf = summary_sta2_inf
    s3list = s3list + summary_sta3_simlist
    s4list = s4list + summary_sta4_simlist
    s3inf = summary_sta3_inf
    s4inf = summary_sta4_inf


# ############   


s1sd = statistics.stdev(s1list)
s3sd = statistics.stdev(s3list)
judgevalue = []
for i in range(len(s1list)):
    judgevalue.append(  abs(s1list[i] - s1inf)/s1inf + abs(s3list[i]- s3inf)/s3inf        )
#################






sind = numpy.argsort(judgevalue)[0:125]

inipaccept = numpy.array(popsimlist)[sind]


popsim = statistics.mean(inipaccept)



genpoint = [0*5/5.2,400*5/5.2,800*5/5.2,1200*5/5.2,1600*5/5.2,2000*5/5.2,2400*5/5.2,2800*5/5.2,\
            3200*5/5.2,3600*5/5.2,4000*5/5.2,8000*5/5.2,12000*5/5.2,16000*5/5.2,20000*5/5.2]
inipopsize =[109000,35500,11500,1900,2300,2900,3500,4200,5800,7700,9600,25800,25800,21200,17300]
genuse = []
puse = []
gruse = []
for i in range(len(genpoint)-1):
    addplist = []
    for j in range(10):
        addplist.append(genpoint[i] + (genpoint[i+1]-genpoint[i])/10*j)
    genuse = genuse + addplist               
for i in range(len(inipopsize)-1):
    addplist = []
    for j in range(10):
        addplist.append(inipopsize[i] + (inipopsize[i+1]-inipopsize[i])/10*j)
    puse = puse + addplist
genuse.append(genpoint[len(genpoint)-1])
puse.append(inipopsize[len(inipopsize)-1])


ind1000 = 0
for i in range(len(genuse)):
    if genuse[i]<1000:
        ind1000 = i
        
ratesim = -numpy.log(puse[ind1000]/popsim)/genuse[ind1000]
        
consta = puse[ind1000] - popsim*numpy.exp(-ratesim*1000)        
for i in range(len(genuse)):
    if genuse[i]<1000:
        puse[i] = popsim*numpy.exp(-ratesim*genuse[i]) 

plt.plot(genuse,numpy.array(puse)*2/10**3,label="BEB")



# ########################################################################################

print(puse[0])








########################################################################################


curpoplist = []
popsimlist = []
ratesimlist = []
s1list = []
s2list = []
s3list = []
s4list = []
for byrun in range(1,26):
    
    name1 = "1M2_ne_realpopsimall_"+"run="+str(byrun) +"CHB_efibd.dat"

    name3 = "1M2_ne_realsta1simall_"+"run="+str(byrun) +"CHB_efibd.dat"
    name4 = "1M2_ne_realsta2simall_"+"run="+str(byrun) +"CHB_efibd.dat"
    name5 = "1M2_ne_realsta1inf_"+"run="+str(byrun) +"CHB_efibd.dat"
    name6 = "1M2_ne_realsta2inf_"+"run="+str(byrun) +"CHB_efibd.dat"
    name7 = "1M2_ne_realsta3simall_"+"run="+str(byrun) +"CHB_efibd.dat"
    name8 = "1M2_ne_realsta4simall_"+"run="+str(byrun) +"CHB_efibd.dat"
    name9 = "1M2_ne_realsta3inf_"+"run="+str(byrun) +"CHB_efibd.dat"
    name10 = "1M2_ne_realsta4inf_"+"run="+str(byrun) +"CHB_efibd.dat"



    file = open(name1,"rb") 
    inipopsimlist= pickle.load(file)
    file.close()


    file = open(name3,"rb") 
    summary_sta1_simlist = pickle.load(file)
    file.close()

    file = open(name4,"rb") 
    summary_sta2_simlist = pickle.load(file)
    file.close()

    file = open(name5,"rb") 
    summary_sta1_inf = pickle.load(file)
    file.close()

    file = open(name6,"rb") 
    summary_sta2_inf = pickle.load(file)
    file.close()

    file = open(name7,"rb") 
    summary_sta3_simlist = pickle.load(file)
    file.close()

    file = open(name8,"rb") 
    summary_sta4_simlist = pickle.load(file)
    file.close()

    file = open(name9,"rb") 
    summary_sta3_inf = pickle.load(file)
    file.close()

    file = open(name10,"rb") 
    summary_sta4_inf = pickle.load(file)
    file.close()
    

    popsimlist = popsimlist + inipopsimlist

    s1list = s1list + summary_sta1_simlist
    s2list = s2list + summary_sta2_simlist
    s1inf = summary_sta1_inf
    s2inf = summary_sta2_inf
    s3list = s3list + summary_sta3_simlist
    s4list = s4list + summary_sta4_simlist
    s3inf = summary_sta3_inf
    s4inf = summary_sta4_inf


#############   



s1sd = statistics.stdev(s1list)
s3sd = statistics.stdev(s3list)
judgevalue = []
for i in range(len(s1list)):
    judgevalue.append(   abs(s1list[i] - s1inf)/s1inf +  abs(s3list[i]- s3inf)/s3inf        )
#################






sind = numpy.argsort(judgevalue)[0:125]

inipaccept = numpy.array(popsimlist)[sind]


popsim = statistics.mean(inipaccept)



genpoint = [0*5/5.2,400*5/5.2,800*5/5.2,1200*5/5.2,1600*5/5.2,2000*5/5.2,2400*5/5.2,2800*5/5.2,\
            3200*5/5.2,3600*5/5.2,4000*5/5.2,8000*5/5.2,12000*5/5.2,16000*5/5.2,20000*5/5.2]
inipopsize =[96800, 22000, 5000, 1540, 2300, 2700,3500,\
              4400,  6250 , 8000, 9600,26000,26500,21200,17300]
genuse = []
puse = []
gruse = []
for i in range(len(genpoint)-1):
    addplist = []
    for j in range(10):
        addplist.append(genpoint[i] + (genpoint[i+1]-genpoint[i])/10*j)
    genuse = genuse + addplist               
for i in range(len(inipopsize)-1):
    addplist = []
    for j in range(10):
        addplist.append(inipopsize[i] + (inipopsize[i+1]-inipopsize[i])/10*j)
    puse = puse + addplist
genuse.append(genpoint[len(genpoint)-1])
puse.append(inipopsize[len(inipopsize)-1])


ind1000 = 0
for i in range(len(genuse)):
    if genuse[i]<1000:
        ind1000 = i
        
ratesim = -numpy.log(puse[ind1000]/popsim)/genuse[ind1000]
        
consta = puse[ind1000] - popsim*numpy.exp(-ratesim*1000)        
for i in range(len(genuse)):
    if genuse[i]<1000:
        puse[i] = popsim*numpy.exp(-ratesim*genuse[i]) 

plt.plot(genuse,numpy.array(puse)*2/10**3,label="CHB")

print(puse[0])

########################################################################################


curpoplist = []
popsimlist = []
ratesimlist = []
s1list = []
s2list = []
s3list = []
s4list = []
for byrun in range(1,26):
    

    name1 = "1M2_ne_realpopsimall_"+"run="+str(byrun) +"ITU_efibd.dat"

    name3 = "1M2_ne_realsta1simall_"+"run="+str(byrun) +"ITU_efibd.dat"
    name4 = "1M2_ne_realsta2simall_"+"run="+str(byrun) +"ITU_efibd.dat"
    name5 = "1M2_ne_realsta1inf_"+"run="+str(byrun) +"ITU_efibd.dat"
    name6 = "1M2_ne_realsta2inf_"+"run="+str(byrun) +"ITU_efibd.dat"
    name7 = "1M2_ne_realsta3simall_"+"run="+str(byrun) +"ITU_efibd.dat"
    name8 = "1M2_ne_realsta4simall_"+"run="+str(byrun) +"ITU_efibd.dat"
    name9 = "1M2_ne_realsta3inf_"+"run="+str(byrun) +"ITU_efibd.dat"
    name10 = "1M2_ne_realsta4inf_"+"run="+str(byrun) +"ITU_efibd.dat"



    file = open(name1,"rb") 
    inipopsimlist= pickle.load(file)
    file.close()



    file = open(name3,"rb") 
    summary_sta1_simlist = pickle.load(file)
    file.close()

    file = open(name4,"rb") 
    summary_sta2_simlist = pickle.load(file)
    file.close()

    file = open(name5,"rb") 
    summary_sta1_inf = pickle.load(file)
    file.close()

    file = open(name6,"rb") 
    summary_sta2_inf = pickle.load(file)
    file.close()

    file = open(name7,"rb") 
    summary_sta3_simlist = pickle.load(file)
    file.close()

    file = open(name8,"rb") 
    summary_sta4_simlist = pickle.load(file)
    file.close()

    file = open(name9,"rb") 
    summary_sta3_inf = pickle.load(file)
    file.close()

    file = open(name10,"rb") 
    summary_sta4_inf = pickle.load(file)
    file.close()
    

    popsimlist = popsimlist + inipopsimlist

    s1list = s1list + summary_sta1_simlist
    s2list = s2list + summary_sta2_simlist
    s1inf = summary_sta1_inf
    s2inf = summary_sta2_inf
    s3list = s3list + summary_sta3_simlist
    s4list = s4list + summary_sta4_simlist
    s3inf = summary_sta3_inf
    s4inf = summary_sta4_inf


#############   

# judgevalue = []
# xreg = numpy.zeros( (len(s1list),2) )
# for i in range(len(s1list)):
#     xreg[i,:] = [s1list[i] - s1inf , s3list[i]- s3inf]
# smat = numpy.cov(numpy.transpose( xreg) )
# for i in range(len(s1list)):
#     ss = numpy.array( [   s1list[i] - s1inf , s3list[i] - s3inf])
#     aaa = numpy.matmul(ss,numpy.linalg.inv(smat))
#     bbb = numpy.matmul(aaa,numpy.transpose(ss))
#     judgevalue.append( numpy.sqrt(bbb) )
# #################

s1sd = statistics.stdev(s1list)
s3sd = statistics.stdev(s3list)
judgevalue = []
for i in range(len(s1list)):
    judgevalue.append(  abs(s3list[i] - s3inf)/s3inf + abs(s4list[i]- s4inf)/s4inf        )
#################






sind = numpy.argsort(judgevalue)[0:125]

inipaccept = numpy.array(popsimlist)[sind]


popsim = statistics.mean(inipaccept)


genpoint = [0*5/5.2,400*5/5.2,800*5/5.2,1200*5/5.2,1600*5/5.2,2000*5/5.2,2400*5/5.2,2800*5/5.2,\
            3200*5/5.2,3600*5/5.2,4000*5/5.2,8000*5/5.2,12000*5/5.2,16000*5/5.2,20000*5/5.2]
inipopsize =[45137, 20200, 9040, 3500, 2700, 3080,3850,\
              4600,  6000 , 8000, 10000,26000,26500,21200,17300]
genuse = []
puse = []
gruse = []
for i in range(len(genpoint)-1):
    addplist = []
    for j in range(10):
        addplist.append(genpoint[i] + (genpoint[i+1]-genpoint[i])/10*j)
    genuse = genuse + addplist               
for i in range(len(inipopsize)-1):
    addplist = []
    for j in range(10):
        addplist.append(inipopsize[i] + (inipopsize[i+1]-inipopsize[i])/10*j)
    puse = puse + addplist
genuse.append(genpoint[len(genpoint)-1])
puse.append(inipopsize[len(inipopsize)-1])


ind1000 = 0
for i in range(len(genuse)):
    if genuse[i]<1000:
        ind1000 = i
        
ratesim = -numpy.log(puse[ind1000]/popsim)/genuse[ind1000]
        
consta = puse[ind1000] - popsim*numpy.exp(-ratesim*1000)        
for i in range(len(genuse)):
    if genuse[i]<1000:
        puse[i] = popsim*numpy.exp(-ratesim*genuse[i]) 

plt.plot(genuse,numpy.array(puse)*2/10**3,label="ITU")

print(puse[0])


# ########################################################################################


curpoplist = []
popsimlist = []
ratesimlist = []
s1list = []
s2list = []
s3list = []
s4list = []
for byrun in range(1,26):
    

    name1 = "1M2_ne_realpopsimall_"+"run="+str(byrun) +"JPT_efibd.dat"

    name3 = "1M2_ne_realsta1simall_"+"run="+str(byrun) +"JPT_efibd.dat"
    name4 = "1M2_ne_realsta2simall_"+"run="+str(byrun) +"JPT_efibd.dat"
    name5 = "1M2_ne_realsta1inf_"+"run="+str(byrun) +"JPT_efibd.dat"
    name6 = "1M2_ne_realsta2inf_"+"run="+str(byrun) +"JPT_efibd.dat"
    name7 = "1M2_ne_realsta3simall_"+"run="+str(byrun) +"JPT_efibd.dat"
    name8 = "1M2_ne_realsta4simall_"+"run="+str(byrun) +"JPT_efibd.dat"
    name9 = "1M2_ne_realsta3inf_"+"run="+str(byrun) +"JPT_efibd.dat"
    name10 = "1M2_ne_realsta4inf_"+"run="+str(byrun) +"JPT_efibd.dat"



    file = open(name1,"rb") 
    inipopsimlist= pickle.load(file)
    file.close()


    file = open(name3,"rb") 
    summary_sta1_simlist = pickle.load(file)
    file.close()

    file = open(name4,"rb") 
    summary_sta2_simlist = pickle.load(file)
    file.close()

    file = open(name5,"rb") 
    summary_sta1_inf = pickle.load(file)
    file.close()

    file = open(name6,"rb") 
    summary_sta2_inf = pickle.load(file)
    file.close()

    file = open(name7,"rb") 
    summary_sta3_simlist = pickle.load(file)
    file.close()

    file = open(name8,"rb") 
    summary_sta4_simlist = pickle.load(file)
    file.close()

    file = open(name9,"rb") 
    summary_sta3_inf = pickle.load(file)
    file.close()

    file = open(name10,"rb") 
    summary_sta4_inf = pickle.load(file)
    file.close()
    

    popsimlist = popsimlist + inipopsimlist

    s1list = s1list + summary_sta1_simlist
    s2list = s2list + summary_sta2_simlist
    s1inf = summary_sta1_inf
    s2inf = summary_sta2_inf
    s3list = s3list + summary_sta3_simlist
    s4list = s4list + summary_sta4_simlist
    s3inf = summary_sta3_inf
    s4inf = summary_sta4_inf


#############   
#############   

# judgevalue = []
# xreg = numpy.zeros( (len(s1list),2) )
# for i in range(len(s1list)):
#     xreg[i,:] = [s1list[i] - s1inf , s3list[i]- s3inf]
# smat = numpy.cov(numpy.transpose( xreg) )
# for i in range(len(s1list)):
#     ss = numpy.array( [ s1list[i] - s1inf , s3list[i] - s3inf])
#     aaa = numpy.matmul(ss,numpy.linalg.inv(smat))
#     bbb = numpy.matmul(aaa,numpy.transpose(ss))
#     judgevalue.append( numpy.sqrt(bbb) )
# #################

s1sd = statistics.stdev(s1list)
s3sd = statistics.stdev(s3list)
judgevalue = []
for i in range(len(s1list)):
    judgevalue.append(  abs(s3list[i] - s3inf)/s3inf + abs(s4list[i]- s4inf)/s4inf        )
#################



sind = numpy.argsort(judgevalue)[0:125]

inipaccept = numpy.array(popsimlist)[sind]


popsim = statistics.mean(inipaccept)



genpoint = [0*5/5.2,400*5/5.2,800*5/5.2,1200*5/5.2,1600*5/5.2,2000*5/5.2,2400*5/5.2,2800*5/5.2,\
            3200*5/5.2,3600*5/5.2,4000*5/5.2,8000*5/5.2,12000*5/5.2,16000*5/5.2,20000*5/5.2]
inipopsize =[76232*5/5.2, 17400*5/5.2, 4000*5/5.2, 1400*5/5.2, 2200*5/5.2, 2800*5/5.2,3600*5/5.2,\
             4400*5/5.2,6000*5/5.2,8000*5/5.2,10000*5/5.2,26800*5/5.2,26800*5/5.2,22000*5/5.2,18000*5/5.2]
genuse = []
puse = []
gruse = []
for i in range(len(genpoint)-1):
    addplist = []
    for j in range(10):
        addplist.append(genpoint[i] + (genpoint[i+1]-genpoint[i])/10*j)
    genuse = genuse + addplist               
for i in range(len(inipopsize)-1):
    addplist = []
    for j in range(10):
        addplist.append(inipopsize[i] + (inipopsize[i+1]-inipopsize[i])/10*j)
    puse = puse + addplist
genuse.append(genpoint[len(genpoint)-1])
puse.append(inipopsize[len(inipopsize)-1])


ind1000 = 0
for i in range(len(genuse)):
    if genuse[i]<1000:
        ind1000 = i
        
ratesim = -numpy.log(puse[ind1000]/popsim)/genuse[ind1000]
        
consta = puse[ind1000] - popsim*numpy.exp(-ratesim*1000)        
for i in range(len(genuse)):
    if genuse[i]<1000:
        puse[i] = popsim*numpy.exp(-ratesim*genuse[i]) 

plt.plot(genuse,numpy.array(puse)*2/10**3,label="JPT")

print(puse[0])

########################################################################################


# ########################################################################################

########################################################################################


curpoplist = []
popsimlist = []
ratesimlist = []
s1list = []
s2list = []
s3list = []
s4list = []
for byrun in range(1,26):
    

    name1 = "1M2_ne_realpopsimall_"+"run="+str(byrun) +"GBR_efibd.dat"

    name3 = "1M2_ne_realsta1simall_"+"run="+str(byrun) +"GBR_efibd.dat"
    name4 = "1M2_ne_realsta2simall_"+"run="+str(byrun) +"GBR_efibd.dat"
    name5 = "1M2_ne_realsta1inf_"+"run="+str(byrun) +"GBR_efibd.dat"
    name6 = "1M2_ne_realsta2inf_"+"run="+str(byrun) +"GBR_efibd.dat"
    name7 = "1M2_ne_realsta3simall_"+"run="+str(byrun) +"GBR_efibd.dat"
    name8 = "1M2_ne_realsta4simall_"+"run="+str(byrun) +"GBR_efibd.dat"
    name9 = "1M2_ne_realsta3inf_"+"run="+str(byrun) +"GBR_efibd.dat"
    name10 = "1M2_ne_realsta4inf_"+"run="+str(byrun) +"GBR_efibd.dat"



    file = open(name1,"rb") 
    inipopsimlist= pickle.load(file)
    file.close()



    file = open(name3,"rb") 
    summary_sta1_simlist = pickle.load(file)
    file.close()

    file = open(name4,"rb") 
    summary_sta2_simlist = pickle.load(file)
    file.close()

    file = open(name5,"rb") 
    summary_sta1_inf = pickle.load(file)
    file.close()

    file = open(name6,"rb") 
    summary_sta2_inf = pickle.load(file)
    file.close()

    file = open(name7,"rb") 
    summary_sta3_simlist = pickle.load(file)
    file.close()

    file = open(name8,"rb") 
    summary_sta4_simlist = pickle.load(file)
    file.close()

    file = open(name9,"rb") 
    summary_sta3_inf = pickle.load(file)
    file.close()

    file = open(name10,"rb") 
    summary_sta4_inf = pickle.load(file)
    file.close()
    

    popsimlist = popsimlist + inipopsimlist

    s1list = s1list + summary_sta1_simlist
    s2list = s2list + summary_sta2_simlist
    s1inf = summary_sta1_inf
    s2inf = summary_sta2_inf
    s3list = s3list + summary_sta3_simlist
    s4list = s4list + summary_sta4_simlist
    s3inf = summary_sta3_inf
    s4inf = summary_sta4_inf


#############   


s1sd = statistics.stdev(s1list)
s3sd = statistics.stdev(s3list)
judgevalue = []
for i in range(len(s1list)):
    judgevalue.append(  abs(s1list[i] - s1inf)/s1inf +  + abs(s4list[i]- s4inf)/s4inf        )
#################






sind = numpy.argsort(judgevalue)[0:125]

inipaccept = numpy.array(popsimlist)[sind]


popsim = statistics.mean(inipaccept)



genpoint = [0*5/5.2,400*5/5.2,800*5/5.2,1200*5/5.2,1600*5/5.2,2000*5/5.2,2400*5/5.2,2800*5/5.2,\
            3200*5/5.2,3600*5/5.2,4000*5/5.2,8000*5/5.2,12000*5/5.2,16000*5/5.2,20000*5/5.2]
inipopsize =[30112*5/5.2, 13000*5/5.2, 5600*5/5.2, 2600*5/5.2, 2600*5/5.2, 2800*5/5.2,3600*5/5.2,\
             4400*5/5.2,6000*5/5.2,8000*5/5.2,10000*5/5.2,26800*5/5.2,26800*5/5.2,22000*5/5.2,18000*5/5.2]
genuse = []
puse = []
gruse = []
for i in range(len(genpoint)-1):
    addplist = []
    for j in range(10):
        addplist.append(genpoint[i] + (genpoint[i+1]-genpoint[i])/10*j)
    genuse = genuse + addplist               
for i in range(len(inipopsize)-1):
    addplist = []
    for j in range(10):
        addplist.append(inipopsize[i] + (inipopsize[i+1]-inipopsize[i])/10*j)
    puse = puse + addplist
genuse.append(genpoint[len(genpoint)-1])
puse.append(inipopsize[len(inipopsize)-1])


ind1000 = 0
for i in range(len(genuse)):
    if genuse[i]<1000:
        ind1000 = i
        
ratesim = -numpy.log(puse[ind1000]/popsim)/genuse[ind1000]
        
consta = puse[ind1000] - popsim*numpy.exp(-ratesim*1000)        
for i in range(len(genuse)):
    if genuse[i]<1000:
        puse[i] = popsim*numpy.exp(-ratesim*genuse[i]) 

plt.plot(genuse,numpy.array(puse)*2/10**3,label="GBR")

print(puse[0])

# ########################################################################################


########################################################################################


curpoplist = []
popsimlist = []
ratesimlist = []
s1list = []
s2list = []
s3list = []
s4list = []
for byrun in range(1,26):
    

    name1 = "1M2_ne_realpopsimall_"+"run="+str(byrun) +"FIN_efibd.dat"

    name3 = "1M2_ne_realsta1simall_"+"run="+str(byrun) +"FIN_efibd.dat"
    name4 = "1M2_ne_realsta2simall_"+"run="+str(byrun) +"FIN_efibd.dat"
    name5 = "1M2_ne_realsta1inf_"+"run="+str(byrun) +"FIN_efibd.dat"
    name6 = "1M2_ne_realsta2inf_"+"run="+str(byrun) +"FIN_efibd.dat"
    name7 = "1M2_ne_realsta3simall_"+"run="+str(byrun) +"FIN_efibd.dat"
    name8 = "1M2_ne_realsta4simall_"+"run="+str(byrun) +"FIN_efibd.dat"
    name9 = "1M2_ne_realsta3inf_"+"run="+str(byrun) +"FIN_efibd.dat"
    name10 = "1M2_ne_realsta4inf_"+"run="+str(byrun) +"FIN_efibd.dat"


    file = open(name1,"rb") 
    inipopsimlist= pickle.load(file)
    file.close()



    file = open(name3,"rb") 
    summary_sta1_simlist = pickle.load(file)
    file.close()

    file = open(name4,"rb") 
    summary_sta2_simlist = pickle.load(file)
    file.close()

    file = open(name5,"rb") 
    summary_sta1_inf = pickle.load(file)
    file.close()

    file = open(name6,"rb") 
    summary_sta2_inf = pickle.load(file)
    file.close()

    file = open(name7,"rb") 
    summary_sta3_simlist = pickle.load(file)
    file.close()

    file = open(name8,"rb") 
    summary_sta4_simlist = pickle.load(file)
    file.close()

    file = open(name9,"rb") 
    summary_sta3_inf = pickle.load(file)
    file.close()

    file = open(name10,"rb") 
    summary_sta4_inf = pickle.load(file)
    file.close()
    

    popsimlist = popsimlist + inipopsimlist

    s1list = s1list + summary_sta1_simlist
    s2list = s2list + summary_sta2_simlist
    s1inf = summary_sta1_inf
    s2inf = summary_sta2_inf
    s3list = s3list + summary_sta3_simlist
    s4list = s4list + summary_sta4_simlist
    s3inf = summary_sta3_inf
    s4inf = summary_sta4_inf


#############   



s1sd = statistics.stdev(s1list)
s3sd = statistics.stdev(s3list)
judgevalue = []
for i in range(len(s1list)):
    judgevalue.append(  abs(s1list[i] - s1inf)/s1inf   + abs(s4list[i] - s4inf)/s4inf )

#################






sind = numpy.argsort(judgevalue)[0:125]

inipaccept = numpy.array(popsimlist)[sind]


popsim = statistics.mean(inipaccept)



genpoint = [0*5/5.2,400*5/5.2,800*5/5.2,1200*5/5.2,1600*5/5.2,2000*5/5.2,2400*5/5.2,2800*5/5.2,\
            3200*5/5.2,3600*5/5.2,4000*5/5.2,8000*5/5.2,12000*5/5.2,16000*5/5.2,20000*5/5.2]
inipopsize =[20000*5/5.2, 10000*5/5.2, 5000*5/5.2, 2600*5/5.2, 2600*5/5.2, 2800*5/5.2, 3600*5/5.2,4400*5/5.2,\
             6000*5/5.2,8000*5/5.2,10000*5/5.2,26800*5/5.2,26800*5/5.2,22000*5/5.2,18000*5/5.2]
genuse = []
puse = []
gruse = []
for i in range(len(genpoint)-1):
    addplist = []
    for j in range(10):
        addplist.append(genpoint[i] + (genpoint[i+1]-genpoint[i])/10*j)
    genuse = genuse + addplist               
for i in range(len(inipopsize)-1):
    addplist = []
    for j in range(10):
        addplist.append(inipopsize[i] + (inipopsize[i+1]-inipopsize[i])/10*j)
    puse = puse + addplist
genuse.append(genpoint[len(genpoint)-1])
puse.append(inipopsize[len(inipopsize)-1])


ind1000 = 0
for i in range(len(genuse)):
    if genuse[i]<1000:
        ind1000 = i
        
ratesim = -numpy.log(puse[ind1000]/popsim)/genuse[ind1000]
        
consta = puse[ind1000] - popsim*numpy.exp(-ratesim*1000)        
for i in range(len(genuse)):
    if genuse[i]<1000:
        puse[i] = popsim*numpy.exp(-ratesim*genuse[i]) 

plt.plot(genuse,numpy.array(puse)*2/10**3,label="FIN")

print(puse[0])

##################################################################################

##################################################################################




# ########################################################################################




curpoplist = []
popsimlist = []
ratesimlist = []
s1list = []
s2list = []
s3list = []
s4list = []
for byrun in range(1,26):
    

    name1 = "1M2_ne_realpopsimall_"+"run="+str(byrun) +"MSL_efibd.dat"

    name3 = "1M2_ne_realsta1simall_"+"run="+str(byrun) +"MSL_efibd.dat"
    name4 = "1M2_ne_realsta2simall_"+"run="+str(byrun) +"MSL_efibd.dat"
    name5 = "1M2_ne_realsta1inf_"+"run="+str(byrun) +"MSL_efibd.dat"
    name6 = "1M2_ne_realsta2inf_"+"run="+str(byrun) +"MSL_efibd.dat"
    name7 = "1M2_ne_realsta3simall_"+"run="+str(byrun) +"MSL_efibd.dat"
    name8 = "1M2_ne_realsta4simall_"+"run="+str(byrun) +"MSL_efibd.dat"
    name9 = "1M2_ne_realsta3inf_"+"run="+str(byrun) +"MSL_efibd.dat"
    name10 = "1M2_ne_realsta4inf_"+"run="+str(byrun) +"MSL_efibd.dat"



    file = open(name1,"rb") 
    inipopsimlist= pickle.load(file)
    file.close()



    file = open(name3,"rb") 
    summary_sta1_simlist = pickle.load(file)
    file.close()

    file = open(name4,"rb") 
    summary_sta2_simlist = pickle.load(file)
    file.close()

    file = open(name5,"rb") 
    summary_sta1_inf = pickle.load(file)
    file.close()

    file = open(name6,"rb") 
    summary_sta2_inf = pickle.load(file)
    file.close()

    file = open(name7,"rb") 
    summary_sta3_simlist = pickle.load(file)
    file.close()

    file = open(name8,"rb") 
    summary_sta4_simlist = pickle.load(file)
    file.close()

    file = open(name9,"rb") 
    summary_sta3_inf = pickle.load(file)
    file.close()

    file = open(name10,"rb") 
    summary_sta4_inf = pickle.load(file)
    file.close()
    

    popsimlist = popsimlist + inipopsimlist

    s1list = s1list + summary_sta1_simlist
    s2list = s2list + summary_sta2_simlist
    s1inf = summary_sta1_inf
    s2inf = summary_sta2_inf
    s3list = s3list + summary_sta3_simlist
    s4list = s4list + summary_sta4_simlist
    s3inf = summary_sta3_inf
    s4inf = summary_sta4_inf


# ############   


s1sd = statistics.stdev(s1list)
s3sd = statistics.stdev(s3list)
judgevalue = []
for i in range(len(s1list)):
    judgevalue.append(  abs(s4list[i] - s4inf)/s4inf        )
#################






sind = numpy.argsort(judgevalue)[0:125]

inipaccept = numpy.array(popsimlist)[sind]


popsim = statistics.mean(inipaccept)



genpoint = [0*5/5.2,400*5/5.2,800*5/5.2,1200*5/5.2,1600*5/5.2,2000*5/5.2,2400*5/5.2,2800*5/5.2,\
            3200*5/5.2,3600*5/5.2,4000*5/5.2,8000*5/5.2,12000*5/5.2,16000*5/5.2,20000*5/5.2]
inipopsize =[12500,12500,11900,11500,10600,9600,9000,9600,10600,11900,13500,28500,26900,20800,17300]
genuse = []
puse = []
gruse = []
for i in range(len(genpoint)-1):
    addplist = []
    for j in range(10):
        addplist.append(genpoint[i] + (genpoint[i+1]-genpoint[i])/10*j)
    genuse = genuse + addplist               
for i in range(len(inipopsize)-1):
    addplist = []
    for j in range(10):
        addplist.append(inipopsize[i] + (inipopsize[i+1]-inipopsize[i])/10*j)
    puse = puse + addplist
genuse.append(genpoint[len(genpoint)-1])
puse.append(inipopsize[len(inipopsize)-1])


ind1000 = 0
for i in range(len(genuse)):
    if genuse[i]<1000:
        ind1000 = i
        
ratesim = -numpy.log(puse[ind1000]/popsim)/genuse[ind1000]
        
consta = puse[ind1000] - popsim*numpy.exp(-ratesim*1000)        
for i in range(len(genuse)):
    if genuse[i]<1000:
        puse[i] = popsim*numpy.exp(-ratesim*genuse[i]) 

plt.plot(genuse,numpy.array(puse)*2/10**3,label="MSL")

print(puse[0])

# ########################################################################################

##################################################################################

curpoplist = []
popsimlist = []
ratesimlist = []
s1list = []
s2list = []
s3list = []
s4list = []
for byrun in range(1,26):
    

    name1 = "1M2_ne_realpopsimall_"+"run="+str(byrun) +"LWK_efibd.dat"

    name3 = "1M2_ne_realsta1simall_"+"run="+str(byrun) +"LWK_efibd.dat"
    name4 = "1M2_ne_realsta2simall_"+"run="+str(byrun) +"LWK_efibd.dat"
    name5 = "1M2_ne_realsta1inf_"+"run="+str(byrun) +"LWK_efibd.dat"
    name6 = "1M2_ne_realsta2inf_"+"run="+str(byrun) +"LWK_efibd.dat"
    name7 = "1M2_ne_realsta3simall_"+"run="+str(byrun) +"LWK_efibd.dat"
    name8 = "1M2_ne_realsta4simall_"+"run="+str(byrun) +"LWK_efibd.dat"
    name9 = "1M2_ne_realsta3inf_"+"run="+str(byrun) +"LWK_efibd.dat"
    name10 = "1M2_ne_realsta4inf_"+"run="+str(byrun) +"LWK_efibd.dat"


    file = open(name1,"rb") 
    inipopsimlist= pickle.load(file)
    file.close()


    file = open(name3,"rb") 
    summary_sta1_simlist = pickle.load(file)
    file.close()

    file = open(name4,"rb") 
    summary_sta2_simlist = pickle.load(file)
    file.close()

    file = open(name5,"rb") 
    summary_sta1_inf = pickle.load(file)
    file.close()

    file = open(name6,"rb") 
    summary_sta2_inf = pickle.load(file)
    file.close()

    file = open(name7,"rb") 
    summary_sta3_simlist = pickle.load(file)
    file.close()

    file = open(name8,"rb") 
    summary_sta4_simlist = pickle.load(file)
    file.close()

    file = open(name9,"rb") 
    summary_sta3_inf = pickle.load(file)
    file.close()

    file = open(name10,"rb") 
    summary_sta4_inf = pickle.load(file)
    file.close()
    

    popsimlist = popsimlist + inipopsimlist

    s1list = s1list + summary_sta1_simlist
    s2list = s2list + summary_sta2_simlist
    s1inf = summary_sta1_inf
    s2inf = summary_sta2_inf
    s3list = s3list + summary_sta3_simlist
    s4list = s4list + summary_sta4_simlist
    s3inf = summary_sta3_inf
    s4inf = summary_sta4_inf


# ############   


s1sd = statistics.stdev(s1list)
s3sd = statistics.stdev(s3list)
judgevalue = []
for i in range(len(s1list)):
    judgevalue.append(    abs(s4list[i]- s4inf)/s4inf        )
#################






sind = numpy.argsort(judgevalue)[0:125]

inipaccept = numpy.array(popsimlist)[sind]


popsim = statistics.mean(inipaccept)



genpoint = [0*5/5.2,400*5/5.2,800*5/5.2,1200*5/5.2,1600*5/5.2,2000*5/5.2,2400*5/5.2,2800*5/5.2,\
            3200*5/5.2,3600*5/5.2,4000*5/5.2,8000*5/5.2,12000*5/5.2,16000*5/5.2,20000*5/5.2]
inipopsize =[15200*5/5.2, 15200*5/5.2, 15600*5/5.2, 14800*5/5.2, 12400*5/5.2, 10400*5/5.2,9800*5/5.2,\
             10400*5/5.2,12000*5/5.2,14400*5/5.2,16800*5/5.2,29400*5/5.2,28000*5/5.2,22000*5/5.2,18000*5/5.2]
genuse = []
puse = []
gruse = []
for i in range(len(genpoint)-1):
    addplist = []
    for j in range(10):
        addplist.append(genpoint[i] + (genpoint[i+1]-genpoint[i])/10*j)
    genuse = genuse + addplist               
for i in range(len(inipopsize)-1):
    addplist = []
    for j in range(10):
        addplist.append(inipopsize[i] + (inipopsize[i+1]-inipopsize[i])/10*j)
    puse = puse + addplist
genuse.append(genpoint[len(genpoint)-1])
puse.append(inipopsize[len(inipopsize)-1])


ind1000 = 0
for i in range(len(genuse)):
    if genuse[i]<1000:
        ind1000 = i
        
ratesim = -numpy.log(puse[ind1000]/popsim)/genuse[ind1000]
        
consta = puse[ind1000] - popsim*numpy.exp(-ratesim*1000)        
for i in range(len(genuse)):
    if genuse[i]<1000:
        puse[i] = popsim*numpy.exp(-ratesim*genuse[i]) 

plt.plot(genuse,numpy.array(puse)*2/10**3,label="LWK")
print(puse[0])

##################################################################################





plt.ylim(0,200)
plt.xlabel('generations ago',fontsize = 16)
plt.ylabel(r"population size $(10^3)$",fontsize = 16)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12,loc = 9,ncol=2)
# plt.title("Estimates with Chromosome 20",fontsize = 16)

 
plt.xlim(0,1000)
plt.legend()