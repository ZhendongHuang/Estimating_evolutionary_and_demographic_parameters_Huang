


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
def find_s(elem, sorted_list):
    'Locate the leftmost value exactly equal to x'
    i = bisect.bisect_left(sorted_list, elem)
    if i != len(sorted_list) and sorted_list[i] == elem:
        return i


def find_s_se(elem, sorted_list):
    'Locate the leftmost value smaller or equal to x'
    i = bisect.bisect_right(sorted_list, elem)
    if i != 0:
        return i-1
    else:
        print("error")


def find_s_s(elem, sorted_list):
    'Locate the leftmost value smaller or equal to x'
    i = bisect.bisect_left(sorted_list, elem)
    if i != 0:
        return i
    else:
        return 0


      
seq_len = 10**6
sam_size = 100
pop_size = 100000
r = 10**(-8)
mu = 1.3*10**(-8)
totaliter=2



#prior = unif(pop_size+-50000 )  unif(0,0.002)



demo_model = msprime.Demography.isolated_model([pop_size], growth_rate=[0.001])



ts = msprime.sim_ancestry(
    samples=sam_size,
    recombination_rate= r, 
    sequence_length= seq_len,
    #population_size = pop_size,
    #random_seed =seed,
    #discrete_genome=False,
    demography = demo_model
    )
# Visualise the simulated ancestral history.
#SVG(ts.draw_svg())

#ts.num_trees



mts = msprime.sim_mutations(ts, rate = mu,
                            #discrete_genome=False,
                            #random_seed=seed
                          )


print("iter: ",iter," done generating")

sample_data = tsinfer.SampleData.from_tree_sequence(mts, use_sites_time=None, use_individuals_time=None)
inferred_ts = tsinfer.infer(sample_data,path_compression=False)

print("iter:",iter," done infering")

###################################################################################
mts = inferred_ts

msites = mts.tables.mutations.site
sitepos = mts.tables.sites.position
mpos = []
for i in range(len(msites)):
    cursite = msites[i]
    mpos.append(sitepos[cursite])
mpos = numpy.array(mpos)
mnodes = mts.tables.mutations.node
ec = mts.tables.edges.child
ep = mts.tables.edges.parent
el = mts.tables.edges.left
er = mts.tables.edges.right


leftmost = min(el)
rightmost = max(er)
for i in range(len(el)):
    if el[i] == leftmost:
        el[i] = 0
    if er[i] == rightmost:
        er[i] = seq_len


if len(mpos) == 0:
    mpos_sn = []
    mnodes_sn = []
else:
    mnp_sn = numpy.transpose(sorted(numpy.transpose(numpy.array([mnodes,mpos])).tolist()))
    mnodes_sn = mnp_sn[0]
    mpos_sn = mnp_sn[1]

eclrp_sc = numpy.transpose(sorted(numpy.transpose(numpy.array([ec,el,er,ep])).tolist()))
ec_sc = eclrp_sc[0]
el_sc = eclrp_sc[1]
er_sc = eclrp_sc[2]
ep_sc = eclrp_sc[3]
enm_sc = [[] for i in range(len(ec_sc))]

curedge = 0
for i in range(len(mnodes_sn)):
    curnodemu = mnodes_sn[i] 
    cursitemu = mpos_sn[i]
    for j in range(curedge,len(ec_sc)):
        if curnodemu == ec_sc[j]:

            if er_sc[j] > cursitemu and el_sc[j] <= cursitemu:

                enm_sc[j].append(cursitemu)  

                curedge = j
                break
            if el_sc[j] > cursitemu:
                curedge = j
                break

        if ec_sc[j] > mnodes_sn[i] or j == len(ec_sc)-1:
            #print("missing mutations")
            curedge = j
            break


print("iter: ",iter," done sorting")


###########################################################################
###########################################################################
##########################################################################

clocate = []
curc = -1
for i in range(len(ec_sc)):
    while ec_sc[i] != curc:
        curc=curc+1
        clocate.append(i)
clocate.append(seq_len)




abdlist = []

for i in range(0,2*sam_size-1):
    for j in range(i+1,i+2):
        canabd = []
        curleft = 0
        curright =seq_len
        #start = 0
        while curleft != seq_len:
            curp1 = i
            curp2 = j 
            curmed = []
            while curp1 != curp2:
                if curp1 < curp2:
                    begin1 = clocate[curp1]
                    end = clocate[curp1+1]
                    s_range = el_sc[begin1:end]
                    begin2 = find_s_se(curleft,s_range)
                    for k in range(begin1+begin2,end):
                        if el_sc[k]<=curleft and er_sc[k] > curleft:
                            curp1 = int(ep_sc[k])
                            curoldr = curright
                            curright = min(er_sc[k],curright)

                            curmed.append(k)

                            break

                else:
                    begin1 = clocate[curp2]
                    end = clocate[curp2+1]
                    s_range = el_sc[begin1:end]
                    begin2 = find_s_se(curleft,s_range)
                    for k in range(begin1+begin2,end):
                        if el_sc[k]<=curleft and er_sc[k] > curleft:
                            curp2 = int(ep_sc[k])
                            curoldr = curright
                            curright = min(er_sc[k],curright) 

                            curmed.append(k)

                            break

            canmution = []
            for k in range(len(curmed)):
                edgemu = enm_sc[curmed[k]]
                startt = find_s_s(curleft,edgemu)
                endd = find_s_s(curright,edgemu)                    
                canmution = canmution + edgemu[startt:endd] 

            canabd.append([i,j,curleft,curright,curp1,canmution])
            curleft = curright
            curright = seq_len





        if len(canabd)==1:
            abdlist = abdlist + canabd
        else:
            curabd = [canabd[0][0],canabd[0][1],canabd[0][2],canabd[0][3],canabd[0][4],canabd[0][5]]
            for k in range(1,len(canabd)):
                if canabd[k][4] != canabd[k-1][4]:
                    abdlist.append(curabd)
                    curabd = [canabd[k][0],canabd[k][1],canabd[k][2],canabd[k][3],canabd[k][4],canabd[k][5]]
                else:
                    curabd[3] = canabd[k][3]
                    curabd[5] = curabd[5] + canabd[k][5]
                if k == len(canabd)-1:
                    abdlist.append(curabd)

#             if i % math.ceil(2*sam_size/100) ==0:
#                 print("finding abd pair propotion ",i/2/sam_size)

lenlist = numpy.zeros(len(abdlist))
for i in range(len(abdlist)):
    lenlist[i] = abdlist[i][3] - abdlist[i][2]
obsmean = statistics.mean(lenlist)
obsdev = statistics.stdev(lenlist)


####################################################################################################
###################################################################################################


acceptlist = [] 

for iter in range(totaliter):
    
    inipop = numpy.random.uniform(50000,150000)
    increrate = numpy.random.uniform(0,0.002)
    
    demo_model = msprime.Demography.isolated_model([inipop], growth_rate=[increrate])



    ts = msprime.sim_ancestry(
        samples=sam_size,
        recombination_rate= r, 
        sequence_length= seq_len,
        #population_size = pop_size,
        #random_seed =seed,
        #discrete_genome=False,
        demography = demo_model
        )
    # Visualise the simulated ancestral history.
    #SVG(ts.draw_svg())
    
    #ts.num_trees
     
    
    
    mts = msprime.sim_mutations(ts, rate = mu,
                                #discrete_genome=False,
                                #random_seed=seed
                              )
    
 
    
    sample_data = tsinfer.SampleData.from_tree_sequence(mts, use_sites_time=None, use_individuals_time=None)
    inferred_ts = tsinfer.infer(sample_data,path_compression=False)
    
    print("iter:",iter," done infering")

###################################################################################
    mts = inferred_ts

    msites = mts.tables.mutations.site
    sitepos = mts.tables.sites.position
    mpos = []
    for i in range(len(msites)):
        cursite = msites[i]
        mpos.append(sitepos[cursite])
    mpos = numpy.array(mpos)
    mnodes = mts.tables.mutations.node
    ec = mts.tables.edges.child
    ep = mts.tables.edges.parent
    el = mts.tables.edges.left
    er = mts.tables.edges.right


    leftmost = min(el)
    rightmost = max(er)
    for i in range(len(el)):
        if el[i] == leftmost:
            el[i] = 0
        if er[i] == rightmost:
            er[i] = seq_len


    if len(mpos) == 0:
        mpos_sn = []
        mnodes_sn = []
    else:
        mnp_sn = numpy.transpose(sorted(numpy.transpose(numpy.array([mnodes,mpos])).tolist()))
        mnodes_sn = mnp_sn[0]
        mpos_sn = mnp_sn[1]

    eclrp_sc = numpy.transpose(sorted(numpy.transpose(numpy.array([ec,el,er,ep])).tolist()))
    ec_sc = eclrp_sc[0]
    el_sc = eclrp_sc[1]
    er_sc = eclrp_sc[2]
    ep_sc = eclrp_sc[3]
    enm_sc = [[] for i in range(len(ec_sc))]

    curedge = 0
    for i in range(len(mnodes_sn)):
        curnodemu = mnodes_sn[i] 
        cursitemu = mpos_sn[i]
        for j in range(curedge,len(ec_sc)):
            if curnodemu == ec_sc[j]:

                if er_sc[j] > cursitemu and el_sc[j] <= cursitemu:

                    enm_sc[j].append(cursitemu)  

                    curedge = j
                    break
                if el_sc[j] > cursitemu:
                    curedge = j
                    break

            if ec_sc[j] > mnodes_sn[i] or j == len(ec_sc)-1:
                #print("missing mutations")
                curedge = j
                break

 


    ###########################################################################
    ###########################################################################
    ##########################################################################

    clocate = []
    curc = -1
    for i in range(len(ec_sc)):
        while ec_sc[i] != curc:
            curc=curc+1
            clocate.append(i)
    clocate.append(seq_len)




    abdlist = []

    for i in range(0,2*sam_size-1):
        for j in range(i+1,i+2):
            canabd = []
            curleft = 0
            curright =seq_len
            #start = 0
            while curleft != seq_len:
                curp1 = i
                curp2 = j 
                curmed = []
                while curp1 != curp2:
                    if curp1 < curp2:
                        begin1 = clocate[curp1]
                        end = clocate[curp1+1]
                        s_range = el_sc[begin1:end]
                        begin2 = find_s_se(curleft,s_range)
                        for k in range(begin1+begin2,end):
                            if el_sc[k]<=curleft and er_sc[k] > curleft:
                                curp1 = int(ep_sc[k])
                                curoldr = curright
                                curright = min(er_sc[k],curright)

                                curmed.append(k)

                                break

                    else:
                        begin1 = clocate[curp2]
                        end = clocate[curp2+1]
                        s_range = el_sc[begin1:end]
                        begin2 = find_s_se(curleft,s_range)
                        for k in range(begin1+begin2,end):
                            if el_sc[k]<=curleft and er_sc[k] > curleft:
                                curp2 = int(ep_sc[k])
                                curoldr = curright
                                curright = min(er_sc[k],curright) 

                                curmed.append(k)

                                break

                canmution = []
                for k in range(len(curmed)):
                    edgemu = enm_sc[curmed[k]]
                    startt = find_s_s(curleft,edgemu)
                    endd = find_s_s(curright,edgemu)                    
                    canmution = canmution + edgemu[startt:endd] 

                canabd.append([i,j,curleft,curright,curp1,canmution])
                curleft = curright
                curright = seq_len





            if len(canabd)==1:
                abdlist = abdlist + canabd
            else:
                curabd = [canabd[0][0],canabd[0][1],canabd[0][2],canabd[0][3],canabd[0][4],canabd[0][5]]
                for k in range(1,len(canabd)):
                    if canabd[k][4] != canabd[k-1][4]:
                        abdlist.append(curabd)
                        curabd = [canabd[k][0],canabd[k][1],canabd[k][2],canabd[k][3],canabd[k][4],canabd[k][5]]
                    else:
                        curabd[3] = canabd[k][3]
                        curabd[5] = curabd[5] + canabd[k][5]
                    if k == len(canabd)-1:
                        abdlist.append(curabd)

#             if i % math.ceil(2*sam_size/100) ==0:
#                 print("finding abd pair propotion ",i/2/sam_size)

                
                
    ####################################################################################################
    ###################################################################################################
 
    lenlist = numpy.zeros(len(abdlist))
    for i in range(len(abdlist)):
        lenlist[i] = abdlist[i][3] - abdlist[i][2]
    simmean = statistics.mean(lenlist)
    simdev = statistics.stdev(lenlist)
    
    if abs(simmean - obsmean) <=1000 and abs(simdev- obsdev) <=2000:
        acceptlist.append([inipop,increrate])
    

run = 1   
name = "acceptlist"+str(run)+"pop100000_exp001.dat"    
    
file = open(name,"wb")
pickle.dump(acceptlist,file)
file.close()    
    
