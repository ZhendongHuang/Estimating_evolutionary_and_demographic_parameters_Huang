import argparse
import os
import re
import numpy as np

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

    
file = open("recombination_map.dat","rb")
rmap = pickle.load(file)
file.close

sam=tsinfer.load("tgp_chr20.samples")

MSLlist = []
for i in range(len(sam.individuals_population)):
    if sam.individuals_population[i]==13:
        MSLlist.append(i)
newdata = sam.subset(individuals=MSLlist)

        
seq_len = 63025522
sam_size = 85
r = rmap
mu = 1.375*10**(-8)
# gene_conver = 10**(-8)
# track_len = 100
# epsilon = 0.0001
totaliter=1

simtime = 100
# priormu = [1.0*10**(-8), 2.0*10**(-8)]
# prioreps = [0, 0.0002]
# accept_threshold = [1000,400]
freqconsider = [1,40]
priorinipop = [60000, 5000]





inipopsimlist = []

summary_sta1_simlist = []
summary_sta2_simlist = []
summary_sta3_simlist = []
summary_sta4_simlist = []


for iter in range(totaliter):

    sitefreq = numpy.sum(newdata.sites_genotypes,1)
    siteuse =  numpy.where( (sitefreq>freqconsider[0] ) * ( sitefreq<freqconsider[1] ) )[0]  
    print("siteuse_rea",len(siteuse))

    datause = newdata.subset(sites = siteuse)

    inferred_ts = tsinfer.infer(datause,path_compression=False)



    mts = inferred_ts

    ###############################################################################
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
    clocate.append(len(ec_sc))




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



#             print("abdlist sample:",i)
    print("numibd",len(abdlist))
    
    abduse = abdlist

#########################################################################

    tolnummu = 0
    for i in range(len(abduse)):
        tolnummu = tolnummu+len(abduse[i][5])


    muloclist = np.zeros(tolnummu)
    indd = 0
    for i in range(len(abduse)):
        abdc = abduse[i]
        lee = len(abdc[5])
        muloclist[indd:(indd+lee)] = abdc[5]
        indd = indd+lee

    smulist = numpy.sort(muloclist)
    count2 = 0
    curloc = -1
    curcount = 1
    for i in range(len(smulist)):
        if smulist[i] == curloc:
            curcount = curcount+1
        else:
            if curcount ==2:
                count2 = count2+1
            curcount = 1
            curloc = smulist[i]

#         summary_sta1_sim = (len(smulist) - 2*count2)/len(abduse)        
#         summary_sta2_sim = (2*count2) / len(abduse)

    summary_sta1_inf = (len(smulist) - 2*count2)        
    summary_sta2_inf = (2*count2)  




########################################################################            



    lenlist = numpy.zeros(len(abduse))
    for i in range(len(abduse)):
        lenlist[i] = abduse[i][3] - abduse[i][2]

    summary_sta3_inf = statistics.mean(lenlist) 
    summary_sta4_inf = statistics.stdev(lenlist) 

#################################################################################################
    stime = time.time()

    for run in range(simtime):
        
        popsim = numpy.random.uniform(priorinipop[0],priorinipop[1])

#         if run == 0:
#             popsim = 109000

#         if run == 1:
#             popsim = 109000

#         if run == 2:
#             popsim = 4000

#         if run == 3:
#             popsim = 4000
           
#         if run == 4:
#             popsim = 50000
         
#         if run == 5:
#             popsim = 50000
   
        
        
        
        inipopsimlist.append(popsim)

        
        
#         demo_model = msprime.Demography.isolated_model([popsim], growth_rate=[ratesim])
#         demo_model.add_population_parameters_change(800,  initial_size=None, growth_rate=0, population=None)


        genpoint = [0*5/5.2,400*5/5.2,800*5/5.2,1200*5/5.2,1600*5/5.2,2000*5/5.2,2400*5/5.2,2800*5/5.2,\
                    3200*5/5.2,3600*5/5.2,4000*5/5.2,8000*5/5.2,12000*5/5.2,16000*5/5.2,20000*5/5.2]
        inipopsize =[12500,12500,11900,11500,10600,9600,9000,9600,10600,11900,13500,28500,26900,20800,17300]
        genuse = []
        inipuse = []
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
            inipuse = inipuse + addplist
        genuse.append(genpoint[len(genpoint)-1])
        inipuse.append(inipopsize[len(inipopsize)-1])

        
        ind1000 = 0
        for i in range(len(genuse)):
            if genuse[i]<1000:
                ind1000 = i
        ratesim =   -numpy.log(inipuse[ind1000]/popsim )/1000   
        for i in range(len(genuse)):
            if genuse[i]<1000:
                inipuse[i] = popsim*numpy.exp(-ratesim*genuse[i]) 
        for i in range(len(genuse)-1):
            gruse.append( numpy.log( inipuse[i]/inipuse[i+1] )/(genuse[i+1] - genuse[i])  )
        gruse.append(0)
        

        demo_model = msprime.Demography.isolated_model([inipuse[0]], growth_rate=[gruse[0]])
        for i in range(1,len(genuse)):
            demo_model.add_population_parameters_change(genuse[i],  initial_size=inipuse[i], growth_rate=gruse[i], population=None)

                
                

#         for i in range(len(genuse)-1):
#             gruse.append( numpy.log( inipuse[i]/inipuse[i+1] )/(genuse[i+1] - genuse[i])  )
#         gruse.append(0)
        

#         demo_model = msprime.Demography.isolated_model([inipuse[0]], growth_rate=[gruse[0]])
#         for i in range(1,len(genuse)):
#             demo_model.add_population_parameters_change(genuse[i],  initial_size=inipuse[i], growth_rate=gruse[i], population=None)





        ts = msprime.sim_ancestry(
            samples=sam_size,
            recombination_rate= r, 
            sequence_length= seq_len,
            record_provenance=False,
#             population_size = pop_size,
#             gene_conversion_rate = gene_conver ,
#             gene_conversion_tract_length = track_len,
            #random_seed =seed,
            #discrete_genome=False
            demography = demo_model
            )
        # Visualise the simulated ancestral history.
        #SVG(ts.draw_svg())

        #ts.num_trees
        
        print("done gene")

        mts = msprime.sim_mutations(ts, rate = mu,
                                    #discrete_genome=False,
                                    #random_seed=seed
                                  )

        print("done mu")

        sample_data = tsinfer.SampleData.from_tree_sequence(mts, use_sites_time=None, use_individuals_time=None)

    
        sitefreq = numpy.sum(sample_data.sites_genotypes,1)
        siteuse =  numpy.where( (sitefreq>freqconsider[0] ) * ( sitefreq<freqconsider[1] ) )[0]  
        print("siteuse_sim",len(siteuse))

        datause = sample_data.subset(sites = siteuse)

        inferred_ts = tsinfer.infer(datause,path_compression=False)



        mts = inferred_ts

        ###############################################################################
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
        clocate.append(len(ec_sc))




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


    #             print("abdlist sample:",i)
        print("numibd",len(abdlist))

        abduse = abdlist


        tolnummu = 0
        for i in range(len(abduse)):
            tolnummu = tolnummu+len(abduse[i][5])


        muloclist = np.zeros(tolnummu)
        indd = 0
        for i in range(len(abduse)):
            abdc = abduse[i]
            lee = len(abdc[5])
            muloclist[indd:(indd+lee)] = abdc[5]
            indd = indd+lee
            
        smulist = numpy.sort(muloclist)
        count2 = 0
        curloc = -1
        curcount = 1
        for i in range(len(smulist)):
            if smulist[i] == curloc:
                curcount = curcount+1
            else:
                if curcount ==2:
                    count2 = count2+1
                curcount = 1
                curloc = smulist[i]
        
#         summary_sta1_sim = (len(smulist) - 2*count2)/len(abduse)        
#         summary_sta2_sim = (2*count2) / len(abduse)

        summary_sta1_sim = (len(smulist) - 2*count2)        
        summary_sta2_sim = (2*count2)  

        summary_sta1_simlist.append(summary_sta1_sim)
        summary_sta2_simlist.append(summary_sta2_sim)
        
        
########################################################################            
        
        

        lenlist = numpy.zeros(len(abduse))
        for i in range(len(abduse)):
            lenlist[i] = abduse[i][3] - abduse[i][2]

        summary_sta3_simlist.append( statistics.mean(lenlist) )
        summary_sta4_simlist.append( statistics.stdev(lenlist) )
        
        print("-----------------------------------")

#         print(inipopsimlist )
#         print(iniratesimlist )
        print(summary_sta1_simlist,summary_sta1_inf)
        print(summary_sta2_simlist,summary_sta2_inf)
        print(summary_sta3_simlist,summary_sta3_inf)
        print(summary_sta4_simlist,summary_sta4_inf)
        print(time.time()-stime, run)
        print("-----------------------------------")




byrun = 1
name1 = "1M2_ne_realpopsimall_"+"run="+str(byrun) +"MSL_efibd.dat"
name3 = "1M2_ne_realsta1simall_"+"run="+str(byrun) +"MSL_efibd.dat"
name4 = "1M2_ne_realsta2simall_"+"run="+str(byrun) +"MSL_efibd.dat"
name5 = "1M2_ne_realsta1inf_"+"run="+str(byrun) +"MSL_efibd.dat"
name6 = "1M2_ne_realsta2inf_"+"run="+str(byrun) +"MSL_efibd.dat"
name7 = "1M2_ne_realsta3simall_"+"run="+str(byrun) +"MSL_efibd.dat"
name8 = "1M2_ne_realsta4simall_"+"run="+str(byrun) +"MSL_efibd.dat"
name9 = "1M2_ne_realsta3inf_"+"run="+str(byrun) +"MSL_efibd.dat"
name10 = "1M2_ne_realsta4inf_"+"run="+str(byrun) +"MSL_efibd.dat"
file = open(name1,"wb") 
pickle.dump(inipopsimlist,file)
file.close()
file = open(name3,"wb") 
pickle.dump(summary_sta1_simlist,file)
file.close()
file = open(name4,"wb") 
pickle.dump(summary_sta2_simlist,file)
file.close()
file = open(name5,"wb") 
pickle.dump(summary_sta1_inf,file)
file.close()
file = open(name6,"wb") 
pickle.dump(summary_sta2_inf,file)
file.close()
file = open(name7,"wb") 
pickle.dump(summary_sta3_simlist,file)
file.close()
file = open(name8,"wb") 
pickle.dump(summary_sta4_simlist,file)
file.close()
file = open(name9,"wb") 
pickle.dump(summary_sta3_inf,file)
file.close()
file = open(name10,"wb") 
pickle.dump(summary_sta4_inf,file)
file.close()