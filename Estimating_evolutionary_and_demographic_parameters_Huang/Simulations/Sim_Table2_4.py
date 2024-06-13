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

        
seq_len = 10**7
sam_size = 5
pop_size = 10000
r = 10**(-8)
mu = 1.3*10**(-8)
# gene_conver = 10**(-8)
# track_len = 100
epsilon = 0
totaliter=1

simtime = 1
priormu = [1.0*10**(-8), 2.0*10**(-8)]
prioreps = [0, 0]
accept_threshold = [1000,400]






muhatlist = []
epshatlist = []
muacceptlist = []
epsacceptlist = []
for iter in range(totaliter):
    
    ts = msprime.sim_ancestry(
        samples=sam_size,
        recombination_rate= r, 
        sequence_length= seq_len,
        population_size = pop_size,
#         gene_conversion_rate = gene_conver ,
#         gene_conversion_tract_length = track_len,
        #random_seed =seed,
        #discrete_genome=False
        )
    # Visualise the simulated ancestral history.
    #SVG(ts.draw_svg())

    #ts.num_trees

    mts = msprime.sim_mutations(ts, rate = mu,
                                #discrete_genome=False,
                                #random_seed=seed
                              )

    sample_data = tsinfer.SampleData.from_tree_sequence(mts, use_sites_time=None, use_individuals_time=None)
    
    print("done1")
        
##########################################################################################

    newdata = tsinfer.SampleData()

    gtypetable1 = sample_data.sites_genotypes 
    gtypetable = gtypetable1[0:len(gtypetable1)]
    gtypetable = gtypetable.astype("int")

    alletable1 = sample_data.sites_alleles 
    alletable = alletable1[0:len(alletable1)]

    posittable1 = sample_data.sites_position 
    posittable = posittable1[0:len(posittable1)]

    for i in range(len(sample_data.sites_position) - 1):
        alle = alletable[i]
        posit = posittable[i]
        posit1 = posittable[i+1]
        gtype = gtypetable[i].tolist()
        newdata.add_site(position = posit, genotypes = gtype,alleles = alle)

        numerror = numpy.random.poisson( (posit1-posit)*2*sam_size*epsilon )
        locerror = numpy.sort(numpy.floor(numpy.random.uniform(posit, posit1-0.000000001,numerror) )  ) .tolist()

        curloc = posit
        for j in range(len(locerror)):
            if locerror[j] > curloc:
                curloc = locerror[j]
                gtype = numpy.zeros(2*sam_size  )
                gtype[int( numpy.random.uniform(0,2*sam_size-0.00001)   ) ] = 1
                gtype = gtype.astype('int').tolist()


                newdata.add_site(position = curloc, genotypes = gtype,alleles = ['A','T'])

    alle = alletable[-1]
    posit = posittable[-1]
    gtype = gtypetable[-1].tolist()
    newdata.add_site(position = posit, genotypes = gtype,alleles = alle)


    newdata.finalise()



################################################################################


    inferred_ts = tsinfer.infer(newdata,path_compression=False)

#############################################################################

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
        for j in range(i+1,2*sam_size):
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


    lllabd = len(abdlist)
    print(lllabd)
###########################################################################
###########################################################################
##########################################################################

#     #find time
#     #nt = mts.tables.nodes.time
#     nt = ts.tables.nodes.time
#     for i in range(len(abdlist)):
#         parent = int(abdlist[i][4])
#         abdlist[i].append(-1)

#         # add typing error
#         abdlist[i].append( len(abdlist[i][5])  )



###########################################################################
###########################################################################
##########################################################################

    
    abduse = abdlist

    allmute = 0
    numabdd = 0
    for i in range(len(abduse) ):        
        abdc = abduse[i]        

 
        if abdc[3] != seq_len:
            allmute = allmute + len(abdc[5])
            numabdd = numabdd+1   
 

    
    summary_sta1_inf = allmute    






##########################################################################
    
    muaccept = []
    epsaccept = []
    summary_sta1_simlist=[]
    summary_sta2_simlist=[]
    stime = time.time()

    for run in range(simtime):
        
        musim = numpy.random.uniform(priormu[0],priormu[1])
        epssim = numpy.random.uniform(prioreps[0],prioreps[1])

 
        ts = msprime.sim_ancestry(
            samples=sam_size,
            recombination_rate= r, 
            sequence_length= seq_len,
            population_size = pop_size,
#             gene_conversion_rate = gene_conver ,
#             gene_conversion_tract_length = track_len,
            #random_seed =seed,
            #discrete_genome=False
            )
        # Visualise the simulated ancestral history.
        #SVG(ts.draw_svg())

        #ts.num_trees


        mts = msprime.sim_mutations(ts, rate = musim,
                                    #discrete_genome=False,
                                    #random_seed=seed
                                  )



        sample_data = tsinfer.SampleData.from_tree_sequence(mts, use_sites_time=None, use_individuals_time=None)



##########################################################################################

        newdata = tsinfer.SampleData()

        gtypetable1 = sample_data.sites_genotypes 
        gtypetable = gtypetable1[0:len(gtypetable1)]
        gtypetable = gtypetable.astype("int")

        alletable1 = sample_data.sites_alleles 
        alletable = alletable1[0:len(alletable1)]

        posittable1 = sample_data.sites_position 
        posittable = posittable1[0:len(posittable1)]

        for i in range(len(sample_data.sites_position) - 1):
            alle = alletable[i]
            posit = posittable[i]
            posit1 = posittable[i+1]
            gtype = gtypetable[i].tolist()
            newdata.add_site(position = posit, genotypes = gtype,alleles = alle)

            numerror = numpy.random.poisson( (posit1-posit)*2*sam_size*epssim )
            locerror = numpy.sort(numpy.floor(numpy.random.uniform(posit, posit1-0.000000001,numerror) )  ) .tolist()

            curloc = posit
            for j in range(len(locerror)):
                if locerror[j] > curloc:
                    curloc = locerror[j]
                    gtype = numpy.zeros(2*sam_size  )
                    gtype[int( numpy.random.uniform(0,2*sam_size-0.00001)   ) ] = 1
                    gtype = gtype.astype('int').tolist()


                    newdata.add_site(position = curloc, genotypes = gtype,alleles = ['A','T'])

        alle = alletable[-1]
        posit = posittable[-1]
        gtype = gtypetable[-1].tolist()
        newdata.add_site(position = posit, genotypes = gtype,alleles = alle)


        newdata.finalise()

        
################################################################################


        inferred_ts = tsinfer.infer(newdata,path_compression=False)


        
        
#############################################################################

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
            for j in range(i+1,2*sam_size):
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



    ###########################################################################
    ###########################################################################
    ##########################################################################

#         #find time
#         #nt = mts.tables.nodes.time
#         nt = ts.tables.nodes.time
#         for i in range(len(abdlist)):
#             parent = int(abdlist[i][4])
#             abdlist[i].append(-1)

#             # total typing error
#             abdlist[i].append( len(abdlist[i][5])  )
            



    ###########################################################################
    ###########################################################################
    ##########################################################################


        abduse = abdlist
        
        
        allmute = 0
        numabdd = 0
        for i in range(len(abduse) ):        
            abdc = abduse[i]        


            if abdc[3] != seq_len:
                allmute = allmute + len(abdc[5])
                numabdd = numabdd+1   

    
        summary_sta1_sim = allmute       


        summary_sta1_simlist.append(summary_sta1_sim)
#         
#         if abs(summary_sta1_sim - summary_sta1_inf) < accept_threshold[0] and abs(summary_sta2_sim - summary_sta2_inf) < accept_threshold[1]:
        muaccept.append(musim)

            
        etime = time.time()

            
        print("run",run,etime-stime)
        
        
    muacceptlist = muacceptlist + muaccept
    
    
#     if len(muaccept)>0:    
#         muhatlist.append(statistics.mean(muaccept))
#     else:
#         muhatlist.append(None)    
        
#     if len(epsaccept)>0:    
#         epshatlist.append(statistics.mean(epsaccept))
#     else:
#         epshatlist.append(None)  
    
    print('iter',iter)
#     print(muhatlist)

byrun = 1
name1 = "zmusimall_seq_len="+str(seq_len)+"trun0"+"run="+str(byrun) +"_sam_size="+str(sam_size)+"_eps(no)=" + str(epsilon)+"_noconver_constne_tiibd.dat"
name3 = "zsta1simall_seq_len="+str(seq_len)+"trun0"+"run="+str(byrun) +"_sam_size="+str(sam_size)+"_eps(no)=" + str(epsilon)+"_noconver_constne_tiibd.dat"
name5 = "zsta1inf_seq_len="+str(seq_len)+"trun0"+"run="+str(byrun) +"_sam_size="+str(sam_size)+"_eps(no)=" + str(epsilon)+"_noconver_constne_tiibd.dat"
name7 = "zabdlen_seq_len="+str(seq_len)+"trun0"+"run="+str(byrun) +"_sam_size="+str(sam_size)+"_eps(no)=" + str(epsilon)+"_noconver_constne_tiibd.dat"

# file = open(name1,"wb") 
# pickle.dump(muacceptlist,file)
# file.close()

# file = open(name3,"wb") 
# pickle.dump(summary_sta1_simlist,file)
# file.close()

# file = open(name5,"wb") 
# pickle.dump(summary_sta1_inf,file)
# file.close()

# file = open(name7,"wb") 
# pickle.dump(lllabd,file)
# file.close()