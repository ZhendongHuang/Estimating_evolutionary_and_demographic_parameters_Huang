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
epsilon = 0.0001
totaliter=1

simtime = 2500
priormu = [1.0*10**(-8), 2.0*10**(-8)]
prioreps = [0.00006, 0.00016]
priorpop = [7000,15000]
priortau = [-2*10**(-5),10**(-5)]



muhatlist = []
epshatlist = []
muacceptlist = []
epsacceptlist = []
pophatlist = []
tauhatlist = []
popacceptlist = []
tauacceptlist = []
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

    ii=0
    for i in range(0,2*sam_size):
        if i==2*sam_size-1:
            ii = 1
        for j in range(i+1,i+2):
            if ii ==1:
                j = 1
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



    
    #compute estimators 
    #cstart store the starting index of the IBD with child i (even number)
    cstart = [-1]*(2*sam_size)
    curcc = -3
    
    # allmute count the total number of mut + type error, in all IBDs.
    allmute = 0
    alll = 0
    allnum = 0
    
    maxnode = -1
    
    lenlist = []
    
    for i in range(len(abduse) ):        
        abdc = abduse[i]    
        
        lenlist.append(abdc[3]-abdc[2])

 
        if abdc[0] > curcc:
            curcc = abdc[0]
            cstart[curcc] = i
           
        
        if abdc[4]>maxnode:
            maxnode = abdc[4]
            
        allmute = allmute + len(abdc[5])
        alll = alll + abdc[3]-abdc[2]
        allnum = allnum+1
    
    maxnode = maxnode+1    
    # minmute compares pair (01,12),(12,23)..., and count only the mute + error
    # in the middle sample, occur before first coalescent, after averaging. If the middle sample
    # does not coalest first, can not be count. but the total will be the counted*3/2. 
   
    minmute = 0
    
    # numminre count the num of recombination happen at middle sammple and before any coalescent. The proportion 
    # of such event and can be detected is 4/9. Thus the total num of recom occur ant middle sample and before 
    # any coalescent is numminre*9/4. 
    
    numminre = 0
    
    # minl compute the length of error happen at the middle sample, the total length will be minl*3/2. 
    
    minl =0 
    
    
    therange = range(2*sam_size-2)
    
    
    for i in therange:
        
        ind1 = cstart[i]
        ind2 = cstart[i+1]
        right = 0
        left = 0
        abd1 = abduse[ind1]
        abd2 = abduse[ind2]
        r1 = abd1[3]
        r2 = abd2[3]
        while right != seq_len:
            if r1<r2:
                ances1 = abd1[4]
                ances2 = abd2[4]               
                minmute = minmute + (ances1< ances2)*len(abd1[5])*(r1-left)/(abd1[3]-abd1[2])/2 \
                               + (ances1 > ances2)*len(abd2[5])*(r1-left)/(abd2[3]-abd2[2])/2   
                
                minl = minl + (r1-left)*(ances1< ances2) + (r1-left)*(ances1 > ances2)
      
                left = r1
                ind1 = ind1+1
                abd1 = abduse[ind1]
                r1 = abd1[3]

            elif r2<r1:
                ances1 = abd1[4]
                ances2 = abd2[4]              
                minmute = minmute + (ances1< ances2)*len(abd1[5])*(r2-left)/(abd1[3]-abd1[2])/2 \
                               + (ances1 > ances2)*len(abd2[5])*(r2-left)/(abd2[3]-abd2[2])/2
                 
                minl = minl + (r2-left)*(ances1< ances2) + (r2-left)*(ances1 > ances2)
                
                left = r2
                ind2 = ind2+1
                abd2 = abduse[ind2]
                r2 = abd2[3]

            else:
                ances1 = abd1[4]
                ances2 = abd2[4]              
                minmute = minmute + (ances1< ances2)*len(abd1[5])*(r2-left)/(abd1[3]-abd1[2])/2 \
                               + (ances1 > ances2)*len(abd2[5])*(r2-left)/(abd2[3]-abd2[2])/2
 
                minl = minl + (r2-left)*(ances1< ances2) + (r2-left)*(ances1 > ances2)
           
                left = r2
                ind2 = ind2+1
                abd2 = abduse[ind2]
                r2 = abd2[3]
                ind1 = ind1+1
                abd1 = abduse[ind1]
                r1 = abd1[3]
                             
                if (ances1 != ances2) or (abd1[4] != abd2[4]):
                    numminre = numminre + 1
                
            
            right = min(r1,r2)
                   
    #C1=3/2*minmute, A1=numminre*9/4, B1=minl*3/2
    #C2=allmute/allnum, A2=3/2, B2=2*alll/allnum
    #Solve: C1=A1(mu/r) + B1 epsilon,  C2= A2(mu/r) + B2 epsilon    
    
    C1 = 3/2*minmute
    A1 = numminre*9/4
    B1 = minl*3/2
    
    C2 = allmute/allnum
    A2 = 3/2
    B2 = 2*alll/allnum

    
    muhat = (B2*C1 -B1*C2)/(A1*B2-A2*B1)*r
    epsilonhat = (A2*C1 - A1*C2)/(A2*B1 - A1*B2)
    

    summary_sta1_inf = allmute  
    summary_sta2_inf = C1
    summary_sta3_inf = statistics.mean(lenlist)
    summary_sta4_inf = statistics.stdev(lenlist)







##########################################################################
    
    muaccept = []
    epsaccept = []
    popaccept = []
    tauaccept = []
    summary_sta1_simlist=[]
    summary_sta2_simlist=[]
    summary_sta3_simlist=[]
    summary_sta4_simlist=[]
    stime = time.time()

    for run in range(simtime):
        
        musim = numpy.random.uniform(priormu[0],priormu[1])
        epssim = numpy.random.uniform(prioreps[0],prioreps[1])
        if epssim<0:
            epssim = 0

        
        popsim = numpy.random.uniform(priorpop[0],priorpop[1])
        tausim = numpy.random.uniform(priortau[0],priortau[1])
        ## put a population cap = 20000*2
        
        if tausim <0:
            turn = numpy.log(20000/popsim)/(-tausim)
            demo_model = msprime.Demography.isolated_model( [popsim],  growth_rate=[tausim] )   
            demo_model.add_population_parameters_change(turn,  initial_size=None, growth_rate=0, population=None)
            
        else:
            demo_model = msprime.Demography.isolated_model([popsim], growth_rate=[tausim])  



        ts = msprime.sim_ancestry(
            samples=sam_size,
            recombination_rate= r, 
            sequence_length= seq_len,
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

        ii=0
        for i in range(0,2*sam_size):
            if i==2*sam_size-1:
                ii = 1
            for j in range(i+1,i+2):
                if ii ==1:
                    j = 1
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





        #compute estimators 
        #cstart store the starting index of the IBD with child i (even number)
        cstart = [-1]*(2*sam_size)
        curcc = -3

        # allmute count the total number of mut + type error, in all IBDs.
        allmute = 0
        alll = 0
        allnum = 0

        maxnode = -1

        lenlist = []

        for i in range(len(abduse) ):        
            abdc = abduse[i]    

            lenlist.append(abdc[3]-abdc[2])


            if abdc[0] > curcc:
                curcc = abdc[0]
                cstart[curcc] = i


            if abdc[4]>maxnode:
                maxnode = abdc[4]

            allmute = allmute + len(abdc[5])
            alll = alll + abdc[3]-abdc[2]
            allnum = allnum+1

        maxnode = maxnode+1    
        # minmute compares pair (01,12),(12,23)..., and count only the mute + error
        # in the middle sample, occur before first coalescent, after averaging. If the middle sample
        # does not coalest first, can not be count. but the total will be the counted*3/2. 

        minmute = 0

        # numminre count the num of recombination happen at middle sammple and before any coalescent. The proportion 
        # of such event and can be detected is 4/9. Thus the total num of recom occur ant middle sample and before 
        # any coalescent is numminre*9/4. 

        numminre = 0

        # minl compute the length of error happen at the middle sample, the total length will be minl*3/2. 

        minl =0 


        therange = range(2*sam_size-2)


        for i in therange:

            ind1 = cstart[i]
            ind2 = cstart[i+1]
            right = 0
            left = 0
            abd1 = abduse[ind1]
            abd2 = abduse[ind2]
            r1 = abd1[3]
            r2 = abd2[3]
            while right != seq_len:
                if r1<r2:
                    ances1 = abd1[4]
                    ances2 = abd2[4]               
                    minmute = minmute + (ances1< ances2)*len(abd1[5])*(r1-left)/(abd1[3]-abd1[2])/2 \
                                   + (ances1 > ances2)*len(abd2[5])*(r1-left)/(abd2[3]-abd2[2])/2   

                    minl = minl + (r1-left)*(ances1< ances2) + (r1-left)*(ances1 > ances2)

                    left = r1
                    ind1 = ind1+1
                    abd1 = abduse[ind1]
                    r1 = abd1[3]

                elif r2<r1:
                    ances1 = abd1[4]
                    ances2 = abd2[4]              
                    minmute = minmute + (ances1< ances2)*len(abd1[5])*(r2-left)/(abd1[3]-abd1[2])/2 \
                                   + (ances1 > ances2)*len(abd2[5])*(r2-left)/(abd2[3]-abd2[2])/2

                    minl = minl + (r2-left)*(ances1< ances2) + (r2-left)*(ances1 > ances2)

                    left = r2
                    ind2 = ind2+1
                    abd2 = abduse[ind2]
                    r2 = abd2[3]

                else:
                    ances1 = abd1[4]
                    ances2 = abd2[4]              
                    minmute = minmute + (ances1< ances2)*len(abd1[5])*(r2-left)/(abd1[3]-abd1[2])/2 \
                                   + (ances1 > ances2)*len(abd2[5])*(r2-left)/(abd2[3]-abd2[2])/2

                    minl = minl + (r2-left)*(ances1< ances2) + (r2-left)*(ances1 > ances2)

                    left = r2
                    ind2 = ind2+1
                    abd2 = abduse[ind2]
                    r2 = abd2[3]
                    ind1 = ind1+1
                    abd1 = abduse[ind1]
                    r1 = abd1[3]

                    if (ances1 != ances2) or (abd1[4] != abd2[4]):
                        numminre = numminre + 1


                right = min(r1,r2)

        #C1=3/2*minmute, A1=numminre*9/4, B1=minl*3/2
        #C2=allmute/allnum, A2=3/2, B2=2*alll/allnum
        #Solve: C1=A1(mu/r) + B1 epsilon,  C2= A2(mu/r) + B2 epsilon    

        C1 = 3/2*minmute
        A1 = numminre*9/4
        B1 = minl*3/2

        C2 = allmute/allnum
        A2 = 3/2
        B2 = 2*alll/allnum


        muhat = (B2*C1 -B1*C2)/(A1*B2-A2*B1)*r
        epsilonhat = (A2*C1 - A1*C2)/(A2*B1 - A1*B2)

    
        summary_sta1_sim = allmute       
        summary_sta2_sim = C1
        summary_sta3_sim = statistics.mean(lenlist)
        summary_sta4_sim = statistics.stdev(lenlist)

        summary_sta1_simlist.append(summary_sta1_sim)
        summary_sta2_simlist.append(summary_sta2_sim)
        summary_sta3_simlist.append(summary_sta3_sim)
        summary_sta4_simlist.append(summary_sta4_sim)

#         
#         if abs(summary_sta1_sim - summary_sta1_inf) < accept_threshold[0] and abs(summary_sta2_sim - summary_sta2_inf) < accept_threshold[1]:
        muaccept.append(musim)
        epsaccept.append(epssim)
        popaccept.append(popsim)
        tauaccept.append(tausim)

            
        etime = time.time()

            
        print("run",run,etime-stime)
        
        
    muacceptlist = muacceptlist + muaccept
    epsacceptlist= epsacceptlist + epsaccept
    popacceptlist = popacceptlist + popaccept
    tauacceptlist = tauacceptlist + tauaccept
    
    
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
name1 = "fmusimall_seq_len="+str(seq_len)+"run="+str(byrun) +"_sam_size="+str(sam_size)+"_eps(est)=" + str(epsilon)+"_noconver_constne(est)_efibd.dat"
name2 = "fepssimall_seq_len="+str(seq_len)+"run="+str(byrun) +"_sam_size="+str(sam_size)+"_eps(est)=" + str(epsilon)+"_noconver_constne(est)_efibd.dat"
name3 = "fsta1simall_seq_len="+str(seq_len)+"run="+str(byrun) +"_sam_size="+str(sam_size)+"_eps(est)=" + str(epsilon)+"_noconver_constne(est)_efibd.dat"
name4 = "fsta2simall_seq_len="+str(seq_len)+"run="+str(byrun) +"_sam_size="+str(sam_size)+"_eps(est)=" + str(epsilon)+"_noconver_constne(est)_efibd.dat"
name5 = "fsta1inf_seq_len="+str(seq_len)+"run="+str(byrun) +"_sam_size="+str(sam_size)+"_eps(est)=" + str(epsilon)+"_noconver_constne(est)_efibd.dat"
name6 = "fsta2inf_seq_len="+str(seq_len)+"run="+str(byrun) +"_sam_size="+str(sam_size)+"_eps(est)=" + str(epsilon)+"_noconver_constne(est)_efibd.dat"
name7 = "fabdlen_seq_len="+str(seq_len)+"run="+str(byrun) +"_sam_size="+str(sam_size)+"_eps(est)=" + str(epsilon)+"_noconver_constne(est)_efibd.dat"
name8 = "fpopsimall_seq_len="+str(seq_len)+"run="+str(byrun) +"_sam_size="+str(sam_size)+"_eps(est)=" + str(epsilon)+"_noconver_constne(est)_efibd.dat"
name9 = "ftausimall_seq_len="+str(seq_len)+"run="+str(byrun) +"_sam_size="+str(sam_size)+"_eps(est)=" + str(epsilon)+"_noconver_constne(est)_efibd.dat"
name10 = "fsta3simall_seq_len="+str(seq_len)+"run="+str(byrun) +"_sam_size="+str(sam_size)+"_eps(est)=" + str(epsilon)+"_noconver_constne(est)_efibd.dat"
name11 = "fsta4simall_seq_len="+str(seq_len)+"run="+str(byrun) +"_sam_size="+str(sam_size)+"_eps(est)=" + str(epsilon)+"_noconver_constne(est)_efibd.dat"
name12 = "fsta3inf_seq_len="+str(seq_len)+"run="+str(byrun) +"_sam_size="+str(sam_size)+"_eps(est)=" + str(epsilon)+"_noconver_constne(est)_efibd.dat"
name13 = "fsta4inf_seq_len="+str(seq_len)+"run="+str(byrun) +"_sam_size="+str(sam_size)+"_eps(est)=" + str(epsilon)+"_noconver_constne(est)_efibd.dat"





file = open(name1,"wb") 
pickle.dump(muacceptlist,file)
file.close()

file = open(name2,"wb") 
pickle.dump(epsacceptlist,file)
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
pickle.dump(lllabd,file)
file.close()

file = open(name8,"wb") 
pickle.dump(popacceptlist,file)
file.close()

file = open(name9,"wb") 
pickle.dump(tauacceptlist,file)
file.close()

file = open(name10,"wb") 
pickle.dump(summary_sta3_simlist,file)
file.close()

file = open(name11,"wb") 
pickle.dump(summary_sta4_simlist,file)
file.close()

file = open(name12,"wb") 
pickle.dump(summary_sta3_inf,file)
file.close()

file = open(name13,"wb") 
pickle.dump(summary_sta4_inf,file)
file.close()