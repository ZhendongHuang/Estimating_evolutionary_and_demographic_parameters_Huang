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
    


seq_len = 10**8
sam_size = 5
pop_size = 10000
r = 10**(-8)
mu = 1.3*10**(-8)
epsilon = 0.001
totaliter=30

name1 = "muhat_true_" + "eps=" + str(epsilon) + "_sam_size=" + str(sam_size) + "constne_len108.dat"
name2 = "epsilonhat_true_" + "eps=" + str(epsilon) + "_sam_size=" + str(sam_size) + "constne_len108.dat"


muhatlist = []
epsilonhatlist = []


for iter in range(totaliter):
    ts = msprime.sim_ancestry(
        samples=sam_size,
        recombination_rate= r, 
        sequence_length= seq_len,
        population_size = pop_size,
        #random_seed =seed,
        #discrete_genome=False,
        #demography = demo_model
        )
    # Visualise the simulated ancestral history.
    #SVG(ts.draw_svg())
    
    #ts.num_trees
     
    
    mts = msprime.sim_mutations(ts, rate = mu,
                                #discrete_genome=False,
                                #random_seed=seed
                              )

    print("iter: ",iter," done generating")
    
    ###########################################################################
    ###########################################################################
    ##########################################################################

 

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
    enm_sc = numpy.zeros(len(ec_sc)) 

    curedge = 0
    for i in range(len(mnodes_sn)):
        curnodemu = mnodes_sn[i]
        for j in range(curedge,len(ec_sc)):
            if curnodemu == ec_sc[j]:
                if er_sc[j] > mpos_sn[i] and el_sc[j] <= mpos_sn[i]:
                    enm_sc[j] = enm_sc[j] + 1
                    curedge = j
                    break
                if el_sc[j] > mpos_sn[i]:
                    curedge = j
                    break
    
            if ec_sc[j] > mnodes_sn[i] or j == len(ec_sc)-1:
                #print("missing mutations")
                curedge = j
                break
            
            
            
    #epclrnm_sp = numpy.transpose(sorted(numpy.transpose(numpy.array([ep_sc,ec_sc,el_sc,er_sc,enm_sc])).tolist()))
    #ep_sp = epclrnm_sp[0]
    #ec_sp = epclrnm_sp[1]
    #el_sp = epclrnm_sp[2]
    #er_sp = epclrnm_sp[3]
    #enm_sp = epclrnm_sp[4]

    
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
                
                
    ################################################                
 
    ################################################                

    abdlist = []
 
    for i in range(0,2*sam_size-1):
        for j in range(i+1,i+2):
            canabd = []
            curleft = 0
            curright =seq_len
            #start = 0
            while curleft != seq_len:
                curnm = 0
                curp1 = i
                curp2 = j 
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
                                curnm=curnm*(curright-curleft)/(curoldr-curleft)+enm_sc[k]*(curright-curleft)/(er_sc[k]-el_sc[k])
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
                                curnm=curnm*(curright-curleft)/(curoldr-curleft)+enm_sc[k]*(curright-curleft)/(er_sc[k]-el_sc[k])    
                                break
                                
                   
          

                canabd.append([i,j,curleft,curright,curp1,curnm])
                curleft = curright
                curright = seq_len
                
          
                
                

            if len(canabd)==1:
                abdlist = abdlist + canabd
            else:
                curabd = [canabd[0][0],canabd[0][1],canabd[0][2],canabd[0][3],canabd[0][4],canabd[0][5] ]
                for k in range(1,len(canabd)):
                    if canabd[k][4] != canabd[k-1][4]:
                        abdlist.append(curabd)
                        curabd = [canabd[k][0],canabd[k][1],canabd[k][2],canabd[k][3],canabd[k][4],canabd[k][5] ]
                    else:
                        curabd[3] = canabd[k][3]
                        curabd[5] = curabd[5] + canabd[k][5]
                    if k == len(canabd)-1:
                        abdlist.append(curabd)
            

    

    print("iter: ",iter," done finding abdlist")


###########################################################################
###########################################################################
##########################################################################




    #find time
    #nt = mts.tables.nodes.time
    nt = ts.tables.nodes.time
    for i in range(len(abdlist)):
        parent = int(abdlist[i][4])
        abdlist[i].append(nt[parent])

        # generate typing error
        muplustype = abdlist[i][5] + numpy.random.poisson((abdlist[i][3]-abdlist[i][2])*2*epsilon )
        abdlist[i].append(muplustype)





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
    
    for i in range(len(abduse) ):        
        abdc = abduse[i]        

 
        if abdc[0] > curcc:
            curcc = abdc[0]
            cstart[curcc] = i
           
        
        if abdc[4]>maxnode:
            maxnode = abdc[4]
            
        allmute = allmute + abdc[7]
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
                minmute = minmute + (ances1< ances2)*abd1[7]*(r1-left)/(abd1[3]-abd1[2])/2 \
                               + (ances1 > ances2)*abd2[7]*(r1-left)/(abd2[3]-abd2[2])/2   
                
                minl = minl + (r1-left)*(ances1< ances2) + (r1-left)*(ances1 > ances2)
      
                left = r1
                ind1 = ind1+1
                abd1 = abduse[ind1]
                r1 = abd1[3]

            elif r2<r1:
                ances1 = abd1[4]
                ances2 = abd2[4]              
                minmute = minmute + (ances1< ances2)*abd1[7]*(r2-left)/(abd1[3]-abd1[2])/2 \
                               + (ances1 > ances2)*abd2[7]*(r2-left)/(abd2[3]-abd2[2])/2
                 
                minl = minl + (r2-left)*(ances1< ances2) + (r2-left)*(ances1 > ances2)
                
                left = r2
                ind2 = ind2+1
                abd2 = abduse[ind2]
                r2 = abd2[3]

            else:
                ances1 = abd1[4]
                ances2 = abd2[4]              
                minmute = minmute + (ances1< ances2)*abd1[7]*(r2-left)/(abd1[3]-abd1[2])/2 \
                               + (ances1 > ances2)*abd2[7]*(r2-left)/(abd2[3]-abd2[2])/2
 
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
    
    print(iter)
    
    muhatlist.append(muhat)
    epsilonhatlist.append(epsilonhat)
    
    
    ####################################################################################################
    ###################################################################################################
    ################################################################################################


file = open(name1,"wb")
pickle.dump(muhatlist,file)
file.close()

file = open(name2,"wb")
pickle.dump(epsilonhat,file)
file.close()
    
    
    

####################################################################################################
###################################################################################################
################################################################################################




summ = numpy.zeros(maxnode)
summul = numpy.zeros(maxnode)
gitruenode = numpy.zeros(maxnode)
gihatnode = numpy.zeros(maxnode)
for i in range(len(abduse)):
    abdc = abduse[i]
    node = abdc[4]
    summ[node] = summ[node] + abdc[7]
    summul[node] = summul[node] + (abdc[3]-abdc[2])*muhat
    gitruenode[node] = abdc[6]

for i in range(len(summ)):
    if summul[i] != 0:
        gihatnode[i] = summ[i]/2/summul[i] - epsilonhat/muhat


gihatnodem0 = []
nonzeroind = []
for i in range(len(gihatnode)):
    if gihatnode[i] !=0:
        gihatnodem0.append(gihatnode[i])
        nonzeroind.append(1)
    else:
        nonzeroind.append(0)



Qmatrix = 2*spmatrix(1, range(len(gihatnodem0)), range(len(gihatnodem0)))
Rmatrix = -2*matrix(gihatnodem0)
lengi = len(gihatnodem0)


Gmatrix =  spmatrix([1]*(lengi-1) + [-1]*(lengi-1)  , [*range(lengi-1)] + [*range(lengi-1)],\
                    [*range(lengi-1)] + [*range(1,lengi)], (lengi-1,lengi))



Hmatrix = matrix( [0.] * (len(gihatnodem0)-1)  )


opt.solvers.options['show_progress'] = False
sol = qp(Qmatrix, Rmatrix,Gmatrix,Hmatrix)["x"]

gihatconm0=[]
for i in range(len(sol)):
    gihatconm0.append(sol[i]) 



gihatnew = savgol_filter(gihatconm0, 1000, 3) 

giind = 0
gihatcon = []
gitrue = []
for i in range(len(nonzeroind)):
    if nonzeroind[i] !=0:
        gihatcon.append(gihatnew[giind])
        gitrue.append(gitruenode[i])
        giind = giind+1


name = "gihatcon_" + "eps=" + str(epsilon) + "_sam_size=" + str(sam_size) + "constne_len108.dat"
file = open(name,"wb")
pickle.dump(gihatcon,file)
file.close()

name = "gitrue_" + "eps=" + str(epsilon) + "_sam_size=" + str(sam_size) + "constne_len108.dat"
file = open(name,"wb")
pickle.dump(gitrue,file)
file.close()



 
