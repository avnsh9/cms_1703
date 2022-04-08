import numpy as np 
import matplotlib.pyplot as plt
from sympy import *
import math
from tqdm import tqdm
import concurrent.futures as cf
import time

#################
# initialisation of values
signal_temp=np.loadtxt('./data1703.dat', usecols=(2),skiprows=1)
events=np.loadtxt('./data1703.dat', usecols=0,skiprows=1)
bkg=np.loadtxt('./data1703.dat', usecols=1,skiprows=1)
event=np.loadtxt('./data1703.dat', dtype=int, usecols=0,skiprows=1)
bin_size=np.loadtxt('./data1703.dat', usecols=3,skiprows=1)
bkg_err=np.loadtxt('./data1703.dat', usecols=4,skiprows=1)

signal=bin_size*signal_temp

a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v=symbols('theta_1 theta_2 theta_3 theta_4 theta_5 theta_6 theta_7 theta_8 theta_9 theta_10 theta_11 theta_12 theta_13 theta_14 theta_15 theta_16 theta_17 theta_18 theta_19 theta_20 theta_21 theta_22')

theta=Matrix([a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v])

#######Initialising matrix########

corr_mat=np.loadtxt('./correlation_MONOJ.dat',usecols=range(22))
corr_mat=np.flipud(corr_mat)


#########Getting covariance matirx from correlation matrix############

std_dev_mat=np.diag(bkg_err)

covar_mat_temp=np.matmul(std_dev_mat,corr_mat)
covar_mat=np.matmul(covar_mat_temp,std_dev_mat)

covar_mat=Matrix(covar_mat)
covar_mat_inv=covar_mat**(-1)

#generating F[\theta]
uu=symbols('\mu')
kk=Matrix(uu*signal+bkg)+theta
kk=diag(kk[0],kk[1],kk[2],kk[3],kk[4],kk[5],kk[6],kk[7],kk[8],kk[9],kk[10],kk[11],kk[12],kk[13],kk[14],kk[15],kk[16],kk[17],kk[18],kk[19],kk[20],kk[21])

# F[\theta]
F_theta=kk*(covar_mat_inv*theta+Matrix([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]))-Matrix(events)

#Jacobian for F[\theta]
JF_theta=F_theta.jacobian(theta)

# initial guess for solution
guess=Matrix([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])



def calculate(jj):
    ll=jj*0.01
    guess=Matrix([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    
    for ii in range(10):
        nn=(JF_theta.subs({uu:-2.00+ll,a:guess[0],b:guess[1],c:guess[2],d:guess[3],e:guess[4],f:guess[5],g:guess[6],h:guess[7],i:guess[8],j:guess[9],k:guess[10],l:guess[11],m:guess[12],n:guess[13],o:guess[14],p:guess[15],q:guess[16],r:guess[17],s:guess[18],t:guess[19],u:guess[20],v:guess[21]})**-1)*(F_theta.subs({uu:-2.00+ll,a:guess[0],b:guess[1],c:guess[2],d:guess[3],e:guess[4],f:guess[5],g:guess[6],h:guess[7],i:guess[8],j:guess[9],k:guess[10],l:guess[11],m:guess[12],n:guess[13],o:guess[14],p:guess[15],q:guess[16],r:guess[17],s:guess[18],t:guess[19],u:guess[20],v:guess[21]}))
        guess=guess-nn

    return jj,[i for i in guess]
    



with cf.ProcessPoolExecutor() as exe:
    results = [exe.submit(calculate, jj) for jj in tqdm(range(231))]
with open('thetas.dat', 'w') as b:
    for f in cf.as_completed(results):
        b.write(str(f.result()[0])+' '+str(f.result()[1])+'\n')
        
        





