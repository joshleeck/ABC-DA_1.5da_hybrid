# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import scipy.linalg
import matplotlib.pyplot as plt
import numpy as np
import math

def Fn_GC_corr(diff, hScale_alpha):
    c = hScale_alpha
    if diff == 0:
        corr = 1.0
    elif abs(diff) <= 2*c:
        if abs(diff) <= c:
            temp5 = -0.25
            temp4 = 0.5
            temp3 = 0.625
            temp2 = -5.0/3.0
            temp1 = 0.0
            temp0 = 1.0
            temp_inv1 = 0.0
        else:
            temp5 = 1.0/12.0
            temp4 = -0.5
            temp3 = 0.625
            temp2 = 5.0/3.0
            temp1 = -5.0
            temp0 = 4.0
            temp_inv1 = -2.0/3.0
            
        adc = abs(diff/c)
        corr = temp5 * adc**5 \
                + temp4 * adc**4 \
                + temp3 * adc**3 \
                + temp2 * adc**2 \
                + temp1 * adc \
                + temp0 * 1.0 \
                + temp_inv1 * 1.0/adc
    else:
        corr=0.0    
    
    return corr

nlongs = 364
hScale_alpha = 250/1.5
Lh = np.zeros((nlongs,nlongs))
Lh_new = np.zeros((nlongs,nlongs))

for i in range(nlongs):
    for j in range(nlongs):
        xdiff1=abs(j-i)
        if xdiff1 <= 0.5*nlongs:
            corr = Fn_GC_corr(xdiff1, hScale_alpha)
            Lh[i,j] = corr
        else:
            if j>i:
                xdiff2 = i + (nlongs-j)
            else:
                xdiff2 = j + (nlongs-i)
            corr = Fn_GC_corr(xdiff2, hScale_alpha)
            Lh[i,j] = corr
            

            
#for extreme localisation
#for i in range(nlongs):
#    for j in range(nlongs):
#        if i==j:
#            Lh[i,j] = 1

#trying to smooth curve
#for i in range(nlongs):
#    min_val = np.mean(Lh[i,:])/2
#    for j in range(nlongs):
#       if Lh[i,j] < min_val:
#            Lh[i,j] = min_val

#half width
smooth = 0
for i in range(nlongs):
    for j in range(nlongs):
        v = 1 
        Lh_new[i,j] += Lh[i,j]
        while v < smooth+1:
            Lh_new[i,j] += Lh[i,(j-v)%nlongs] + Lh[i,(j+v)%nlongs]
            v+=1
        Lh_new[i,j] = round(Lh_new[i,j]/(2*smooth+1),8)
                            
#checks for symmetry
for i in range(nlongs):
    for j in range(nlongs):
        if Lh_new[i,j] != Lh_new[j,i]:
            print(Lh_new[i,j], Lh_new[j,i])


ev, evec = scipy.linalg.eigh(Lh_new)
#print(ev)

#check negative ev
negative=False
for i in ev:
    if i <=0:
        negative=True
if negative:
    print('Got negative values')

sum1=sum(ev)

#floor eigenvalues
#for i in range(len(ev)):
#    if ev[i]<0:
#       ev[i] = 0

#no relationship between scale and magnitude of min ev    
#print(np.min(ev))
sum2=sum(ev)
#print(sum(ev))

#restore initial total variance
#fac=sum1/sum2
#ev = ev*fac

ev_mat = np.zeros((nlongs,nlongs), dtype=complex)
for i in range(nlongs):
    ev_mat[i,i] = ev[i]

tmp = np.matmul(ev_mat,np.transpose(evec),dtype=complex)
tmp2 = np.matmul(evec,tmp,dtype=complex)
#print(tmp2)


plt.plot(Lh_new[50,:], label='original')
plt.plot(tmp2[50,:], label='reconstructed')
plt.title('hScale_alpha=250km')
plt.ylim(0,1.1)
plt.ylabel('Correlation')
plt.xlabel('Longitudinal gridpoint index')
plt.axhline(y=1,color='grey',linestyle='--')
plt.legend()
plt.savefig('loc1',dpi=300)
plt.show()


#print(evec)