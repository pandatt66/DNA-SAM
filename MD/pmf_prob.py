# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 17:02:59 2016

@author: Xinjie Tong
"""

import numpy as np
import sys
from shift_pmf import substract
import matplotlib.pyplot as plt

shift_pmf=sys.argv[1]
pmf_base=sys.argv[2]
pmf_SH1=sys.argv[3]
pmf_SH2=sys.argv[4]
dnmp=sys.argv[5]

org_data=np.loadtxt(shift_pmf,float)
lines =len(org_data)
start=org_data[(lines-1),:]
i=start[0]+0.01
add_wids=np.arange(i,2.76,0.01)
c=np.zeros(len(add_wids))
add_data=np.column_stack((add_wids,c))
all_data=np.vstack((org_data,add_data))
Q=np.sum(np.exp(-all_data[:,1]))   
pmf_prob=np.exp(-all_data[:,1])/Q

    
traj_data1=substract(pmf_SH1,pmf_base)
traj_data2=substract(pmf_SH2,pmf_base)
dist1=traj_data1[:,1]
dist2=traj_data2[:,1]
min_base_SH=[]
a=0
while a<len(traj_data1):
    min_base_SH.append(min(abs(dist1[a]), abs(dist2[a])))
    a+=1
exp=np.histogram(min_base_SH,bins=200,range=(1.51,2.76))

plt.plot(exp[1][:-1],exp[0]/float(len(traj_data1)),label='non_equil')
plt.plot(all_data[:,0],pmf_prob,label='pmf')
plt.xlabel('Distance from dNMP rings to Sulfur(nm)')
plt.ylabel('Probability')
plt.title(dnmp)
plt.legend()
plt.show()
plt.savefig(dnmp+'_pmf_compare.tif')

#j2 = [a for a in dist1 if abs(a) <= 0.35]
#j3=  [b for b in dist2 if abs(b) <= 0.35]
#
#ads=float(len(j2)+len(j3))
#all_mv=float(len(traj_data1)+len(traj_data2))
#exp_prob= ads/all_mv
#np.savetxt(dnmp+'_pmf_'+'prob'+ '.dat',(pmf_prob.reshape(1,),exp_prob),header="pmf_prob|non-equil_prob")
