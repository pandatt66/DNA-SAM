# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 18:18:07 2016

@author:Xinjie Tong
"""

from read_xvg import read_xvg
import sys
import numpy as np

sh_traj = sys.argv[1]
base_traj = sys.argv[2]
dnmp = sys.argv[3]
wids=sys.argv[4]

sh=read_xvg(sh_traj)
base=read_xvg(base_traj)

diff = base[:,1]-sh[:,1]
sbstct=np.copy(sh)
sbstct[:,1]=np.copy(diff)
np.savetxt('pmf_'+dnmp+ '_sbstct.dat', sbstct)

##choose half side of the box:distance d>0 or d <0
numt=len(diff)
dist=[]
time=[]
if abs(min(diff)) > max(diff):
    for i in range (numt):
        if diff[i]>0:
           dist.append(diff[i])
           time.append(sbstct[i,0])               
    hfsys=np.column_stack((time,dist))

elif abs(min(diff)) < max(diff):
    for i in range (numt):
        if diff[i]<0:
           dist.append(diff[i])
           time.append(sbstct[i,0])               
    hfsys=np.column_stack((time,dist))

mindis=min(hfsys[:,1])
mindis_index=np.argmin(dist)

t0=hfsys[mindis_index,0]
data=np.zeros([int(wids),2])
data[0,:]=[t0,mindis]
for i in range(1,int(wids)):
    data[i,1]=(data[i-1,1]+0.1)

nrest_index=[]
def find_nearest(): 
    for i in range(1,int(wids)):
        index =(np.abs(hfsys[:,1]-data[i,1])).argmin()
        (nrest_index).append(index)

find_nearest()

find_time=np.concatenate([np.array([data[0,:]]),(hfsys[nrest_index,:])])
np.savetxt('pmf_'+dnmp+ '_time.dat', find_time)

