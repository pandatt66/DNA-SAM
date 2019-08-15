# -*- coding: utf-8 -*-
"""
Created on Sat Jan  2 13:09:12 2016

@author: Xinjie Tong
"""
#from read_xvg import read_xvg
import sys
import numpy as np


coord_fixed = sys.argv[1]
dnmp = sys.argv[2]

coord = np.loadtxt(coord_fixed)

coord2=coord[:,:2]

countrow= len(coord2)-1

new=np.empty((countrow,2))
for i in range(countrow):
    a=0.5*(coord2[i,0]+coord2[i+1,0])
    b=((coord2[i+1,1]-coord2[i,1])/(coord2[i+1,0]-coord2[i,0]))
    new[i,:]= (a*1000,b)
np.savetxt(dnmp+'_1_veloc_traj.xvg', new)
