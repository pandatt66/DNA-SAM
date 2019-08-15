# -*- coding: utf-8 -*-
"""
Created on Mon May  9 12:24:41 2016

@author:Xinjie Tong
"""

import my_stats
import numpy as np
import sys

#states could be ads, des,ads_upper, ads_lower,des_upper, des_lower
#Should be  consist with input data file
veldata=sys.argv[1]
dnmp=sys.argv[2]
states=sys.argv[3]


data=np.loadtxt(veldata)

mean_vel, se=my_stats.wghtdunc_variance(data[:,0],data[:,1])

np.savetxt(dnmp+'_mean_'+ states+'.dat',(mean_vel,se),header="mean_vel|se")
