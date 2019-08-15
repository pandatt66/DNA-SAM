# -*- coding: utf-8 -*-
"""
Created on Sat Feb 13 12:30:15 2016

@author: xtong3
"""
#run shift_pmf.py 1pca_dopc_profile.xvg 'PCA in DOPC bilayer' pmf-pca-dopc

import sys
import numpy as np
from read_xvg import read_xvg
import matplotlib.pylab as plt
import my_plot_settings_article as mpsa

def substract(traj1,traj2):
    lower=read_xvg(traj1)
    upper=read_xvg(traj2)
    diff = upper[:,1]-lower[:,1]
    sbstct=np.copy(lower)
    sbstct[:,1]=np.copy(diff)
    return sbstct

if __name__ == "__main__":

    profile=sys.argv[1]
    outfile_prefix = sys.argv[2]
    #traj1=sys.argv[3]
    #traj2=sys.argv[4]

    kbT=read_xvg(profile)[:195,1]
    wids= read_xvg(profile)[:195,0]
    j=-1
    for i in wids >=3.4:
        if i == False:

	        j+=1
        ave=np.mean(kbT[j:])
        kbT[:]=kbT[:]-ave
        data=np.column_stack((wids,kbT))
        np.savetxt('pmf_shift_'+outfile_prefix+ '.dat', data)

#sbstct=substract(traj1,traj2)
#ave_dist=np.mean(sbstct[100:,1])

#np.savetxt(dnmp+'_edgrp_SH_'+'dist'+ '.dat',ave_dist.reshape(1,))

    plt.plot(wids, kbT, 'r',linewidth=2)
    plt.axhline(y=0, ls='--')
    plt.xlim(0, 4.2)
    plt.ylim(-25, 10)

    mpsa.axis_setup('x')
    mpsa.axis_setup('y')
    plt.xlabel("$d_w$ (nm)", labelpad=mpsa.axeslabelpad)
	#plt.tick_params(labelsize=15)
    plt.ylabel('PMF ($k_B$T)',labelpad=mpsa.axeslabelpad)
    mpsa.save_figure(outfile_prefix + '.png')
  #  plt.show()
