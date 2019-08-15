# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 19:58:49 2019

@author: Xinjie
"""
import my_plot_settings_article as mpsa
import numpy as np
import matplotlib.pylab as plt
from read_xvg import read_xvg
import argparse
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import sys
from scipy import integrate

#run plot_pmf_same_dnmp.py pmf_shift dat pmf_all dGMP

def plot_pmfs_same_dnmp(infile_prefix, infile_extension,nuc_name, outfile_prefix,sys_list):
    
    nsysnames = len(sys_list)
    for isysname in range(nsysnames):
    
        
        sysname = sys_list[isysname]
        
        infile = infile_prefix + '_' + nuc_name +'_'+ sysname + '.' + infile_extension

        data=read_xvg(infile)
        kbT = data[:, 1]
        ax = plt.subplot(111)
        if sysname == 'ch3':
                surf = data[:, 0]-1.83 + 0.388
                plt.plot(surf, kbT, 'g', label='CH3-SAM')
                with open('pmf_table_ch3.dat','ab') as f:
            
                    pmf_min = np.min(kbT)
                    pmf_min_idx = np.argmin(kbT[:-10])
                    surf_idx = int(np.argwhere(surf>=0)[0])
                    surf_idx_end = int(np.argwhere(surf>=0.8)[0])
                    Fads = integrate.simps(kbT[surf_idx:surf_idx_end],surf[surf_idx:surf_idx_end])
                    dw_max = surf[pmf_min_idx:surf_idx_end][np.argmax(kbT[pmf_min_idx:surf_idx_end])]

                    table = np.column_stack((pmf_min,Fads,dw_max))
                    np.savetxt(f,table,header = "pmf_min|Fads|dw_max||--->" +nuc_name)
                
        
        
        elif sysname == 'cooch3':
                surf = data[:, 0]-2.05 + 0.388
                plt.plot(surf, data[:, 1], 'b', label='COOCH3-SAM' )
                with open('pmf_table_cooch3.dat','ab') as f:
            
                    pmf_min = np.min(kbT)
                    pmf_min_idx = np.argmin(kbT[:-10])
                    surf_idx = int(np.argwhere(surf>=0)[0])
                    surf_idx_end = int(np.argwhere(surf>=0.8)[0])
                    Fads = integrate.simps(kbT[surf_idx:surf_idx_end],surf[surf_idx:surf_idx_end])
                    dw_max = surf[pmf_min_idx:surf_idx_end][np.argmax(kbT[pmf_min_idx:surf_idx_end])]

                    table = np.column_stack((pmf_min,Fads,dw_max))
                    np.savetxt(f,table,header = "pmf_min|Fads|dw_max||--->" +nuc_name)
                
                
        else:
                surf = data[:, 0]-2.18 + 0.388
                plt.plot(surf, data[:, 1], 'r', label='OC6H5-SAM')
                with open('pmf_table_poxy.dat','ab') as f:
            
                    pmf_min = np.min(kbT)
                    pmf_min_idx = np.argmin(kbT[:-10])
                    surf_idx = int(np.argwhere(surf>=0)[0])
                    surf_idx_end = int(np.argwhere(surf>=0.8)[0])
                    Fads = integrate.simps(kbT[surf_idx:surf_idx_end],surf[surf_idx:surf_idx_end])
                    dw_max = surf[pmf_min_idx:surf_idx_end][np.argmax(kbT[pmf_min_idx:surf_idx_end])]

                    table = np.column_stack((pmf_min,Fads,dw_max))
                    np.savetxt(f,table,header = "pmf_min|Fads|dw_max||--->" +nuc_name)
            
            
            
            
    plt.xlim(0.0, 0.8)
    plt.ylim(-8, 10)        
            
    #mpsa.axis_setup('x')
   # mpsa.axis_setup('y')    
    plt.xlabel("$d_w$ (nm)")
    plt.ylabel('PMF ($ k_B$T)')
    plt.legend()
    xmajorLocator = MultipleLocator(0.2)
    xminorLocator = MultipleLocator(0.05)

    ymajorLocator = MultipleLocator(4)
    yminorLocator = MultipleLocator(1)

    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_minor_locator(xminorLocator)

    ax.yaxis.set_major_locator(ymajorLocator)
    ax.yaxis.set_minor_locator(yminorLocator)
  #  plt.rcParams['xtick.labelsize']=12
   # plt.rcParams['ytick.labelsize']=12
    
    mpsa.save_figure(outfile_prefix + '_'+nuc_name +'.png')
  #  plt.show()
    
if __name__ == "__main__":

    # Command line arguments
    parser = argparse.ArgumentParser()
    
    parser.add_argument('infile_prefix', help='Prefix for input file name.')
    parser.add_argument('infile_extension', help='Extension for input file name.')
    parser.add_argument('outfile_prefix', help='Prefix for output file name.')
    parser.add_argument('nuc_name')
    parser.add_argument('-sysname', '--sys_list', nargs='+', 
                        default=['ch3', 'cooch3', 'poxy'], 
                        choices=['ch3', 'cooch3', 'poxy'],
                        help='sysname of SAMs. Defalut is ch3,cooch3,poxy.')
    
    args = parser.parse_args()
    
    infile_prefix = args.infile_prefix
    infile_extension = args.infile_extension
    outfile_prefix = args.outfile_prefix
    nuc_name = args.nuc_name
    sys_list = args.sys_list
    
    plot_pmfs_same_dnmp(infile_prefix, infile_extension,nuc_name, outfile_prefix,sys_list)
    
