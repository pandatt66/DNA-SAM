# -*- coding: utf-8 -*-
"""
Created on Thu May 12 11:09:49 2016

@author: Xinjie Tong
"""
##The default is to use all nucleotides. For all nucleotides use the script like this:

##python plot.py pmf_shift dat pmf_plot

#For a subset of the nucleotides, use the script like this:

#run plot.py pmf_shift dat pmf_plot cooch3

import my_plot_settings_article as mpsa
import numpy as np
import matplotlib.pylab as plt
from read_xvg import read_xvg
import argparse
from matplotlib.ticker import MultipleLocator, FormatStrFormatter



#1.83 - 0.388
#cooch3 2.05 - 0.388
#poxy 2.18 - 0.388
def plot_pmfs(infile_prefix, infile_extension,sysname, nuc_list, outfile_prefix):

    nnucs = len(nuc_list)

    for inuc in range(nnucs):

        nuc = nuc_list[inuc]

        infile = infile_prefix + '_' + nuc +'_'+ sysname + '.' + infile_extension

        data=read_xvg(infile)
        ax = plt.subplot(111)
        if sysname == 'ch3':
            dist = 1.83-0.388
        elif sysname == 'cooch3':
            dist = 2.05 - 0.388
        else:
            dist = 2.18 - 0.388

        if nuc == 'dCMP':
            plt.plot(data[:, 0] -dist, data[:, 1], 'g', label= nuc)
        elif nuc == 'dTMP':
            plt.plot(data[:, 0] -dist , data[:, 1], 'b', label= nuc)
        elif nuc == 'dGMP':
            plt.plot(data[:, 0] -dist, data[:, 1], 'm', label= nuc)
        else:
            plt.plot(data[:, 0] -dist, data[:, 1], 'r', label= nuc)
     #   print(data[np.argmin(data[:,1]), 0]-dist)
    mpsa.axis_setup('x')
    mpsa.axis_setup('y')
    plt.xlim(0.0, 0.8)
    plt.ylim(-8, 10)
    plt.xlabel("$d_w$ (nm)", labelpad= mpsa.axeslabelpad)
    plt.ylabel('PMF ($ k_B$T)', labelpad= mpsa.axeslabelpad)
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

    mpsa.save_figure(outfile_prefix + '_'+sysname +'.png')
    plt.show()



def plot_pmfs_SH(infile_prefix, infile_extension,sysname, nuc_list, outfile_prefix):

    nnucs = len(nuc_list)

    for inuc in range(nnucs):
        ax.yaxis.set_major_locator(ymajorLocator)
        ax.yaxis.set_minor_locator(yminorLocator)

        nuc = nuc_list[inuc]

        infile = infile_prefix + '_' + nuc +'_'+ sysname + '.' + infile_extension

        data=read_xvg(infile)

        if nuc == 'dCMP':
            plt.plot(data[:, 0], data[:, 1], 'g', label=sysname +'_' + nuc)
        elif nuc == 'dTMP':
            plt.plot(data[:, 0], data[:, 1], 'b', label=sysname +'_' + nuc)
        elif nuc == 'dGMP':
            plt.plot(data[:, 0], data[:, 1], 'm', label=sysname +'_' + nuc)
        else:
            plt.plot(data[:, 0], data[:, 1], 'r', label=sysname +'_' + nuc)
    mpsa.axis_setup('x')
    mpsa.axis_setup('y')
    plt.xlabel("SH-base Distance(nm)", labelpad= mpsa.axeslabelpad)
    plt.ylabel('PMF(kT)', labelpad= mpsa.axeslabelpad)
    plt.legend()
    plt.savefig(outfile_prefix +'_'+sysname+ '.png', bbox_inches ='tight', pad_inches = 0.02, dpi =300)
    plt.show()


def plot_density(infile_prefix, infile_extension, nuc_list, outfile_prefix):

    nnucs = len(nuc_list)

    for inuc in range(nnucs):

        nuc = nuc_list[inuc]

        infile = infile_prefix + '_' + nuc + '.' + infile_extension

        data=read_xvg(infile)

        if nuc == 'dCMP':
            plt.plot(data[:, 0], data[:, 1], 'g', label='dCMP')
        elif nuc == 'dTMP':
            plt.plot(data[:, 0], data[:, 1], 'b', label='dTMP')
        elif nuc == 'dGMP':
            plt.plot(data[:, 0], data[:, 1], 'm', label='dGMP')
        else:
            plt.plot(data[:, 0], data[:, 1], 'r', label='dAMP')

    plt.xlabel("zcoord(nm)")
    plt.ylabel('suface-atom density(kg/m^3)')
    plt.legend()
    plt.savefig(outfile_prefix + '.png')
    plt.show()

#run plot.py 1_assoa_numcont_fixed dat assoa_dCMP_plot -nucs dCMP
def plot_assonums(infile_prefix, infile_extension, nuc_list, outfile_prefix):
    nnucs = len(nuc_list)
    for inuc in range(nnucs):
        nuc = nuc_list[inuc]

        infile = nuc + '_' + infile_prefix  +  '.' + infile_extension

        data=read_xvg(infile)

        if nuc == 'dCMP':
            plt.plot(data[:, 0]/1000, data[:, 1], 'g.', label='dCMP')
            avg = np.mean(data[:,1],axis=0)
        elif nuc == 'dTMP':
            plt.plot(data[:, 0]/1000, data[:, 1], 'b.', label='dTMP')
            avg=np.mean(data[:,1],axis=0)
        elif nuc == 'dGMP':
            plt.plot(data[:, 0]/1000, data[:, 1], 'm.', label='dGMP')
            avg=np.mean(data[:,1],axis=0)
        else:
            plt.plot(data[:, 0]/1000, data[:, 1], 'r.', label='dAMP')
            avg=np.mean(data[:,1],axis=0)
    plt.xlabel("tims(ns)")
    plt.ylabel('Association numbers')
    plt.legend()
    plt.savefig(outfile_prefix + '.png')
    plt.show()

    print(avg)
    return avg

    #run plot.py 1_coord_fixed dat xcoord poxy -nucs dCMP
def plot_xcoord(infile_prefix, infile_extension, nuc_list, outfile_prefix,sysname):

    nnucs = len(nuc_list)
    for inuc in range(nnucs):
        nuc = nuc_list[inuc]

        infile = nuc + '_' + infile_prefix  +  '.' + infile_extension

        data=read_xvg(infile)
        ax = plt.subplot(111)
        if nuc == 'dCMP':
            plt.plot(data[:, 0], data[:, 1], 'g', label='dCMP')
        elif nuc == 'dTMP':
            plt.plot(data[:, 0], data[:, 1], 'b', label='dTMP')
        elif nuc == 'dGMP':
            plt.plot(data[:, 0], data[:, 1], 'm', label='dGMP')
        else:
            plt.plot(data[:, 0], data[:, 1], 'r', label='dAMP')


    plt.xlim(0, 2000)
   # plt.ylim(-5000, 0)
    xmajorLocator = MultipleLocator(500)
    xminorLocator = MultipleLocator(100)

    ymajorLocator = MultipleLocator(1000)
    yminorLocator = MultipleLocator(500)

    ax.xaxis.set_major_locator(xmajorLocator)
    ax.yaxis.set_major_locator(ymajorLocator)
    ax.yaxis.set_minor_locator(yminorLoplot_zcoordcator)


    plt.legend()
    plt.xlabel("Time (ns)")
    plt.ylabel('x (nm)')
    plt.legend()
    mpsa.save_figure(outfile_prefix +'_' +sysname + '.png')
    plt.show()

#run plot.py 1_coord_fixed dat zcoord -nucs dAMP

def plot_zcoord(infile_prefix, infile_extension, nuc_list, outfile_prefix):

    nnucs = len(nuc_list)
    for inuc in range(nnucs):
        nuc = nuc_list[inuc]

        infile = nuc + '_' + infile_prefix  +  '.' + infile_extension

        data=read_xvg(infile)

        def plot_zcoord_setting(outfile_prefix):
            ax = plt.subplot(111)
            plt.xlim(0, 1888)
           # plt.ylim(-1.0, 1.2)

            plt.xlabel("Time (ns)")
            plt.ylabel('z (nm)')

            xmajorLocator = MultipleLocator(500)
           # xminorLocator = MultipleLocator(1)

            ymajorLocator = MultipleLocator(0.5)
            #yminorLocator = MultipleLocator(0.1)

            ax.xaxis.set_major_locator(xmajorLocator)
           # ax.xaxis.set_minor_locator(xminorLocator)

            ax.yaxis.set_major_locator(ymajorLocator)
            #ax.yaxis.set_minor_locator(yminorLocator)


           # plt.legend()
            mpsa.save_figure(outfile_prefix + '_'+nuc +'.png')
            plt.close()


        if nuc == 'dCMP':
            plt.plot(data[:, 0], data[:, 3], 'g', label='dCMP')
            plot_zcoord_setting(outfile_prefix)
        elif nuc == 'dTMP':
            plt.plot(data[:, 0], data[:, 3], 'b', label='dTMP')
            plot_zcoord_setting(outfile_prefix)
        elif nuc == 'dGMP':
            plt.plot(data[:, 0], data[:, 3], 'm', label='dGMP')
            plot_zcoord_setting(outfile_prefix)
        else:
            plt.plot(data[:, 0], data[:, 3], 'r', label='dAMP')
            plot_zcoord_setting(outfile_prefix)


#run plot.py hbnum_cooch3 xvg hbnum_dAMP_plot -nucs dAMP
def plot_hbnums(infile_prefix, infile_extension, nuc_list, outfile_prefix):

    nnucs = len(nuc_list)
    for inuc in range(nnucs):
        nuc = nuc_list[inuc]

        infile = infile_prefix + '_' + nuc + '.' + infile_extension

        data=read_xvg(infile)

        if nuc == 'dCMP':
            plt.plot(data[:, 0]/1000, data[:, 1], 'g.', label='dCMP')
            avg = np.mean(data[:,1],axis=0)
        elif nuc == 'dTMP':
            plt.plot(data[:, 0]/1000, data[:, 1], 'b.', label='dTMP')
            avg = np.mean(data[:,1],axis=0)
        elif nuc == 'dGMP':
            plt.plot(data[:, 0]/1000, data[:, 1], 'm.', label='dGMP')
            avg = np.mean(data[:,1],axis=0)
        else:
            plt.plot(data[:, 0]/1000, data[:, 1], 'r.', label='dAMP')
            avg = np.mean(data[:,1],axis=0)

    plt.xlabel("tims(ns)")
    plt.ylabel('Hydrogen bond numbers')
    plt.legend()
    plt.savefig(outfile_prefix + '.png')
    plt.show()
    print(avg)
############################################################
#run plot.py hbnum_cooch3_dw_dist xvg hbnum_distr_plot
###
def plot_dw_distr_hbnums(infile_prefix, infile_extension, nuc_list, outfile_prefix):

    nnucs = len(nuc_list)
    for inuc in range(nnucs):
        nuc = nuc_list[inuc]

        infile = infile_prefix + '_' + nuc + '.' + infile_extension

        data=read_xvg(infile)

        if nuc == 'dCMP':
            plt.plot(data[:, 0], data[:, 1], 'g', label='dCMP')
        elif nuc == 'dTMP':
            plt.plot(data[:, 0], data[:, 1], 'b', label='dTMP')
        elif nuc == 'dGMP':
            plt.plot(data[:, 0], data[:, 1], 'm', label='dGMP')
        else:
            plt.plot(data[:, 0], data[:, 1], 'r', label='dAMP')

    plt.xlabel("Dw(nm))")
    plt.ylabel('Distance distribution of H-bond')
    plt.legend(loc='upper left')
    plt.savefig(outfile_prefix + '.png')
    plt.show()


if __name__ == "__main__":

    # Command line arguments
    parser = argparse.ArgumentParser()

    parser.add_argument('infile_prefix', help='Prefix for input file name.')
    parser.add_argument('infile_extension', help='Extension for input file name.')
    parser.add_argument('outfile_prefix', help='Prefix for output file name.')
    parser.add_argument('sysname', help='system name.')
    parser.add_argument('-nucs', '--nuc_list', nargs='+',
                        default=['dAMP', 'dCMP', 'dGMP', 'dTMP'],
                        choices=['dAMP', 'dCMP', 'dGMP', 'dTMP'],
                        help='Names of nucleotides. Defalut is dAMP dCMP dGMP dTMP.')

    args = parser.parse_args()

    infile_prefix = args.infile_prefix
    infile_extension = args.infile_extension
    outfile_prefix = args.outfile_prefix
    sysname = args.sysname
    nuc_list = args.nuc_list

    #plot_pmfs(infile_prefix, infile_extension,sysname, nuc_list, outfile_prefix)
    #plot_assonums(infile_prefix, infile_extension, nuc_list, outfile_prefix)
    #plot_xcoord(infile_prefix, infile_extension, nuc_list, outfile_prefix)
   # plot_xcoord(infile_prefix, infile_extension, nuc_list, outfile_prefix,sysname)
    #plot_zcoord(infile_prefix, infile_extension, nuc_list, outfile_prefix)
    #plot_density(infile_prefix, infile_extension, nuc_list, outfile_prefix)
    #plot_hbnums(infile_prefix, infile_extension, nuc_list, outfile_prefix)
   # plot_dw_distr_hbnums(infile_prefix, infile_extension, nuc_list, outfile_prefix)
