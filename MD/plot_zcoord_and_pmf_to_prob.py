# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 17:02:59 2016

@author: Xinjie Tong
"""
#Compare equilibrium free energy  with nonequilibrium one, converted to frequency#
import my_plot_settings_article as mpsa
import numpy as np
import matplotlib.pyplot as plt
import argparse
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


def cmpr_noneq_eq(infile_extension,sysname, nuc_list, outfile_prefix):

    nnucs = len(nuc_list)

    for inuc in range(nnucs):

        nuc = nuc_list[inuc]

        infile = 'pmf_shift' + '_' + nuc +'_'+ sysname + '.' + infile_extension
        infile2 =nuc + '_zcoord_to_prob.dat'
        nonequil_zcoord = np.loadtxt(infile2)

        shift_pmf = np.loadtxt(infile)

        energy = shift_pmf[:,1]
        dw = shift_pmf[:,0] - dsurf +r_probe
        Q=np.sum(np.exp(-energy))

        ads_prob= np.exp(-energy)/Q
        density = ads_prob/(dw[20]-dw[19])
    #    print(len(ads_prob))

        ax = plt.subplot(111)

        if nuc == 'dCMP':
            #plt.plot(dw,ads_prob,'g',label='equil')
            plt.plot(dw,density,'g',label='equil')
            plt.plot(nonequil_zcoord[:,0], nonequil_zcoord[:,1], 'k', label='non-equil')
            if sysname == 'cooch3':
                smooth= np.loadtxt('dCMP_zcoord_to_prob_smoothcooch3.dat')
                plt.plot(smooth[:,0],smooth[:,1],'c--',label='smooth')
        elif nuc == 'dTMP':
            #plt.plot(dw,ads_prob,'b',label='equil')
            plt.plot(dw,density,'b',label='equil')
            plt.plot(nonequil_zcoord[:,0], nonequil_zcoord[:,1], 'k', label='non-equil')
            if sysname == 'cooch3':
                smooth= np.loadtxt('dTMP_zcoord_to_prob_smoothcooch3.dat')
                plt.plot(smooth[:,0],smooth[:,1],'c--',label='smooth')
        elif nuc == 'dGMP':
            #plt.plot(dw,ads_prob,'m',label='equil')
            plt.plot(dw,density,'m',label='equil')
            plt.plot(nonequil_zcoord[:,0], nonequil_zcoord[:,1], 'k', label='non-equil')
            if sysname == 'cooch3':
                smooth= np.loadtxt('dGMP_zcoord_to_prob_smoothcooch3.dat')
                plt.plot(smooth[:,0],smooth[:,1],'c--',label='smooth')
        else:
            #plt.plot(dw,ads_prob,'r',label='equil')
            plt.plot(dw,density,'r',label='equil')
            plt.plot(nonequil_zcoord[:,0], nonequil_zcoord[:,1], 'k', label='non-equil')
            if sysname == 'cooch3':
                smooth= np.loadtxt('dAMP_zcoord_to_prob_smoothcooch3.dat')
                plt.plot(smooth[:,0],smooth[:,1],'c--',label='smooth')


        plt.xlim(0.0, 0.6)

        if sysname == 'ch3':
          plt.ylim(0, 12)
          ymajor =4
          yminor = 1
        elif sysname == 'cooch3':
          plt.ylim(0, 2.8)
          ymajor =0.5
          yminor = 0.1
        else:
          plt.ylim(0, 18)
          ymajor =3
          yminor = 1




        xmajorLocator = MultipleLocator(0.2)
        xminorLocator = MultipleLocator(0.05)

        ymajorLocator = MultipleLocator(ymajor)
        yminorLocator = MultipleLocator(yminor)

        ax.xaxis.set_major_locator(xmajorLocator)
        ax.xaxis.set_minor_locator(xminorLocator)

        ax.yaxis.set_major_locator(ymajorLocator)
        ax.yaxis.set_minor_locator(yminorLocator)



        plt.xlabel("$d_w $(nm)")
        plt.ylabel('Pd ($nm^{-1}$)')
        plt.legend()
        mpsa.save_figure(nuc +'_'+ outfile_prefix +'_'+sysname+ '.png')
        plt.close()


if __name__ == "__main__":

    # Command line arguments
    parser = argparse.ArgumentParser()

    parser.add_argument('infile_extension', help='Extension for input file name.')
    parser.add_argument('outfile_prefix', help='Prefix for output file name.')
    parser.add_argument('sysname', help='system name.')
    parser.add_argument('-nucs', '--nuc_list', nargs='+',
                        default=['dCMP', 'dTMP', 'dGMP', 'dAMP'],
                        choices=['dCMP', 'dTMP', 'dGMP', 'dAMP'],
                        help='Names of nucleotides. Defalut is dCMP dTMP dGMP dAMP.')

    args = parser.parse_args()

    infile_extension = args.infile_extension
    outfile_prefix = args.outfile_prefix
    sysname = args.sysname
    nuc_list = args.nuc_list

    r_probe = 0.388

    if sysname == 'ch3':
       dsurf = 1.83
    elif sysname == 'cooch3':
       dsurf = 2.05
    else:
       dsurf = 2.18


    cmpr_noneq_eq(infile_extension,sysname, nuc_list, outfile_prefix)
