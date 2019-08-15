"""
Created on Mon Feb 18 19:58:49 2019

@author: Xinjie
"""
import my_plot_settings_article as mpsa
import numpy as np
import matplotlib.pylab as plt
from read_xvg import read_xvg
import argparse
import uncertainties
import pandas as pd


parser = argparse.ArgumentParser()

parser.add_argument('sysname', help='system name.')
parser.add_argument('-nucs', '--nuc_list', nargs='+',
                    default=['dCMP', 'dTMP', 'dGMP', 'dAMP'],
                    choices=['dCMP', 'dTMP', 'dGMP', 'dAMP'],
                    help='Names of nucleotides. Defalut is dCMP dTMP dGMP dAMP.')

args = parser.parse_args()


sysname = args.sysname
nuc_list = args.nuc_list


vel_output_folder = "./velocities/"

nnucs = len(nuc_list)

df=pd.DataFrame()


for inuc in range(nnucs):
    nuc = nuc_list[inuc]

    ads_vel_upper = np.loadtxt(vel_output_folder+'ads_upper_all_velocity_' + nuc + '_block_error_extrapolation_0.dat')
    ads_vel_lower = np.loadtxt(vel_output_folder+'ads_lower_all_velocity_' + nuc + '_block_error_extrapolation_0.dat')
    des_vel_lower = np.loadtxt(vel_output_folder+'des_all_velocity_' + nuc + '_block_error_extrapolation_0.dat')
    
    x= uncertainties.ufloat(ads_vel_upper[0],ads_vel_upper[2])
    y= uncertainties.ufloat(ads_vel_lower[0],ads_vel_lower[2])

    propagate_error = (x+y)/2

    ads_mean_with_error=[propagate_error.nominal_value, propagate_error.std_dev]
    


    with open(vel_output_folder+'ads_vel_mean_with_propagated_error'+'.dat','ab') as f:
             np.savetxt(f,ads_mean_with_error,header=" ads_vel_mean+-propagated_error" +"||" + \
             nuc)
 
    # round its positive value to two decimal
    def rp(x):
        return abs(round(x,2))

    data = {'ads_mean_vel':rp(ads_mean_with_error[0]),'uncertainty1':rp(ads_mean_with_error[1]),\
            'des_mean_vel':rp(des_vel_lower[0]),'uncertainty2':rp(des_vel_lower[2])}

    df_data = pd.DataFrame(data,index = [nuc])
    df_data['text1']=df_data.ads_mean_vel.astype(str) +'+-'+ df_data.uncertainty1.astype(str)
    df_data['text2']=df_data.des_mean_vel.astype(str) +'+-'+ df_data.uncertainty2.astype(str)
    df = df.append(df_data)

df.to_excel(excel_writer=sysname + "_ads_des_vel_mean_with_error.xlsx")







