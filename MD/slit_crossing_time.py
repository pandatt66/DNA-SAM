# -*- coding: utf-8 -*-
"""
Created on Mon July 22 19:58:49 2019

@author: Xinjie
"""
import my_plot_settings_article as mpsa
import numpy as np
import matplotlib.pylab as plt
from read_xvg import read_xvg
import argparse
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import sys
from scipy import integrate,stats
import pandas as pd



#Output: slit crossing time, adsorption percentage of time are saved in .xlsx


def slit_cross_time(infile_prefix,infile_extension,outfile_prefix, sysname,nuc_list):

    df = pd.DataFrame()
    nnucs = len(nuc_list)
    for inuc in range(nnucs):
        nuc = nuc_list[inuc]

        infile = nuc +'_' + infile_prefix  + '.' + infile_extension

# Throw 10ns for equilibrium.
        data=read_xvg(infile)[99:]
        zcoord = data[:,3]
        xcoord = data[:,1]
        time = data[:,0]

        i = 0
      # record each time when dNMP cross the cutoff distance
        ads_cross_time =[]
        #trunc_ndx=[]
        found = False

        if zcoord[0]>=cutoff:
           n=1
           while i < len(zcoord)-1:
                 if n%2!=0:
                    if zcoord[i]<=cutoff*(-1)**n:
                       ads_cross_time.append(time[i])
                       if found == False:
                           zcoord_trunc = zcoord[i:]
                           trunc_ndx=i
                           print(trunc_ndx)

                           time_trunc = time[i:]
                           found = True
                       i+=1
                       n+=1
                    else:
                       i+=1
                 elif n%2==0:
                    if zcoord[i]>=cutoff*(-1)**n:
                       ads_cross_time.append(time[i])
                       i+=1
                       n+=1
                    else:
                       i+=1



        elif zcoord[0]<=-cutoff:
             n=2
             while i < len(zcoord)-1:
                   if n%2==0:
                        if zcoord[i]>=cutoff*(-1)**n:
                           ads_cross_time.append(time[i])
                           if found == False:
                              zcoord_trunc = zcoord[i:]
                              trunc_ndx=i
                              print(trunc_ndx)
                              time_trunc = time[i:]
                              found = True
                           i+=1
                           n+=1
                        else:
                           i+=1
                   elif n%2!=0:
                        if zcoord[i]<=cutoff*(-1)**n:
                           ads_cross_time.append(time[i])
                           i+=1
                           n+=1
                        else:
                           i+=1

        else:
               idx1=np.where(zcoord>=cutoff)[0][0]
               idx2= np.where(zcoord<=-cutoff)[0][0]
               if idx1>idx2:
                  ads_cross_time.append(time[idx1])
                  if found == False:
                     zcoord_trunc = zcoord[idx1:]
                     trunc_ndx=idx1
                     print(trunc_ndx)
                     time_trunc = time[idx1:]
                     found = True

                  n=2

                  while i < len(zcoord_trunc)-1:
                         if n%2==0:
                            if zcoord_trunc[i]>=cutoff*(-1)**n:
                               ads_cross_time.append(time_trunc[i])
                               i+=1
                               n+=1
                            else:
                               i+=1
                         elif n%2!=0:
                            if zcoord_trunc[i]<=cutoff*(-1)**n:
                               ads_cross_time.append(time_trunc[i])
                               i+=1
                               n+=1
                            else:
                               i+=1
             #     ads_time= (len(zcoord_trunc[np.where(np.abs(zcoord_trunc)>=cutoff)])-1)*0.0125  # cnts* timestep

               else:
                 # import pdb; pdb.set_trace()
                  ads_cross_time.append(time[idx2])

                  if found == False:
                     zcoord_trunc = zcoord[idx2:]
                     trunc_ndx=idx2
                     print(trunc_ndx)
                     time_trunc = time[idx2:]
                     found = True
                  n=1
                  while i < len(zcoord_trunc)-1:
                        if n%2!=0:
                           if zcoord_trunc[i]<=cutoff*(-1)**n:
                               ads_cross_time.append(time_trunc[i])
                               i+=1
                               n+=1
                           else:
                               i+=1
                        elif n%2==0:
                           if zcoord_trunc[i]>=cutoff*(-1)**n:
                              ads_cross_time.append(time_trunc[i])
                              i+=1
                              n+=1
                           else:
                              i+=1




        j = 1

        ads_passage=[]
        while j < len(ads_cross_time)-2:
              ads_passage.append(ads_cross_time[j+1]-ads_cross_time[j])
              j+=1
        ads_total = sum(ads_passage)


        ave_ads = ads_total/len(ads_passage)
        ads_passage_sem = stats.sem(ads_passage)

        # abadon the first record
        time_len= time[-1]-ads_cross_time[1]
        ads_time = len(zcoord[np.where(np.abs(zcoord)>=cutoff)]-1)*0.0125 - ads_cross_time[1]
        ads_percent = ads_time /time_len

# calculate adsoprtion velocity

        def xcoord_diff(consec):
            xdiff = []

            for ndx in consec:
                xdiff.append(xcoord_trunc[ndx+1] - xcoord_trunc[ndx])
            return xdiff

        # upper wall ads index
        xcoord_trunc = xcoord[trunc_ndx:]
        ads_upper_ndx = np.where(zcoord_trunc>=cutoff)[0]

        #If they are consecutive, caculate the xcoord difference
        ads_ndx_consec_upper = ads_upper_ndx[np.where(np.diff(ads_upper_ndx)==1)]

    #    ads_xdiff_upper = np.diff(xcoord[ads_ndx_consec_upper])
        ads_xdiff_upper=xcoord_diff(ads_ndx_consec_upper)
        # save all velocity for each time

    #    import pdb; pdb.set_trace()
        ads_lower_ndx = np.where(zcoord_trunc<=-cutoff)[0]
        ads_ndx_consec_lower = ads_lower_ndx[np.where(np.diff(ads_lower_ndx)==1)]
        #ads_xdiff_lower = np.diff(xcoord_trunc[ads_ndx_consec_lower])
        ads_xdiff_lower=xcoord_diff(ads_ndx_consec_lower)

        #   time (step) used in v/t
        vel_dt = time[2]-time[1]
        # ads_shape_upper = ads_xdiff_upper.shape[0]%10
        # ads_shape_lower = ads_xdiff_lower.shape[0]%10
        #
        # ads_vel_all_upper =np.mean(ads_xdiff_upper[:-ads_shape_upper].reshape(-1,10),axis=1)/vel_dt
        # ads_vel_all_lower =np.mean(ads_xdiff_lower[:-ads_shape_lower].reshape(-1,10),axis=1)/vel_dt
        #
        # ads_vel_all = np.hstack((ads_vel_all_upper,ads_vel_all_lower))
        ads_upper_vel_all = ads_xdiff_upper/vel_dt
        ads_lower_vel_all = ads_xdiff_lower/vel_dt
        ads_each_time_upper = time_trunc[ads_ndx_consec_upper]
        ads_each_time_lower = time_trunc[ads_ndx_consec_lower]
        ads_vel_mean = (np.mean(ads_xdiff_upper)+np.mean(ads_xdiff_lower))/(2*vel_dt)  #m/s
    #    import pdb; pdb.set_trace()


        #desorption is the others
        des_ndx = np.where(np.abs(zcoord_trunc)<cutoff)[0]
        des_ndx_consec = des_ndx[np.where(np.diff(des_ndx)==1)]


        des_xdiff = xcoord_diff(des_ndx_consec)
        des_vel_all = des_xdiff/vel_dt
        des_vel_mean = np.mean(des_xdiff)/vel_dt
        des_each_time = time[des_ndx_consec]
    #    import pdb; pdb.set_trace()


    #create dataframe to save in excel

    # Format to percentage with 2 decimal
        def perc(x):
            return "{:.2%}".format(x)
        # round to 2 decimal
        def rd(x):
            return "{:.2f}".format(x)

        data = {'ads_percent':perc(ads_percent),'cross_time':rd(ave_ads),'cross_time_sem':rd(ads_passage_sem)}
        df_data = pd.DataFrame(data,index = [nuc])
        df_data['text_ct']=df_data.cross_time.astype(str) +'+-'+ df_data.cross_time_sem.astype(str)
        df = df.append(df_data)
        df.to_excel(excel_writer=sysname + "_slit_cross_time.xlsx")


## Save adsoprtion percentage, average slit crossing time (adsoprtion passage time) and standard error of the mean
        with open(vel_output_folder+outfile_prefix+ '.dat','ab') as f:
             np.savetxt(f,np.array((ads_percent,ave_ads,ads_passage_sem)),header="ads_percent" +"_" + \
             nuc + ' |ave_ads_pass_time' + ' |standard error of mean')



## Save adsoprtion and desorption velocities on average
#  use the mean velocity calculated after error progation from ads and des :ads_vel_mean_with_propagated_error.dat

        with open(vel_output_folder+'ads_des_velocity'+ '.dat','ab') as f:
             np.savetxt(f,np.array((ads_vel_mean,des_vel_mean)),header="ads_vel_mean" +"_" + \
             nuc + " |des_vel_mean" +"_" +nuc )

## Save adsoprtion and desorption velocities with time
## Later use this to calculate error progation from ads and des with confidence_interval.py

        with open(vel_output_folder+'ads_upper_all_velocity_'+nuc + '.xvg','ab') as f:
             np.savetxt(f,np.column_stack((ads_each_time_upper,ads_upper_vel_all)),header="ads_time" +"_" + \
             nuc + ' |ave_ads_vel')

        with open(vel_output_folder+'ads_lower_all_velocity_'+nuc + '.xvg','ab') as f:
             np.savetxt(f,np.column_stack((ads_each_time_lower,ads_lower_vel_all)),header="ads_time" +"_" + \
             nuc + ' |ave_ads_vel')

        with open(vel_output_folder+'des_all_velocity_'+nuc + '.xvg','ab') as f:
             np.savetxt(f,np.column_stack((des_each_time,des_vel_all)),header="ads_time" +"_" + \
             nuc + ' |ave_ads_vel')


if __name__ == "__main__":

    # Command line arguments
    parser = argparse.ArgumentParser()

    parser.add_argument('infile_prefix', help='Prefix for input file name.')
    parser.add_argument('infile_extension', help='Extension for input file name.')
    parser.add_argument('outfile_prefix', help='Prefix for output file name.')
    parser.add_argument('sysname', help='system name.')
    parser.add_argument('-nucs', '--nuc_list', nargs='+',
                        default=['dCMP', 'dTMP', 'dGMP', 'dAMP'],
                        choices=['dCMP', 'dTMP', 'dGMP', 'dAMP'],
                        help='Names of nucleotides. Defalut is dCMP dTMP dGMP dAMP.')

    args = parser.parse_args()

    infile_prefix = args.infile_prefix
    infile_extension = args.infile_extension
    outfile_prefix = args.outfile_prefix
    sysname = args.sysname
    nuc_list = args.nuc_list

    vel_output_folder = "./velocities/"

    # cutoff dw =0.6 based on PMF
    # convert to zcoord: |z| +dw = channel_sasa_width(sw)/2 +0.388
    # sasa_width box - (2*(sasa-SH)-1.5)
    # sasa_width poxy: 7.48-(2*2.18) -1.5 = 1.62
    # sasa_width cooch3: 7.07-(2*2.05) -1.5 = 1.47
    # sasa_width ch3: 7.12-(2*1.83)-1.5 = 1.96
    r_probe = 0.388

    if sysname == 'ch3':
       sw=1.96
       dw = 0.4
       cutoff = sw/2 + r_probe -dw
    elif sysname == 'cooch3':
       sw = 1.47
       dw=0.4
       cutoff = sw/2 + r_probe -dw
    else:
       sw = 1.62
       dw=0.4                        # dw cutoff from PMF and nonequil probility
       cutoff = sw/2 + r_probe -dw   # convert dw to z cutoff


    slit_cross_time(infile_prefix,infile_extension,outfile_prefix,sysname,nuc_list)
