import numpy as np
import sys
from shift_pmf import substract
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.signal import savgol_filter


# run zcoord_to_prob.py dGMP_base_coord.xvg dGMP_SH1_coord.xvg dGMP_SH2_coord.xvg dGMP cooch3
pmf_base=sys.argv[1]
pmf_SH1=sys.argv[2]
pmf_SH2=sys.argv[3]
dnmp=sys.argv[4]
sysname=sys.argv[5]



traj_data1=substract(pmf_base,pmf_SH1)
traj_data2=substract(pmf_base,pmf_SH2)
dist1=traj_data1[:,1]
dist2=traj_data2[:,1]
min_base_SH=[]
a=0
while a<len(traj_data1):
    min_base_SH.append(min(abs(dist1[a]), abs(dist2[a])))
    a+=1

if sysname == 'ch3':
   dsurf = 1.83
elif sysname == 'cooch3':
   dsurf = 2.05
else:
   dsurf = 2.18




dw = np.subtract((min_base_SH), (dsurf -0.388))
(freq, bin_edges) = np.histogram(dw,bins=np.arange(np.min(dw),np.max(dw),0.007))
#(freq2, bin_edges2) = np.histogram(dw,bins=30)

bin_width = bin_edges[1]-bin_edges[0]
bins_redge = bin_edges[1:]
prob = freq /float(np.sum(freq))
density = prob/bin_width
yy = savgol_filter(density,21,3)
zz = savgol_filter(yy,21,2)

#prob2 = freq2 /float(np.sum(freq2))
#density2 = prob2/(bin_edges2[1]-bin_edges2[0])



#bins_redge2 = bin_edges2[1:]
#xnew = np.linspace(bins_redge.min(),bins_redge.max(),200)
#f= CubicSpline(bins_redge,yy)
#power_smooth = f(xnew)

smooth_file =np.column_stack((bins_redge, yy))

np.savetxt(dnmp+'_zcoord_to_prob_smooth' + sysname+ '.dat',smooth_file )

plt.plot(bins_redge,density)
plt.plot(bins_redge,zz)


np.savetxt(dnmp+'_zcoord_to_prob'+ '.dat',np.column_stack((bins_redge,density)),header="dw|zcoord_prob")






