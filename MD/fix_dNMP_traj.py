import numpy as np
import matplotlib.pyplot as plt
import sys
from read_xvg import read_xvg


grofile = sys.argv[1]
coordfile = sys.argv[2]
outfile_prefix = sys.argv[3]

#read box sizes from gro file to get box size & calculate box center
fid = open(grofile, 'r')
gro_data = fid.readlines()
box_sizes = np.array(gro_data[-1].split()).astype(float)
boxx = box_sizes[0]
boxxo2 = boxx/2.0
boxy = box_sizes[1]
boxyo2 = boxy/2.0
boxz = box_sizes[2]
boxzo2 = boxz/2.0
box_center = boxzo2
fid.close

#read file with time, coordinates & convert time to ns
data = read_xvg(coordfile)
data[:, 0] = data[:, 0]/1000.0

#shift coordinates above box center down by the box size
ind = data[:, 3] > box_center
data[ind, 3] = data[ind, 3] - boxz

#undo periodic boundary conditions in the x and y directions
npts = data.shape[0]
for itime in range(npts-1):
	data[itime+1, 1] = data[itime+1, 1] - boxx*round((data[itime+1, 1] - data[itime, 1])/boxx)
	data[itime+1, 2] = data[itime+1, 2] - boxy*round((data[itime+1, 2] - data[itime, 2])/boxy)

#plot fixed x,y, and z trajectories and save
# coord_names = ['x', 'y', 'z']
# for icoord in range(3):
#     plt.figure(icoord)
#     plt.plot(data[:, 0], data[:, icoord+1])
#     plt.xlabel(r'Time (ns)')
#     plt.ylabel(r'dCMP COM ' + coord_names[icoord] + ' position (nm)')
#     plt.savefig(outfile_prefix + '_' + coord_names[icoord] + '.png')

#save fixed data to outfile
np.savetxt(outfile_prefix + '.dat', data)
