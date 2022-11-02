import numpy as np
import deepdish as dd
from scipy.io import FortranFile
import pickle

# hdf5 data:
suffix=''

#############
# Must run this code for the x, y, and z projections, by commenting out the appropriate sections below:

#z-projections
proj_suffix='xy'
proj_id=0

#y-projections
#proj_suffix='xz'
#proj_id=1

#x-projections
#proj_suffix='yz'
#proj_id=2

#############


def read_maps(num_cosmo=0, base=None):
    #kappa = np.zeros((2048,2048))
    maps = []; zz = []
    for isnap in range(0,90):
        fname = base + "cosmo_map_{0}_{1}.h5".format(num_cosmo, isnap)
        try:
            data = dd.io.load(fname)
            #_rho = data['rho'][0]/np.mean(data['rho'][0])

            #[0,1,2] are for x,y,z projections
            _rho = data['rho'][proj_id]
            maps.append(_rho)
            zz.append(data['time'])
        except:
            continue
    maps = np.array(maps)
    zz   = np.array(zz)
    return maps, zz

num_cosmo = 0
#maps, z1 = read_maps(num_cosmo = num_cosmo, base = "shared/Lgrid/test_grid/ring{0}/0.00/".format(num_cosmo))         # that's the test sims, with ring[0-36]
maps, z1 = read_maps(num_cosmo = num_cosmo, base = "shared/L512/test_grid/ring1/0.00/")                               # that's the rescaled sims, with ring[0-4]



print(maps[7])
print(z1)
print(np.shape(z1)[0])

node=000

# 4096_1 
#for snap in range(0,42):
     #we start with an offset to remove maps with z too large
     #grid = maps[snap+7]
#for snap in range(0,20):
     # 4096_2 
     #grid = maps[snap+1]
 # 4096_4 
for snap in range(0,np.shape(z1)[0]):
     grid = maps[snap]
     #print(grid[0:10,0:10])
     #print(grid)
     #filename  = './projection_'+proj_suffix+'_snap_' + str(snap).zfill(3) + '_Node_' + str(node).zfill(3) + suffix + '.unf'
     #filename  = './shared/Lgrid/test_grid/ring'+str(num_cosmo)+'/0.00/projection_'+proj_suffix+'_snap_' + str(snap).zfill(3) + '_Node_' + str(node).zfill(3) + suffix + '.unf'
     filename  = './shared/L512/test_grid/ring1/0.00/projection_'+proj_suffix+'_snap_' + str(snap).zfill(3) + '_Node_' + str(node).zfill(3) + suffix + '.unf'
     f = FortranFile(filename, 'w')
     f.write_record(np.array(grid, dtype=float))
     f.close()
     print('written and closed'+filename)



