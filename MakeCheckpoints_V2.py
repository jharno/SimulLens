import numpy as np
import deepdish as dd
from scipy.io import FortranFile
import pickle
from astropy import cosmology
from astropy.cosmology import w0waCDM
import astropy.units as u
from astropy.cosmology import z_at_value

pc_to_m=3.0854e+16
H0_100_over_c=1.5673e-26


# hdf5 data:
suffix=''
#suffix='_2048'
#suffix='_2048_hr'
#data = dd.io.load('lightcone_power'+suffix+'.h5')

# pickled data
#suffix='_4096'
#suffix='_4096_hr'
#suffix='_4096_2'
#suffix='_4096_4'
#suffix='_8192_hr'
#with open('lightcone_power'+suffix+'.pickle', 'rb') as handle:
#    data = pickle.load(handle)
#data["cosmology"]
#data["distance"]
#maps = data["rho"]

#z-projections
#proj_suffix='xy'
#proj_id=0

#y-projections
#proj_suffix='xz'
#proj_id=1

#x-projections
proj_suffix='yz'
proj_id=2

def read_maps(num_cosmo=0, base=None):
    #kappa = np.zeros((2048,2048))
    maps = []; zz = []; chi_s = []
    #for isnap in range(0,90):
    for isnap in range(1,90):
        fname = base + "cosmo_map_{0}_{1}.h5".format(num_cosmo, isnap)
        try:
            data = dd.io.load(fname)
            #_rho = data['rho'][0]/np.mean(data['rho'][0])
            # Get cosmology:
            #cosmo_data = data["cosmology"]
            #h=cosmo_data['hubble']
            #H0_in = 100.0*h
            #Omegam_in = cosmo_data['omega_matter']

            #[0,1,2] are for x,y,z projections
            _rho = data['rho'][proj_id]
            maps.append(_rho)
            zz.append(1./data['time'] - 1.0)              
            chi_s.append(data['drange'][1])              
        except:
            continue
    maps = np.array(maps)
    zz   = np.array(zz)
    chi_s  = np.array(chi_s)

    #print(Omegam_in, H0_in)
    #Omegal_in = 1.0 - Omegam_in
    #cosmo = w0waCDM(H0=H0_in, Om0=Omegam_in, Ode0=Omegal_in, w0=-1.0, wa=0.0)
    return maps, zz, chi_s

num_cosmo = 0
maps, z1, chi_s = read_maps(num_cosmo = num_cosmo, base = "../../BACCO/shared/Lgrid/test_grid/ring{0}/0.00/".format(num_cosmo)) # that's the test sims, with ring[0-36]
#maps, z1, chi_s= read_maps(num_cosmo = num_cosmo, base = "../../BACCO/shared/L512/test_grid/ring1/0.00/")                               # that's the rescaled sims, with ring[0-4]


print('z_lenses')
for i in range(0,np.shape(z1)[0]):
   print(*z1[i], i)

H0_in = 67.77
Omegam_in = 0.307112
Omegal_in = 1.0 - Omegam_in

cosmo = w0waCDM(H0=H0_in, Om0=Omegam_in, Ode0=Omegal_in, w0=-1.0, wa=0.0)
print('z_sources')
for i in range(0,np.shape(z1)[0]):
   z_s = z_at_value(cosmo.comoving_distance, chi_s[i]/0.6777*u.Mpc, zmin=0.000)
   print(z_s)
