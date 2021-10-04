from astropy import cosmology
from astropy.cosmology import w0waCDM
import astropy.units as u
from astropy.cosmology import z_at_value
import numpy as np

#print("*** USING PYTHON COSMODIST ***")
cosmo_par_array = np.loadtxt('cosmo_par.tmp', usecols=(0,1,2,3,4), unpack=True)
#file = np.loadtxt("./snapshots.dat",comments='%#',usecols=(0,4,5,6,7,8),unpack=True)

#print(cosmo_par_array)

Om_m = cosmo_par_array[0]
Om_l = cosmo_par_array[1]
w_de = cosmo_par_array[2]
h_0 = cosmo_par_array[3]
z = cosmo_par_array[4]

#cosmo = w0waCDM(H0=68.98, Om0=0.2905, Ode0=0.7095, w0=-1.0, wa=0.0)
cosmo = w0waCDM(H0=h_0*100.0, Om0=Om_m, Ode0=Om_l, w0=w_de, wa=0.0)

#print("***")
#print('cosmo:')
#print(cosmo)
#print("***")

#-------------------
# Example
#Set the width of boxes to collapse:
#Box_Mpc_over_h = (505.0/2.0)*u.Mpc
#Box = Box_Mpc_over_h / cosmo.h # in Mpc

#print("Box size:")
#print(Box*cosmo.h)
#print(Box_Mpc_over_h)
#print("***")

#z_test = z_at_value(cosmo.comoving_distance, Box)
#chi_test = cosmo.comoving_distance(z_test)
#chi_test_Mpc_over_h = chi_test*cosmo.h
#print("comoving distance to the first redshift plane:")
#print(chi_test_Mpc_over_h)
#print("***")
#----------------

chi_Mpc_over_h = cosmo.comoving_distance(z)*cosmo.h
#print(chi_Mpc_over_h)
#print("*** DONE PYTHON COSMODIST ***")

file_out_str= "./chi_from_python.tmp"
data = [chi_Mpc_over_h.value, Om_m, Om_l, w_de, h_0, z]
#print(data)
np.savetxt(file_out_str,data, fmt='%.8e')

