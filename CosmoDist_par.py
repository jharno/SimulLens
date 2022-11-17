#from astropy import cosmology
from astropy.cosmology import w0waCDM
#import astropy.units as u
from astropy.cosmology import z_at_value
from numpy import loadtxt, savetxt, size

#print("*** USING PYTHON COSMODIST ***")
#cosmo_par_array = np.loadtxt('cosmo_par.tmp', usecols=(0,1,2,3,4), unpack=True)
cosmo_par_array = loadtxt('cosmo_par_parallel.tmp', unpack=True)
#file = np.loadtxt("./snapshots.dat",comments='%#',usecols=(0,4,5,6,7,8),unpack=True)

#print(cosmo_par_array)

Om_m = cosmo_par_array[0]
Om_l = cosmo_par_array[1]
w_de = cosmo_par_array[2]
h_0 = cosmo_par_array[3]
z = cosmo_par_array[4:]

N_z = size(z)

cosmo = w0waCDM(H0=h_0*100.0, Om0=Om_m, Ode0=Om_l, w0=w_de, wa=0.0)

chi_Mpc_over_h = cosmo.comoving_distance(z)*cosmo.h
#print(chi_Mpc_over_h)
print("*** DONE PYTHON COSMODIST ***")

file_out_str= "./chi_from_python_par.tmp"
data = [chi_Mpc_over_h.value, Om_m, Om_l, w_de, h_0, z]
#print(data)
savetxt(file_out_str,chi_Mpc_over_h.value, fmt='%.8e')
