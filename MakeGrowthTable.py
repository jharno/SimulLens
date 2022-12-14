# Routine to construct an interpolation table for Growth factor
from pyccl import Cosmology, growth_factor
from numpy import loadtxt, savetxt, zeros

# 1-read in cosmology
cosmo_file = loadtxt('Cosmo4Growth.dat')

Omega_c = cosmo_file[0]
Omega_b = cosmo_file[1]
h       = cosmo_file[2]
sigma_8 = cosmo_file[3]
w0      = cosmo_file[4]
wa      = cosmo_file[5]
n_s     = cosmo_file[6]
# 2 Make the cosmology object
cosmo = Cosmology(Omega_c=Omega_c, Omega_b=Omega_b, h=h, sigma8 = sigma_8, n_s = n_s,w0 = w0,wa = wa, Omega_k=0.0)

# Construct the redshift table
nzbins = 1000
z_min = 0
z_max = 10
dz = (z_max - z_min)/nzbins
#print("z , D")
for i in range(nzbins) : 
   z = dz*i
   a = 1.0/(1.0 + z)
   D = growth_factor(cosmo,a)
   print(z,D) 

