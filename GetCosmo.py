from astropy import cosmology
from astropy.cosmology import w0waCDM
import astropy.units as u
from astropy.cosmology import z_at_value
import numpy as np
import sys

#f(r) from line
H0_str = sys.argv[1]
H0_in = np.float(H0_str)

Omegam_str = sys.argv[2]
Omegam_in = np.float(Omegam_str)
Omegal_in= 1.0-Omegam_in

#model = sys.argv[3]

cosmo = w0waCDM(H0=H0_in, Om0=Omegam_in, Ode0=Omegal_in, w0=-1.0, wa=0.0)

#f(r) node00
#cosmo = w0waCDM(H0=67.37, Om0=0.31315, Ode0=0.68685, w0=-1.0, wa=0.0)

#f(r) node049
#cosmo = w0waCDM(H0=72.861, Om0=0.36187, Ode0=0.63813, w0=-1.0, wa=0.0)

#FID
#cosmo = w0waCDM(H0=68.98, Om0=0.2905, Ode0=0.7095, w0=-1.0, wa=0.0)

#00
#cosmo = w0waCDM(H0=67.66, Om0=0.3282, Ode0=0.6718, w0=-1.2376, wa=0.0)

#01
#cosmo = w0waCDM(H0=71.04, Om0=0.1019, Ode0=0.8981, w0=-1.6154, wa=0.0)

#02
#cosmo = w0waCDM(H0=62.38, Om0=0.2536, Ode0=0.7464, w0=-1.7698, wa=0.0)

#03
#cosmo = w0waCDM(H0=65.84, Om0=0.1734, Ode0=0.8266, w0=-0.5223, wa=0.0)

#04
#cosmo = w0waCDM(H0=60.34, Om0=0.3759, Ode0=0.6241, w0=-0.9741, wa=0.0)

#05
#cosmo = w0waCDM(H0=74.59, Om0=0.4758, Ode0=0.5242, w0=-1.3046, wa=0.0)

#06
#cosmo = w0waCDM(H0=80.31, Om0=0.1458, Ode0=0.8542, w0=-1.4498, wa=0.0)

#07
#cosmo = w0waCDM(H0=69.40, Om0=0.3099, Ode0=0.6901, w0=-1.8784, wa=0.0)

#08
#cosmo = w0waCDM(H0=63.74, Om0=0.4815, Ode0=0.5185, w0=-0.7737, wa=0.0)

#09
#cosmo = w0waCDM(H0=80.06, Om0=0.3425, Ode0=0.6575, w0=-1.5010, wa=0.0)

#10
#cosmo = w0waCDM(H0=76.45, Om0=0.5482, Ode0=0.4518, w0=-1.9127, wa=0.0)

#11
#cosmo = w0waCDM(H0=65.05, Om0=0.2898, Ode0=0.7102, w0=-0.6649, wa=0.0)

#12
#cosmo = w0waCDM(H0=68.19, Om0=0.4247, Ode0=0.5753, w0=-1.1986, wa=0.0)

#13
#cosmo = w0waCDM(H0=78.33, Om0=0.3979, Ode0=0.6021, w0=-1.1088, wa=0.0)

#14
#cosmo = w0waCDM(H0=78.90, Om0=0.1691, Ode0=0.8309, w0=-1.6903, wa=0.0)

#15
#cosmo = w0waCDM(H0=75.67, Om0=0.1255, Ode0=0.8745, w0=-0.9878, wa=0.0)

#16
#cosmo = w0waCDM(H0=66.91, Om0=0.5148, Ode0=0.4852, w0=-1.3812, wa=0.0)

#17
#cosmo = w0waCDM(H0=62.85, Om0=0.1928, Ode0=0.8072, w0=-0.8564, wa=0.0)

#18
#cosmo = w0waCDM(H0=71.51, Om0=0.2784, Ode0=0.7216, w0=-1.0673, wa=0.0)

#19
#cosmo = w0waCDM(H0=73.88, Om0=0.2106, Ode0=0.7894, w0=-0.5667, wa=0.0)

#20
#cosmo = w0waCDM(H0=61.61, Om0=0.4430, Ode0=0.5570, w0=-1.7037, wa=0.0)

#21
#cosmo = w0waCDM(H0=81.29, Om0=0.4062, Ode0=0.5938, w0=-1.9866, wa=0.0)

#22
#cosmo = w0waCDM(H0=77.06, Om0=0.2294, Ode0=0.7706, w0=-0.8602, wa=0.0)

#23
#cosmo = w0waCDM(H0=69.88, Om0=0.5095, Ode0=0.4905, w0=-0.7164, wa=0.0)

#24
#cosmo = w0waCDM(H0=72.71, Om0=0.3652, Ode0=0.6348, w0=-1.5414, wa=0.0)

#Set the width of boxes to collapse:
Box_Mpc_over_h = (500.0/2.0)*u.Mpc
#Box_Mpc_over_h = (505.0/2.0)*u.Mpc
Box = Box_Mpc_over_h / cosmo.h # in Mpc

#print("***")
#print('cosmo:')
#print(cosmo)
#print("***")
#print("Projection size:")
#print(Box*cosmo.h)
#print(Box_Mpc_over_h)
#print("***")
z_test = z_at_value(cosmo.comoving_distance, Box)
chi_test = cosmo.comoving_distance(z_test)
chi_test_Mpc_over_h = chi_test*cosmo.h
#print("comoving distance to the first redshift plane:") 
#print(chi_test_Mpc_over_h)
#print("***")


#maxNumLenses = 33 # set to 29 for model 01
maxNumLenses = 25 # set to 29 for model 01
z_limit = 3.0
chi_0 = 0.0*u.Mpc*cosmo.h
chi = np.zeros(maxNumLenses)
z_slices = np.zeros(maxNumLenses)


#print("***")
#print('z_lens:')
#for n in range(0,maxNumLenses):
for n in range(maxNumLenses,-1,-1):
     chi_slice = ((chi_0 + Box*(n+0.5)))
     chi_slice_Mpc_over_h = chi_slice*cosmo.h
     z_slice = z_at_value(cosmo.comoving_distance, chi_slice, zmin=0.01)
     #print(np.str("{:5.4f}".format(z_slice)))
     #print(chi_slice_Mpc_over_h.value, z_slice)
     #if z_slice > z_limit:
     #    break

#print("***")
#print('z_sources:')
#for n in range(0,maxNumLenses):
for n in range(maxNumLenses,-1,-1):
     chi_slice = ((chi_0 + Box*(n+1)))
     chi_slice_Mpc_over_h = chi_slice*cosmo.h
     z_slice = z_at_value(cosmo.comoving_distance, chi_slice, zmin=0.01)
     #print(np.str("{:5.4f}".format(z_slice)))
     #print(chi_slice_Mpc_over_h.value,z_slice)
     #if z_slice > z_limit:
     #    break

#print("***")
#print('z_lens      z_sources:')
#print("***")
#for n in range(0,maxNumLenses):
for n in range(maxNumLenses,-1,-1):
     chi_slice = ((chi_0 + Box*(n+0.5)))
     chi_slice_Mpc_over_h = chi_slice*cosmo.h
     z_slice = z_at_value(cosmo.comoving_distance, chi_slice, zmin=0.01)
     chi_slice = ((chi_0 + Box*(n+1)))
     chi_slice_Mpc_over_h = chi_slice*cosmo.h
     z_slice_s = z_at_value(cosmo.comoving_distance, chi_slice, zmin=0.01)
     if z_slice <  z_limit:
        print(np.str("{:5.4f}".format(z_slice))+"      "+np.str("{:5.4f}".format(z_slice_s)))


#%---------------------------------
#z_slices = z_slices(1:sum(n_max));
#chi = chi(1:sum(n_max));

#chi_s = zeros(1,sum(n_max));
#z_s = zeros(1,sum(n_max));
#for n = 1:n_max(1)
#   chi_s(n) = chi(n) + lbox(1)/2;
#   z_s(n) = InverseChi(chi_s(n), Om);
#end

#z_output = [z_slices' z_s' ];
#chi_output = [chi' chi_s'];

