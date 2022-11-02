# SimulLens

# Installation guide:

#1-download SimulLens from git:

#2-descend  into the main directory and edit Makefile:
>cd SimulLens
>cp Makefile_iris Makefile

#In particular, edit lines 9 & 10, currently:
#FFTWLIB=/home/jharno/data/lib/lib/
#FFTWINC=/home/jharno/data/lib/include/
#Make sure these instead point to the FFTW3 libraries in your system. These libraries must be compiled with gfortran with float precision

#3-create required soft links:
>ln -s Lens.fh  Lens_Lgrid.fh
>ln -s checkpoints_Lgrid checkpoints
>ln -s checkpoints_Lgrid_zs checkpoints_zs
>cd source; ln -s ../Lens.fh ./

#4- Compile: 
make; cp SimulLens SimulLens_Lgrid 

#5- Make sure your projection maps are in the correct format. You can do this by running this python code:
python ExtractProj.py
# This will create either the xy, xz and yz projections. See inside, there are lines to comment/uncomment to make the three cases. 

#5- Run
./SimulLens_Lgrid '' 1 10 shared/Lgrid/test_grid/ring0/0.00/ shared/Lgrid/test_grid/ring0/0.00/ './random_shifts/' '' 000

#6- When changing simulations, you to need:
#a- Generate new projection files with the correct format. See ExtractProj.py for an example.
#b- Copy the Lens.fh file, modify the cosmology, resolution nc, box size, nslices... therein
#c- Modify the cosmology in MakeCheckpoints_V2.py
#d- Generate new checkpoints and checkpoints_zs files with the same structure using MakeCheckpoints_V2.py. (Keep both columns for the checkpoints.)
#e- Update soft links pointing to the checkpoints, checkpoints_zs, Lens.fh and data
#f- Compile and run



