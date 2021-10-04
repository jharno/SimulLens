run_flag=''       #kappa+delta+C_ell, MapsNz
#run_flag='_gamma' #shear




for cosmo in {1..1} ;  do
   let counter=$cosmo+1
   if [ $cosmo -lt 10 ]; then
      model_str='0'$cosmo
   else
      model_str=$cosmo
   fi
   echo Compiling FORGE model $cosmo $model_str
   #nsnap=$(head -n $counter snap_table.dat |tail -n 1 |gawk {'print $2'})
   let nsnap=$(wc -l checkpoints_$model_str)
   echo nsnap=$nsnap
   let nslice=$nsnap-1
   echo Nslice=$nslice

   # Extract cosmo and nslice from file
   omegam=$(head -n $counter Nodes_Omm-S8-h-fR0-sigma8-As-B0_LHCrandommaximin_Seed1_Nodes50_Dim4_AddFidTrue_extended.dat | tail -n 1| gawk {'print $1'})
   h=$(head -n $counter Nodes_Omm-S8-h-fR0-sigma8-As-B0_LHCrandommaximin_Seed1_Nodes50_Dim4_AddFidTrue_extended.dat | tail -n 1| gawk {'print $3'})
   echo para: omega = $omegam, h =  $h nslice = $nslice;
  
   rm Lens.fh ; 
   cp -v Lens_fr_template.fh Lens_fr.tmp0
   sed '3 s/OMEGAM/'"$omegam"'/' Lens_fr.tmp0  > Lens_fr.tmp1
   sed '5 s/HUBBLE/'"$h"'/' Lens_fr.tmp1  > Lens_fr.tmp2
   sed '20 s/NPLANES_MINUS_ONE/'"$nslice"'/' Lens_fr.tmp2 > Lens_fr.tmp3

   if [ $cosmo -lt 10 ]; then
      sed '49 s/NODE/0'"$cosmo"'/' Lens_fr.tmp3 > Lens_fr.tmp4
      sed '50 s/NODE/0'"$cosmo"'/' Lens_fr.tmp4 > Lens.fh
   else
      sed '49 s/NODE/'"$cosmo"'/' Lens_fr.tmp3 > Lens_fr.tmp4
      sed '50 s/NODE/'"$cosmo"'/' Lens_fr.tmp4 > Lens.fh
   fi
   #cat Lens.fh
   rm Lens_fr.tmp?

   # SimulLens & random shift files:
   #make clean; make; mkdir -p 'f_R/'$model_str; mv -v SimulLens 'f_R/'$model_str'/SimulLens'$run_flag; 
   #cp CosmoDist.py 'f_R/'$model_str'/';cp checkpoints_$model_str* 'f_R/'$model_str/
   #mkdir 'f_R/'$model_str/random_shifts/
   #for cone in {1..25}; do 
   #   tail -n $nsnap random_shifts/random_shift_subvol_LOS_cone$cone > 'f_R/'$model_str/random_shifts/random_shift_subvol_LOS_cone$cone 
   #done
   
   #MapsNz
   make clean; make MapsNz; mkdir -p 'f_R/'$model_str; mv -v MapsNz 'f_R/'$model_str'/MapsNz'$run_flag; 
   cp CosmoDist.py 'f_R/'$model_str'/'; cp CosmoDist_par.py 'f_R/'$model_str'/'; cp checkpoints_$model_str* 'f_R/'$model_str/

   # tidal fields:
   #ln -s Lens_cosmoSLICS_$cosmo.fh Lens.fh; make clean; make; mkdir 'tidal/'$model; mv -v SimulLens 'tidal/'$model'/SimulLens_tidal_nosmooth';
done
