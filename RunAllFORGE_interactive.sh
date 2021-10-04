for cosmo in {0..0} ;  do
   let counter=$cosmo+1
   seed='_a'
   model=$cosmo$seed
   echo Compiling FORGE model $cosmo $model
   nsnap=$(head -n $counter snap_table.dat |tail -n 1 |gawk {'print $2'})
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

   make clean; make; mkdir -p 'fr/'$model; mv -v SimulLens 'fr/'$model'/'; cp CosmoDist.py 'fr/'$model'/'

   #ln -s Lens_cosmoSLICS_$cosmo.fh Lens.fh; make clean; make; mkdir 'tidal/'$model; mv -v SimulLens 'tidal/'$model'/SimulLens_tidal_nosmooth';
done
