for cosmo in {1..1} ; do 
   if [ $cosmo -lt 10 ]; then
      model_str='0'$cosmo
   else
      model_str=$cosmo
   fi
   cp -v ./Submit.sh f_R/$model_str/
   cd f_R/$model_str 
   sbatch Submit.sh 0$model_str
   #bash Submit.sh 0$model_str
   cd ../../ 
done
