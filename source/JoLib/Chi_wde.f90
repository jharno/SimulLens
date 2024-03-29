! Subroutine that calls python cosmological distances, and outputs the comoving distance chi, in Mpc/h
subroutine chi_python(z,omegam,omegal,w_de,h, chi_wde)
!
  implicit none
  real(4) z,omegam,omegal, w_de, h,  chi_wde

  !print *, 'Inside Chi_wde'

!#==============================================================

!#And then copy the data from the RAID to the local worker
!#Use lockfiles so we don't have an I/O jam when lots of workers
!#are working simultaneously



!#==========
!#Start checking for the I/O lockfile before transferring the data
!#This command staggers the time when different jobs check for the existence
!#of the all important lockfile

!#ran=`echo  $[RANDOM%300]| awk '{print $1/100}'`
!#sleep $ran
!# Look for a lockfile and wait until one doesn't exist
!#WORKERID=`hostname`
!#LOCKFILEDIR=/home/$USER/IO_lockfile
!#locked=1
!#while [ $locked = 1 ]
!#do
!#    sleep 1
!#    locked=0
!#    if mkdir $LOCKFILEDIR ; then
!#       echo "Locking succeeded" > $LOCKFILEDIR/locked_$WORKERID_$FIELD
!#    else
!#       locked=1
!#    fi
!#done
!# Unlocked, now time to transfer the data
!#============


  call system("echo Computing cosmology distances")
  call system("while [ -f ./cosmo_par.tmp ]; do echo cosmo file exists, will sleep ; sleep 1 ; done")

  !call system("rm cosmo_par.tmp")
  !call system("rm chi_from_python.tmp")
  open(30,file='cosmo_par.tmp')
  write(30,*) omegam, omegal, w_de, h, z
  close(30)
  !call system("cat cosmo_par.tmp")
  call system("python CosmoDist.py")
  !call system("cat chi_from_python.tmp")
  open(30,file='chi_from_python.tmp')
  read(30,*) chi_wde
  close(30)

  call system("rm cosmo_par.tmp")
  call system("rm chi_from_python.tmp")
  call system("echo Done!")

!# Once copying is completed delete the lockfile so another
!# process can begin

!#rm -rf $LOCKFILEDIR
!#echo "removed lockfile"


  !chi_wde = -999.0
  !chi_wde = chi
  return 
end subroutine chi_python
!end function chi_wde
