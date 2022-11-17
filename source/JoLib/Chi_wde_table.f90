! Subroutine that calls python cosmological distances, and outputs the
! comoving distance chi, in Mpc/h
subroutine chi_table_python(z,omegam,omegal,w_de,h, chi_wde,Np)
!
  implicit none
  integer Np
  real(4) z(Np),omegam,omegal, w_de, h,  chi_wde(Np)


  call system("echo IN THE CHI_TABLE subroutine")
  call system("while [ -f ./cosmo_par_parallel.tmp ]; do echo cosmo file exists, will sleep ; sleep 1 ; done")

  !call system("rm cosmo_par.tmp")
  !call system("rm chi_from_python.tmp")
  open(30,file='cosmo_par_parallel.tmp')
  !write(30,*) omegam, omegal, w_de, h, z(:)
  write(30,'(f8.5)', advance='no') omegam, omegal, w_de, h, z(:)
  close(30)
  !call system("cat cosmo_par.tmp")
  call system("python CosmoDist_par.py")
  !call system("python3.7 CosmoDist_par.py")
  !call system("cat chi_from_python.tmp")
  open(30,file='chi_from_python_par.tmp')
  read(30,*) chi_wde
  close(30)

  call system("rm cosmo_par_parallel.tmp")
  call system("rm chi_from_python_par.tmp")

  return
end subroutine chi_table_python
!end function chi_wde
