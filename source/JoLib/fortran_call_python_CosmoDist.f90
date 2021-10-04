program SystemTest
  !real, parameter :: Om_m=0.2905, Om_l = 0.7095, w_de = -1.0
  real Om_m, Om_l, w_de, h_0, z
  character(len=7) Om_m_str, Om_l_str, w_de_str, h_0_str, z_str

  call GetArg(1,Om_m_str)
  read (Om_m_str,'(f7.3)') Om_m
  call GetArg(2,Om_l_str)
  read (Om_l_str,'(f7.3)') Om_l
  call GetArg(3,w_de_str)
  read (w_de_str,'(f7.3)') w_de
  call GetArg(4,h_0_str)
  read (h_0_str,'(f7.3)') h_0
  call GetArg(5,z_str)
  read (z_str,'(f7.3)') z
  
  print *, "FORTRAN: Calling python cosmological distance"
  print *, "with Om_m = ", Om_m, ", Om_l = ", Om_l, ", w_de = ", w_de, " and h_0 = ", h_0 
  print *, "at redshift z = ", z

  call system("rm cosmo_par.tmp")
  open(30,file='cosmo_par.tmp')
  write(30,*) Om_m, Om_l, w_de, h_0, z
  close(30)  
  call system("cat cosmo_par.tmp")

  call system("python CosmoDist.py ")
  print *, "FORTRAN: EXIT"
end program SystemTest
