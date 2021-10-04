subroutine pcor_1d(a,b,ab,n)
  implicit none
  include 'fftw3.f'

  integer n
  complex, dimension(0:n/2) :: ca,cb,cab
  real, dimension(0:n-1) :: a,b,ab
  integer(8) :: plan,iplan

  call sfftw_plan_dft_r2c_1d(plan,n,a,ca,FFTW_ESTIMATE)
  call sfftw_execute(plan)
  call sfftw_destroy_plan(plan)
  call sfftw_plan_dft_r2c_1d(plan,n,b,cb,FFTW_ESTIMATE)
  call sfftw_execute(plan)
  call sfftw_destroy_plan(plan)
  ca=ca*conjg(cb)/real(n)**1
  call sfftw_plan_dft_c2r_1d(iplan,n,ca,ab,FFTW_ESTIMATE)
  call sfftw_execute(iplan)
  ab = ab / real(n)**1
  call sfftw_destroy_plan(iplan)

  return
end subroutine pcor_1d

subroutine npcor_1d(a,b,ab,n)
  implicit none

  integer n,n2
  real,dimension(0:2*n-1) :: aa,bb,cc
  real,dimension(0:n-1) :: a,b,ab

  aa=0
  bb=0
  aa(0:n-1)=a
  bb(0:n-1)=b
  call pcor_1d(aa,bb,cc,2*n)
  ab=cc(0:n-1)
  
  return
end subroutine npcor_1d

subroutine pcor_2d(a,b,ab,n)
  implicit none
  include 'fftw3.f'

  integer n,sfft_2d
  complex, dimension(0:n/2,0:n-1) :: ca,cb,cab
  real, dimension(0:n-1,0:n-1) :: a,b,ab
  integer(8) :: plan,iplan

  call sfftw_plan_dft_r2c_2d(plan,n,n,a,ca,FFTW_ESTIMATE)
  call sfftw_execute(plan)
  call sfftw_destroy_plan(plan)
  call sfftw_plan_dft_r2c_2d(plan,n,n,b,cb,FFTW_ESTIMATE)
  call sfftw_execute(plan)
  call sfftw_destroy_plan(plan)
  cab=ca*conjg(cb)/real(n)**2
  call sfftw_plan_dft_c2r_2d(iplan,n,n,cab,ab,FFTW_ESTIMATE)
  call sfftw_execute(iplan)
  ab = ab / real(n)**2
  call sfftw_destroy_plan(iplan)

  return
end subroutine pcor_2d

subroutine npcor_2d(a,b,ab,n)
  implicit none

  integer n,n2
  real,dimension(0:2*n-1,0:2*n-1) :: aa,bb,cc
  real,dimension(0:n-1,0:n-1) :: a,b,ab

  aa=0
  bb=0
  aa(0:n-1,0:n-1)=a
  bb(0:n-1,0:n-1)=b
  call pcor_2d(aa,bb,cc,2*n)
  ab=cc(0:n-1,0:n-1)
  
  return
end subroutine npcor_2d

subroutine pcor_3d(a,b,ab,n)
  implicit none
  include 'fftw3.f'

  integer n
  complex, dimension(0:n/2,0:n-1,0:n-1) :: ca,cb,cab
  real, dimension(0:n-1,0:n-1,0:n-1) :: a,b,ab
  integer(8) :: plan,iplan

  call sfftw_plan_dft_r2c_3d(plan,n,n,n,a,ca,FFTW_ESTIMATE)
  call sfftw_execute(plan)
  call sfftw_destroy_plan(plan)
  call sfftw_plan_dft_r2c_3d(plan,n,n,n,b,cb,FFTW_ESTIMATE)
  call sfftw_execute(plan)
  call sfftw_destroy_plan(plan)
  ca=ca*conjg(cb)/real(n)**3
  call sfftw_plan_dft_c2r_3d(iplan,n,n,n,ca,ab,FFTW_ESTIMATE)
  call sfftw_execute(iplan)
  ab = ab / real(n)**3
  call sfftw_destroy_plan(iplan)

  return
end subroutine pcor_3d

subroutine npcor_3d(a,b,ab,n)
  implicit none

  integer n,n2
  real, dimension(0:2*n-1,0:2*n-1,0:2*n-1) :: aa,bb,cc
  real,dimension(0:n-1,0:n-1,0:n-1) :: a,b,ab

  aa=0
  bb=0
  aa(0:n-1,0:n-1,0:n-1)=a
  bb(0:n-1,0:n-1,0:n-1)=b
  call pcor_3d(aa,bb,cc,2*n)
  ab=cc(0:n-1,0:n-1,0:n-1)

  return
end subroutine npcor_3d

