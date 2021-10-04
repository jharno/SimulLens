

program CrossCorrelateMaps
  implicit none
  !include 'Lens.fh'
  

  real, parameter ::  Area = 60.85
  integer, parameter::  npc=407

  real, dimension(npc,npc) :: map1
  real, dimension(npc,npc) :: map2
  real, dimension(npc*2, npc*2) :: zeropad_map1, zeropad_map2
  real, dimension(npc) :: power

  real angle, pi
  integer i, ii!, zmin, cur_z, cur_z_s
  character(len=512)Map_str1, Map_str2, fn

  character(*), parameter :: output_form=   'binary'
  character(*), parameter :: output_access= 'sequential' 
  character(*), parameter :: input_form=    'binary' 


  call GetArg(1,Map_str1)
  call GetArg(1,Map_str2)

  pi=acos(-1.)
  angle = sqrt(Area)/180*pi

  ! Get maps
  write(*,*) 'Reading maps...'        

        open(10,file=trim(Map_str1), access = output_access, form = output_form, status = 'old')
        read(10) map1(:,:)
        close(10)
!----------------
        open(11,file=trim(Map_str2), access = output_access, form = output_form, status = 'old')
        read(11) map2(:,:)
        close(11)


  write(*,*) 'Done reading maps'

  zeropad_map1 = 0
  zeropad_map1(npc/2+1:3*npc/2, npc/2+1:3*npc/2) = map1
  zeropad_map2 = 0
  zeropad_map2(npc/2+1:3*npc/2, npc/2+1:3*npc/2) = map2
 
  !open(12,file='../PlanckCFHTMaps/testMap.dat', access = output_access, form = output_form, status = 'replace')
  !write(12) zeropad_map1(:,:)
  !close(12)
  write(*,*) 'zero padded the two maps'

  !---------------------------
  write(*,*) 'Autopower CFHT'
  power = 0

  call ps2_deconvolve(zeropad_map1,zeropad_map1,power,npc*2)
  
  fn='../PlanckCFHTMaps/l2cl_kappa_CFHT.dat'

  open(30,file=trim(fn))
  do i=1,npc!/2
     write(30,*) i*(2*pi/angle - 1),power(i)
  enddo
  close(30)

 !---------------------------
  write(*,*) 'Autopower CMB'
  power = 0

  call ps2_deconvolve(zeropad_map2,zeropad_map2,power,npc*2)
  
  fn='../PlanckCFHTMaps/l2cl_kappa_Planck.dat'

  open(30,file=trim(fn))
  do i=1,npc!/2
     write(30,*) i*(2*pi/angle - 1),power(i)
  enddo
  close(30)

 !---------------------------
  write(*,*) 'Cross-power CFHT/CMB'
  power = 0

  call ps2_deconvolve(zeropad_map1,zeropad_map2,power,npc*2)
  
  fn='../PlanckCFHTMaps/l2cl_kappa_Planck_CFHT.dat'

  open(30,file=trim(fn))
  do i=1,npc!/2
     write(30,*) i*(2*pi/angle - 1),power(i)
  enddo
  close(30)


end program CrossCorrelateMaps
