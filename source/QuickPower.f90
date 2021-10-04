

program GetPower2D
  implicit none
  !include 'Lens.fh'
  

  real, parameter ::   Area = 100.
  integer, parameter:: nc = 12288,  npc=7745


  !real, dimension(npc,npc,nslice) :: map_3D
  real, dimension(npc,npc) :: map_in_out
  !real, dimension(npc*2, npc*2) :: zeropad_map
  real, dimension(npc) :: power

  real angle, pi
  integer i
  character(len=500) :: fn_in, fn_out
  !real, dimension(70) :: z_file, Nz, dchi

  character(*), parameter :: output_form=   'unformatted'
  character(*), parameter :: output_access= 'stream'
  character(*), parameter :: input_form=    'unformatted'
  character(*), parameter :: input_access=  'stream'

  !equivalence (map_3D, zeropad_map)

  call GetArg(1,fn_in)
  call GetArg(2,fn_out)

  pi=acos(-1.)
  angle = sqrt(Area)/180*pi 
   
  ! Get maps
  write(*,*) 'Reading maps...'        
  open(14,file=trim(fn_in), access = output_access, form = output_form, status = 'old')
  read(14) map_in_out(:,:)
  close(14)
  write(*,*) 'Done reading maps'

  power = 0.0

  !---------
  ! Zero-Pad
  !---------

  !zeropad_map = 0
  !zeropad_map(npc/2+1:3*npc/2, npc/2+1:3*npc/2) = map_in_out
  !call ps2_deconvolve(zeropad_map,zeropad_map,power,npc*2)
  
  !fn='l2cl_kappa_nz.dat_zeropad_LOS'

  !-----------
  !No zero-pad
  !-----------

  call ps2_deconvolve(map_in_out,map_in_out,power,npc)
  
  open(30,file=trim(fn_out))
  do i=1,npc/2  ! Remove the '/2' for zeropad maps
     write(30,*) real(i)*(2*pi/angle - 1),power(i)
  enddo
  close(30)
  write(*,*) 'Got l2Cl'


end program GetPower2D
