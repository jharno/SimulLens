

program GetGauss
  implicit none
  !include 'Lens.fh'
  

  real, parameter :: omegam = 0.2905, h = 0.6898, box = 505., Area = 60.
  integer, parameter:: nc = 12288, nslice=18, npc=6000, nk = 6000 !100000

  character(*), parameter :: Lens_output_path='/VN90/jharno/KiDS_6000_finer/'

  real, dimension(npc,npc) :: map_in_out!, ran1, ran2
  !real, dimension(npc+2,npc) :: map_in_out, ran1, ran2
  !real, dimension(npc*2, npc*2) :: zeropad_map
  real, dimension(npc) :: power, power_gauss

  real, dimension(nk) :: ell
  real, dimension(nk) :: C_ell
 
  real angle, pi, z_lens(nslice)
  integer i, cur_z !, ii, zmin, cur_z, cur_z_s
  character(len=7) z_string,LOS_str
  character(len=180) :: fn

  character(*), parameter :: output_form=   'binary'
  character(*), parameter :: output_access= 'sequential' 
  character(*), parameter :: input_form=    'binary' 

  call GetArg(1,LOS_str)
  pi=acos(-1.)
  angle = sqrt(Area)/180.*pi ! in rad
  
  cur_z = 6

  !Read lens/source distributions
  open (11,file='List_redshifts', status = 'old')
  do i = 1,nslice
     read (11, *) z_string
     read (z_string,'(f7.3)') z_lens(i)
  end do
  close(11)




  write (fn,'(f5.3,"kappa_weight.dat")') z_lens(cur_z)
  open(10,file=Lens_output_path//trim(fn)//'_LOS'//trim(LOS_str), access = output_access, form = output_form, status = 'old')
  read(10) map_in_out(1:npc,1:npc)
  close(10)
  write(*,*) 'Read ', trim(fn)  
  !map_in_out = map_in_out*64

  !ran1 = map_in_out

  power = 0.0

  !call ps2_deconvolve(map_in_out,map_in_out,power,npc)
  !write(*,*) 'Got Power spectrum'
  
  !power = power*2048.

  do i = 1,npc
    ell(i) = i*(2*pi/angle -1)
  !  C_ell(i)  = power(i)*2*pi/ell(i)**2 
   !write(*,*) ell(i), C_ell(i)
  enddo


  !map_in_out = ran1

  !open(11,file='./elCl_par_WMAP9_SN_BAO_z0_z0.582_nokcut.dat', status = 'old')  
  !do i=1,nk
  !   read(11,*) ell_camb(i), C_ell_camb(i)
  !enddo
  !close(11)

  !do i = 1,nk
  !  write(*,*) ell_camb(i), C_ell_camb(i)
  !enddo
  
  !stop

  !call GaussRandomField_2d_r2c(map_in_out, npc, ell, C_ell, npc,ran1,ran2 )
  

  !--------
  !write(*,*) 'Got Gaussian Field, computing power spectrum...'
  !power_gauss = 0

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

  call ps2_deconvolve(map_in_out,map_in_out,power, npc)
  !call ps2_camb(map_in_out,map_in_out,power,power_gauss, npc)
 

  !------------

 
  !fn='l2cl_kappa_cic.dat'
  fn='l2cl_kappa_ngp.dat'
  !fn='l2cl_kappa_camb.dat'

  read (z_lens(6), *) z_string

  open(30,file=trim(z_string)//Lens_output_path//trim(fn))
  do i=1,npc/2  ! Remove the '/2' for zeropad maps
     write(30,*) ell(i),power(i), power_gauss(i)
  enddo
  close(30)
  write(*,*) 'Got l2Cl'


end program GetGauss
