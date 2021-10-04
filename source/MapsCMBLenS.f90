

program GetNzMaps
  implicit none
  !include 'Lens.fh'
  

  !real, parameter :: omegam = 0.2905, h = 0.6898, box = 505., Area = 60.
  real, parameter :: omegam = 0.2905, h = 0.6898, box = 505., Area = 100.


  !---------
  ! SLICS-HR
  !integer, parameter:: nc = 12288, nslice=18, npc=6000
  !character(*), parameter :: Lens_input_path= '/VN90/jharno/KiDS_6000_finer_highres/'
  !character(*), parameter :: Lens_output_path='/VN90/jharno/KiDS_6000_finer_highres/Maps_Nz_CFHT_Mass_Map_paper/'
  !character(*), parameter :: Lens_output_path='/VN90/jharno/KiDS_6000_finer_highres/Maps_Nz_KiDS/'
  !character(*), parameter :: Lens_output_path='/VN90/jharno/KiDS_6000_finer_highres/Maps_Nz_baryon_paper/'
  
  !---------
  ! SLICS-LE
  !integer, parameter:: nc = 12288, nslice=18, npc=6000
  !character(*), parameter :: Lens_input_path= '/VN90/jharno/KiDS_6000_finer/'
  !character(*), parameter :: Lens_output_path='/VN90/jharno/KiDS_6000_finer/Maps_Nz_CFHT_Mass_Map_paper/'
  !character(*), parameter :: Lens_output_path='/VN90/jharno/KiDS_6000_finer/Maps_Nz_KiDS/'
  !character(*), parameter :: Lens_output_path='/VN90/jharno/KiDS_6000_finer/Maps_Nz_baryon_paper/'
  !character(*), parameter :: Lens_output_path='/VN90/jharno/KiDS_6000_finer/Maps_Nz_baryon_paper_V2/'
  !---------

  !---------
  ! SLICS-LR
  !integer, parameter:: nc = 3072, nslice=18, npc=6000
  !character(*), parameter :: Lens_input_path= '/VN70-2/jharno/KiDS_6000_LR/'
  !character(*), parameter :: Lens_output_path='/VN70-2/jharno/KiDS_6000_LR/Maps_Nz_CFHT_Mass_Map_paper/'
  !character(*), parameter :: Lens_output_path='/VN70-2/jharno/KiDS_6000_LR/Maps_Nz_KiDS/'
  !character(*), parameter :: Lens_output_path='/VN70-2/jharno/KiDS_6000_LR/Maps_Nz_baryon_paper/'
  !----------
 
  !---------
  ! SLICS-LE 10deg
  !integer, parameter:: nc = 12288, nslice=18, npc=7745 ! For KiDSLenS
  integer, parameter:: nc = 12288, nslice=28, npc=7745  ! For CMBLenS
  character(*), parameter :: Lens_input_path2= './CMBLenS/'
  character(*), parameter :: Lens_input_path= './KiDSLenS/'
  character(*), parameter :: Lens_output_path= './CMBLenS/'
  !character(*), parameter :: Lens_output_path='/VN90/jharno/KiDS_6000_finer/Maps_Nz_CFHT_Mass_Map_paper/'
  !character(*), parameter :: Lens_output_path='/VN90/jharno/KiDS_6000_finer/Maps_Nz_KiDS/'
  !character(*), parameter :: Lens_output_path='/VN90/jharno/KiDS_6000_finer/Maps_Nz_baryon_paper/'
  !character(*), parameter :: Lens_output_path='/VN90/jharno/KiDS_6000_finer/Maps_Nz_baryon_paper_V2/'
  

  !character(*), parameter :: Lens_output_path='/VN90/jharno/KiDS_6000_finer/'

  real, dimension(npc,npc,nslice) :: map_3D
  real, dimension(npc,npc) :: map_in_out
  real, dimension(npc*2, npc*2) :: zeropad_map
  real, dimension(npc) :: power

  real kernel, z_write, z_source(nslice), z_lens(nslice), chi, Nz_CFHT, dz, gw, angle, pi
  integer i, ii, zmin, cur_z, cur_z_s
  character(len=7) z_string,LOS_str, LOS_str_CMB
  character(len=180) :: fn
  real, dimension(70) :: z_file, Nz, dchi

#ifdef gfort  
  character(*), parameter :: output_form=   'unformatted'
  character(*), parameter :: output_access= 'stream' 
  character(*), parameter :: input_form=    'unformatted' 
#else
  character(*), parameter :: output_form=   'binary'
  character(*), parameter :: output_access= 'sequential' 
  character(*), parameter :: input_form=    'binary' 
#endif

  !equivalence (map_3D, zeropad_map)

  call GetArg(1,LOS_str)
  call GetArg(2,LOS_str_CMB)

  pi=acos(-1.)
  angle = sqrt(Area)/180*pi  ! For UBC Lens


  !Read N(z) from C. Heymans N(z) KiDS & CFHT (Kilbinger?):
  !open (11,file='NZ_KIDS_v1.1.dat', status = 'old')
  !do i = 1,70
  !   read (11, *) z_file(i), Nz(i), Nz_CFHT
     !write(*,*) z_file(i), Nz(i)
  !end do
  !close(11)

  ! Select at most one of the following to ***overwrite***:
  !Use instead KiDS Fit function from Joachimi?
  !open (11,file='Nz_Fit2KiDS.dat', status = 'old')

  !Or the CFHT mass map  N(z):
  !open (11,file='Nz_CFHT_Mass.dat', status = 'old')

  !or the baryon paper N(z):
  !open (11,file='Nz_Fit2_CFHT_baryons.dat', status = 'old')
  !open (11,file='Nz_CFHT_Baryons.dat', status = 'old')
  !open (11,file='nofz_ZB_0.4_1.3.dat', status = 'old')
  !do i = 1,70
  !   read (11, *) z_file(i), Nz(i) ! Baryons and CFHTLenS_mass
  !   !read (11, *)  Nz(i) ! Fit2KiDS
  !   !write(*,*) z_file(i), Nz(i)
  !end do
  !close(11)

  !dz = 0.05
  !dchi = 0;

  !z_file = z_file - dz/2.0 ! z_file is upper edge of the bin, not centre
  !write(*,*)'z_file = ', z_file
  !write(*,*) 'sum of N(z) = ', sum(Nz*dz)
  !stop

  ! first bin goes from zero to between first and second point:
  !dchi(1) = 0.5*(chi(z_file(2), omegam, h)*h + chi(z_file(1), omegam,h)*h)
  !do i = 2,69
  !  !other bins go from mid points between each points:
  !  dchi(i) = 0.5*(chi(z_file(i+1), omegam, h)*h - chi(z_file(i-1), omegam, h)*h)
  !enddo
  ! last bin is given the distance between the last two points...
  !dchi(70) = (chi(z_file(70), omegam, h)*h - chi(z_file(69), omegam,h)*h)

  !write(*,*) dchi

  ! transform with Jacobian
  !Nz = Nz*dz/dchi

  ! Normalize
  !Nz = Nz/sum(Nz*dchi)

  !write(*,*) '**********'
  !write(*,*) 'Nz stats: '
  !write(*,*) '**********'
  !write(*,*) sum(Nz*dchi)
  !write(*,*) Nz, dchi
  


  !Read lens/source distributions
  open (11,file='List_redshifts_CMBLenS', status = 'old')
  do i = 1,nslice
     read (11, *) z_string
     read (z_string,'(f7.3)') z_lens(i)
  end do
  close(11)
  open (11,file='List_redshifts_source_CMBLenS', status = 'old')
  do i = 1,nslice
     read (11, *) z_string
     read (z_string,'(f7.3)') z_source(i)
  end do
  close(11)
  write(*,*)z_lens 
  write(*,*)z_source 

  ! Get maps
  write(*,*) 'Reading maps...'        

  do i = 1,nslice

        !Get Redshift
        z_write = z_lens(i)
        !write(*,*) 'Reading z=', z_write,'and LOS ', LOS_str
 
        if(i.le.22) then
           write (fn,'(f5.3,"delta.dat_bicubic")') z_write
        elseif(i.le.27)then
           write (fn,'(f6.3,"delta.dat_bicubic")') z_write
        else
           write (fn,'(f7.3,"delta.dat_bicubic")') z_write
        endif

        if(i.le.18) then
           open(10,file=Lens_input_path//trim(fn)//'_LOS'//trim(LOS_str), access = output_access, form = output_form, status = 'old')
        else
           open(10,file=Lens_input_path2//trim(fn)//'_LOS'//trim(LOS_str_CMB), access = output_access, form = output_form, status = 'old')
        endif
        read(10) map_in_out(:,:)
        close(10)
        write(*,*)'Read delta Map for z=',z_write


        map_3D(:,:,i)= map_in_out(:,:)

        !write(*,*) '**********'
        !write(*,*) 'Map stats: '
        !write(*,*) '**********'
        !write(*,*) sum(map_3D(:,:,i))
        !write(*,*) minval(map_3D(:,:,i)), maxval(map_3D(:,:,i))
  enddo

  write(*,*) 'Done reading delta maps'


  map_in_out = 0

  !----------------
  ! Loop over lenses
  do i = 1,nslice

     cur_z_s = 1

     !find minimum for the Nz integration
     zmin = 0
     do ii = 1,70
        if(z_file(ii) > z_lens(i))then
           zmin = ii
           exit
        endif
     enddo
     

     ! Get integal over sources
     !gw = 0
     !do ii=zmin, 70
     !   gw = gw + Nz(ii)*(1.0 - chi(z_lens(i),omegam,h)/chi(z_file(ii), omegam, h))*dchi(ii)
     !enddo
     gw = (1-chi(z_lens(i),omegam,h)/chi(z_source(28),omegam,h))

     !write(*,*) z_lens(i), z_file(zmin), Nz(zmin), gw
    
     ! Get Kernel
     kernel = 3./2. * omegam * box/real(nc)/(3e3**2)*chi(z_lens(i), omegam,h)*h * (1+z_lens(i)) * gw
     !write(*,*) 'Kernel4 = ', kernel


     ! Integrate
     map_in_out(:,:) = map_in_out(:,:) + map_3D(:,:,i)*kernel

  enddo 

  !--------------
  ! write to file


  fn = 'kappa_CMBLenS.dat'
  !fn = 'kappa_nz.dat'
  !fn = 'kappa_nz_mass.dat'




  open(11,file=Lens_output_path//trim(fn)//'_LOS'//trim(LOS_str), access = output_access, form = output_form)
  write(11) map_in_out(:,:)
  write(*,*)'Wrote Map for LOS', trim(LOS_str)
  close(11)
  !write(*,*) 'Commented out the writing of maps to file, doing simply l2Cl...'
  

#ifdef power_spectrum

  power = 0

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

  write(*,*)'Computing l2Cl ...'
  call ps2_deconvolve(map_in_out,map_in_out,power,npc)
  write(*,*)'Done'

  fn='l2cl_kappa_CMBLenS.dat_LOS'
  !fn='l2cl_kappa_nz_mass_ngp.dat_LOS'
  !fn='l2cl_kappa_nz_ngp.dat_LOS'
  !fn='l2cl_kappa_nz.da_LOS'

  open(30,file=Lens_output_path//trim(fn)//trim(LOS_str))
  do i=1,npc/2  ! Remove the '/2' for zeropad maps
     write(30,*) i*(2*pi/angle - 1),power(i)
  enddo
  close(30)
  write(*,*) 'Got l2Cl'

#endif

end program GetNzMaps
