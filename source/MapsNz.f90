

program GetNzMaps
  implicit none
  include 'Lens.fh'
  

  ! read these parameters:
  !real, parameter :: omegam = 0.2905, h = 0.6898, lbox = 505., omegav = 0.7095, w_de = -1.0

 
  !---------
  ! SLICS-LE 10deg
  !integer, parameter:: nc = 12288, nslice=18, npc=7745 ! For KiDSLenS
  character(*), parameter :: Maps_input_path2= './MapsDir/'
  character(*), parameter :: Maps_input_path= './MapsDir/'
  character(*), parameter :: Maps_output_path= './MapsDir/'
  
  real, dimension(npc,npc,nslice) :: map_3D
  real, dimension(npc,npc) :: map_in_out
  real, dimension(npc*2, npc*2) :: zeropad_map
  real, dimension(npc) :: power

  real kernel, z_write, z_source(nslice), z_lens(nslice), chi_lens(nslice), chi, chi_wde, Nz_CFHT, dz, gw, angle, z_source_plane, pi, Nz_buffer(1)
  real z_tmp, chi_l, chi_h, z_mean
  integer i, ii, zmin, cur_z, cur_z_s, Nz_col
  character(len=7) z_string,LOS_str,z_source_plane_str, z_lens_plane_str, Nz_col_str
  character(len=180) :: fn, Nz_filename, fn_out
  integer, parameter :: nzbins = 300 !501
  real, parameter :: Norm = 1.5
  real, dimension(nzbins) :: z_file, Nz, dchi, chi_file



#ifdef gfort  
  character(*), parameter :: output_form=   'unformatted'
  character(*), parameter :: output_access= 'stream' 
  character(*), parameter :: input_form=    'unformatted' 
#else
  character(*), parameter :: output_form=   'binary'
  character(*), parameter :: output_access= 'sequential' 
  character(*), parameter :: input_form=    'binary' 
#endif

  !external chi

  !equivalence (map_3D, zeropad_map)

  call GetArg(1,LOS_str)

#ifdef z_slice
  call GetArg(2,z_source_plane_str) 
  call GetArg(3,z_lens_plane_str) 
  read (z_source_plane_str,'(f7.3)') z_source_plane
  write(*,*) 'Placing z_source plane at ' , z_source_plane
  ! specify a N(z) file to get the z_array and the dz/dchi...
  Nz_filename = '../PopulateHalos/pz_SRD/pz_true.txt'
  !stop
#else
  call GetArg(2,Nz_filename)
  call GetArg(3,Nz_col_str)
  call GetArg(4,fn_out)
  read (Nz_col_str,'(i2)') Nz_col
  print *, 'Reading tomo bin column ', Nz_col_str, Nz_col 
#endif
 
  pi=acos(-1.)
  angle = sqrt(Area)/180*pi  ! For UBC Lens

  !---------------------------------


  !Read N(z) from C. Heymans N(z) KiDS & CFHT (Kilbinger?):
  !open (11,file='NZ_KIDS_v1.1.dat', status = 'old')
  !do i = 1,70
  !   read (11, *) z_file(i), Nz(i), Nz_CFHT
  !   !write(*,*) z_file(i), Nz(i)
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
  !open(11,file='nofz_spec_w_01_03.dat', status='old')
  !open(11,file='nofz_spec_w_03_05.dat', status='old')
  !open(11,file='nofz_spec_w_05_07.dat', status='old')
  !open(11,file='nofz_spec_w_07_09.dat', status='old')
  !open(11,file='nofz_spec_w_01_09.dat', status='old')

  !or the LSST-SRD n(z) :
  !# z full tomo0 tomo1 tomo2 tomo3 tomo4
  !open (11,file='../PopulateHalos/pz_SRD/pz_true.txt', status = 'old')
  open (11,file=Nz_filename, status = 'old')
  do i = 1,nzbins
     read (11, *) z_file(i), Nz_buffer(:)
     Nz(i) = Nz_buffer(Nz_col)
     !read (11, *)  Nz(i) ! Fit2KiDS
     write(*,*) z_file(i), Nz(i)
  end do
  close(11)

  !stop

  dz = z_file(2) - z_file(1)
  dchi = 0;

  !******
  ! z_file needs to be the center of the bin.
  !z_file = z_file + dz/2.0 ! z_file is lower edge of the bin, not centre
  !z_file = z_file - dz/2.0 ! z_file is upper edge of the bin, not centre
  !write(*,*)'z_file = ', z_file
  !write(*,*) 'sum of N(z) = ', sum(Nz*dz)
  !stop

  !---------------
  ! Test chi code:
  z_tmp = 3.0
  write(*,*) 'Old Chi Code:',chi(z_tmp,omegam,h)*h
  write(*,*) 'New Chi Code:'
  call chi_python(z_tmp, omegam, omegav, w_de, h, chi_wde, 1)
  print *, chi_wde
  !stop
  !--------------

  write(*,*) 'Filling the chi_file array:'
  call chi_python(z_file(:),omegam, omegav, w_de, h, chi_file(:), nzbins)
  write(*,*)'Done'
  !print *, chi_file
  
  !do i = 1,nzbins
  !   call chi_python(z_file(i), omegam, omegav, w_de, h, chi_file(i),1)
  !enddo
  

  ! first bin goes from zero to between first and second point:
  !call chi_python(z_file(1), omegam, omegav, w_de, h, chi_l)  
  !call chi_python(z_file(2), omegam, omegav, w_de, h, chi_h)  
  !dchi(1) = 0.5*(chi_l+ chi_h)

  do i = 1,nzbins-1
  !  !other bins go from mid points between each points:
  !  call chi_python(z_file(i), omegam, omegav, w_de, h, chi_l)  
  !  call chi_python(z_file(i+1), omegam, omegav, w_de, h, chi_h)  
    dchi(i) = chi_file(i+1) - chi_file(i)
  !  !dchi(i) = 0.5*(chi(z_file(i+1), omegam, h)*h - chi(z_file(i-1), omegam, h)*h)
  !   dchi(i) = chi_file(i+1) - chi_file(i)
  enddo
  ! last bin is given the distance between the last two points...
  !call chi_python(z_file(nzbins), omegam, omegav, w_de, h, chi_l)  
  call chi_python(z_file(nzbins)+dz, omegam, omegav, w_de, h, chi_h,1)  
  dchi(nzbins) = chi_h - chi_file(nzbins)

  !write(*,*) dchi
  !stop


#ifdef z_slice  
  Nz(:)=0
  do i = 1,nzbins
    if(z_file(i) .ge. z_source_plane) then  
       !Nz(i-1)=1
       Nz(i)=1
       exit
    endif
  enddo
#endif

  ! transform with Jacobian
  Nz = Nz*dz/dchi

  ! Normalize
  Nz = Nz/sum(Nz*dchi)

  z_mean = 0.0
  do i=1,nzbins
    z_mean = z_mean + z_file(i)*Nz(i)*dchi(i) 
  enddo 

  write(*,*) '**********'
  write(*,*) 'Nz stats: '
  write(*,*) '**********'
  write(*,*) 'Sum N(chi)dchi=', sum(Nz*dchi)
  write(*,*) 'Mean z =', z_mean


  !Read lens/source distributions
  !open (11,file='List_redshifts', status = 'old')
  !do i = 1,nslice
  open (11,file='checkpoints_cosmo', status = 'old')
  do i = nslice,1,-1
     read (11, *) z_string
     read (z_string,'(f7.4)') z_lens(i)
     !read (z_string,'(f7.3)') z_lens(i)
     chi_lens(i) = lbox/2.0*i - lbox/4.0
     ! round up
     z_lens(i) = (nint(z_lens(i)*1000.0)/1000.0)
     write(*,*) 'checkpoints: ', z_lens(i), chi_lens(i)
  end do
  close(11)
  !open (11,file='List_redshifts_source', status = 'old')
  !do i = 1,nslice
  open (11,file='checkpoints_cosmo_zs', status = 'old')
  do i = nslice,1,-1
     read (11, *) z_string
     read (z_string,'(f7.3)') z_source(i)
  end do
  close(11)

  !stop 


  ! Get maps
  write(*,*) 'Reading maps...'        

  do i = 1,nslice

        !Get Redshift
        z_write = z_lens(i)
        !write(*,*) 'Reading z=', z_write,'and LOS ', LOS_str
 
#ifdef read_delta 
        write (fn,'(f5.3,"delta.dat_bicubic")') z_write
        open(10,file=Maps_input_path//trim(fn)//'_LOS'//trim(LOS_str), access = output_access, form = output_form, status = 'old')
        read(10) map_in_out(:,:)
        close(10)
        !write(*,*)'Read delta Map for z=',z_write
#endif

#ifdef read_kappa
        write (fn,'(f5.3,"kappa_weight.dat")') z_write
        open(10,file=Maps_input_path//trim(fn)//'_LOS'//trim(LOS_str), access = output_access, form = output_form, status = 'old')
        read(10) map_in_out(:,:)
        close(10)
        !write(*,*)'Read kappa Map for z=',z_write
#endif
!----------------
#ifdef write_shear1
        write (fn,'(f5.3,"gamma1_weight.dat")') z_write
        open(11,file=Maps_input_path//trim(fn)//'_LOS'//trim(LOS_str), access = output_access, form = output_form, status = 'old')
        read(11) map_in_out(:,:)
        !write(*,*)'Read shear Maps for z=',z_write
        close(11)
#endif

#ifdef write_shear2
        write (fn,'(f5.3,"gamma2_weight.dat")') z_write
        open(12,file=Maps_input_path//trim(fn)//'_LOS'//trim(LOS_str), access = output_access, form = output_form, status = 'old')
        read(12) map_in_out(:,:)
        close(12)
        !write(*,*)'Read shear Maps for z=',z_write
#endif
#ifdef write_defl1
        write (fn,'(f5.3,"deflx_weight.dat")') z_write
        open(13,file=Maps_input_path//trim(fn)//'_LOS'//trim(LOS_str), access = output_access, form = output_form, status = 'old')
        read(13) map_in_out(:,:)
        close(13)
        !write(*,*)'Read shear Maps for z=',z_write
#endif
#ifdef write_defl2
        write (fn,'(f5.3,"defly_weight.dat")') z_write
        open(14,file=Maps_input_path//trim(fn)//'_LOS'//trim(LOS_str), access = output_access, form = output_form, status = 'old')
        read(14) map_in_out(:,:)
        close(14)
        !write(*,*)'Read shear Maps for z=',z_write
#endif

        map_3D(:,:,i)= map_in_out(:,:)

        !write(*,*) '**********'
        !write(*,*) 'Map stats: '
        !write(*,*) '**********'
        !write(*,*) sum(map_3D(:,:,i))
        !write(*,*) minval(map_3D(:,:,i)), maxval(map_3D(:,:,i))
  enddo

  write(*,*) 'Done reading maps'

#ifndef read_delta 
  !-------------------------
  ! Get the un-weighted maps

  ! start with first plane
  cur_z = 1
  cur_z_s = 1

  kernel = 3./2. * omegam * lbox/real(nc)/(3e3**2)*chi(z_lens(cur_z), omegam,h)*h * (1+z_lens(cur_z)) * &
        (1. - chi(z_lens(cur_z),omegam, h )/chi(z_source(cur_z_s), omegam, h))
  map_3D(:,:,cur_z_s) = map_3D(:,:,cur_z_s)/kernel

  !write(*,*) '**********'
  !write(*,*) 'Map stats: '
  !write(*,*) '**********'
  !write(*,*) sum(map_3D(:,:,cur_z_s))
  !write(*,*) minval(map_3D(:,:,cur_z_s)), maxval(map_3D(:,:,cur_z_s))
  !write(*,*) 'Kernel1 = ', kernel

  do cur_z_s = 2,nslice
     do cur_z = 1,cur_z_s - 1

        ! for each source, undo the sum over all lenses:
        kernel = 3./2. * omegam * lbox/real(nc)/(3e3**2)*chi(z_lens(cur_z), omegam,h)*h * (1+z_lens(cur_z)) * &
           (1. - chi(z_lens(cur_z),omegam, h )/chi(z_source(cur_z_s), omegam, h))
        !write(*,*) 'Kernel2 = ', kernel
    
        !write(*,*) '**********'
        !write(*,*) 'Map stats: '
        !write(*,*) '**********'
        !write(*,*) sum(map_3D(:,:,cur_z_s))
        !write(*,*) minval(map_3D(:,:,cur_z_s)), maxval(map_3D(:,:,cur_z_s))

        map_3D(:,:,cur_z_s) = map_3D(:,:,cur_z_s) - map_3D(:,:,cur_z)*kernel

        !write(*,*) '**********'
        !write(*,*) 'Map stats: '
        !write(*,*) '**********'
        !write(*,*) sum(map_3D(:,:,cur_z_s))
        !write(*,*) minval(map_3D(:,:,cur_z_s)), maxval(map_3D(:,:,cur_z_s))
!stop
     enddo

     ! then divide off by the final kernel, closest to the source
     cur_z = cur_z_s
     kernel = 3./2. * omegam * lbox/real(nc)/(3e3**2)*chi(z_lens(cur_z), omegam,h)*h * (1+z_lens(cur_z)) * &
         (1. - chi(z_lens(cur_z),omegam, h )/chi(z_source(cur_z_s), omegam, h))
     !write(*,*) 'Kernel3 = ', kernel
     map_3D(:,:,cur_z_s) = map_3D(:,:,cur_z_s)/kernel

     !write(*,*) '**********'
     !write(*,*) 'Map stats: '
     !write(*,*) '**********'
     !write(*,*) sum(map_3D(:,:,cur_z_s))
     !write(*,*) minval(map_3D(:,:,cur_z_s)), maxval(map_3D(:,:,cur_z_s))

  enddo
  write(*,*)' Got the delta planes'  
#else
  write(*,*)' Had the delta planes already'  
#endif

  map_in_out = 0

  !----------------
  ! Loop over lenses
  do i = 1,nslice

     cur_z_s = 1

     !find minimum for the Nz integration
     zmin = 0
     do ii = 1,nzbins
        if(z_file(ii) > z_lens(i))then
           zmin = ii
           exit
        endif
     enddo
     

     ! Get integal over sources
     gw = 0
     do ii=zmin, nzbins
        gw = gw + Nz(ii)*(1.0 - chi_lens(i)/chi_file(ii))*dchi(ii)
        !gw = gw + Nz(ii)*(1.0 - chi(z_lens(i),omegam,h)/chi(z_file(ii), omegam, h))*dchi(ii)
     enddo

     write(*,*) z_lens(i), zmin, z_file(zmin), Nz(zmin), gw
    
     ! Get Kernel
     !kernel = 3./2. * omegam * lbox/real(nc)/(3e3**2)*chi(z_lens(i), omegam,h)*h * (1+z_lens(i)) * gw
     kernel = 3./2. * omegam * lbox/real(nc)/(3e3**2)*chi_lens(i) * (1+z_lens(i)) * gw * Norm
     !write(*,*) 'Kernel4 = ', kernel


     ! Integrate
     map_in_out(:,:) = map_in_out(:,:) + map_3D(:,:,i)*kernel

  enddo 

  !--------------
  ! write to file


#ifdef write_kappa
  !write(*,*) 'Opening  map at z_lens  = ' , trim(z_lens_plane_str)
  !write (fn,'(f5.3,"kappa_weight.dat")') z_write
  !fn = 'kappa_nz_mass.dat'
  !fn = 'kappa_LSST-SRD_tomo1.dat'
  !fn = 'kappa_KiDS450_tomo2.dat'
  !fn = 'kappa_KiDS450_tomo3.dat'
  !fn = 'kappa_KiDS450_tomo4.dat'
  !fn = 'kappa_KiDS450_tomo0.dat'
#ifdef z_slice 
  fn=trim(z_lens_plane_str)//'kappa_weight.dat'
  !write (fn,'(f5.3,"kappa_weight.dat")') z_write
#else
  fn = fn_out
#endif
#endif


#ifdef write_shear1
  fn = 'gamma1_nz_mass.dat'
  !fn = 'gamma1_nz.dat'
#endif
#ifdef write_shear2
  fn = 'gamma2_nz_mass.dat'
  !fn = 'gamma2_nz.dat'
#endif

write(*,*) 'Opening map file ', fn

#if defined(write_shear1) || defined(write_shear2) || defined(write_kappa)
  open(11,file=Maps_output_path//trim(fn)//'_LOS'//trim(LOS_str), access = output_access, form = output_form)
  write(11) map_in_out(:,:)
  write(*,*)'Wrote Map for LOS', trim(LOS_str)
  close(11)
  !write(*,*) 'Commented out the writing of maps to file, doing simply l2Cl...'
#else
write(*,*)'No maps written'
#endif  


!--------------------
#ifdef power_spectrum

  write(*,*) 'starting l2Cl...'

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

  write(*,*)'calling ps2'
  call ps2_deconvolve(map_in_out,map_in_out,power,npc)
  
  write(*,*)'Done ps2'

  fn='l2cl_'//fn_out
  !fn='l2cl_kappa_nz_LSST-SRD_tomo'//trim(Nz_col-1)//'.dat_LOS'

  !open(30,file=Maps_output_path//trim(z_lens_plane_str)//trim(fn)//trim(LOS_str))
  !open(30,file=Maps_output_path//trim(z_source_plane_str)//trim(fn)//trim(LOS_str))
  open(30,file=Maps_output_path//trim(fn)//trim(LOS_str))
  do i=1,npc/2  ! Remove the '/2' for zeropad maps
     write(30,*) i*(2*pi/angle - 1),power(i)
  enddo
  close(30)
  write(*,*) 'Got l2Cl'

#endif

end program GetNzMaps
