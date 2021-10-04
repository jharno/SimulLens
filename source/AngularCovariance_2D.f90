!! Angular Covariance written by Joachim Harnois- Deraps
!! Power-spectrum modified on Hy Trac's + Ting Ting Lu's code.
program cicpow
  implicit none
#include 'BAO.fh'

  !! np should be set to nc (1:1) or hc (1:2)
  integer, parameter :: nt=1

#ifndef HALO
#ifdef subsample
  integer(kind=8), parameter :: np_max = hc**3 
  integer(kind=8) :: np
  real :: mp
  real, dimension(6,np_max) :: xv
  integer, dimension(np_max) :: ll
#else
  integer(kind=8), parameter :: np=hc**3
  real, parameter :: mp=int(nc,kind=8)**3/np !nc**3/np
  real, dimension(6,np) :: xv
  integer, dimension(np) :: ll
#endif

#ifdef GETPID
  integer(kind=8), dimension(np) :: PID
#endif

#ifdef GAUSS
    real, dimension(hc) :: ps_gauss
#endif

#else
  integer(kind=8), parameter :: np_max = hc**3 
  integer(kind=8) :: np  
  real :: mp
  real, dimension(6,np_max) :: xv
  integer, dimension(np_max) :: ll
#ifdef HALOPID  
  integer,parameter :: N_p = 10
  integer(kind=8), dimension(N_p,10000) :: halo_pid
#endif
#endif



  !! cubepm
!  integer, parameter :: nodes_dim=2
  real, parameter  :: ncc=nc/nodes_dim !! number of cells / cubic 
  real, parameter  :: rnc = nc

  !! Dark matter arrays
  integer, dimension(nc,nc) :: hoc
  integer, dimension(2,nc,nc,nt) :: htoc

  !! Power spectrum arrays
  real, dimension(nc) :: ps,err
  real, dimension(nc) :: PS_Ave!,  CrossCorrAve, N_CrossCorrAve
  real, dimension(3*(nc/2)**2+1) :: CrossCorrAve, N_CrossCorrAve, kr2weight
  real, dimension(nc+2,nc,nc) :: d, CrossCorr, N_CrossCorr
  real, dimension(nc+2,nc,nc) ::  psShell_Inner, psShell_Outer, N_psShell_Inner, N_psShell_Outer
  complex, dimension(nc/2 + 1,nc,nc) ::  d_COMPLEX, psV_COMPLEX, CrossCorr_COMPLEX, N_CrossCorr_COMPLEX 
  complex, dimension(nc/2 + 1,nc,nc) ::  psShell_Inner_COMPLEX, psShell_Outer_COMPLEX, N_psShell_Inner_COMPLEX,N_psShell_Outer_COMPLEX
  
  !equivalence (psShell_Inner,psShell_Inner_COMPLEX,CrossCorr, CrossCorr_COMPLEX)
  !equivalence (psShell_Outer,psShell_Outer_COMPLEX)
  !equivalence (N_psShell_Inner,N_psShell_Inner_COMPLEX)
  !equivalence (d, d_COMPLEX, psV_COMPLEX, N_CrossCorr,N_psShell_Outer,N_psShell_Outer_COMPLEX, N_CrossCorr_COMPLEX) 
  
  equivalence (psShell_Inner,psShell_Inner_COMPLEX,CrossCorr, CrossCorr_COMPLEX)
  equivalence (psShell_Outer,psShell_Outer_COMPLEX)
  equivalence (N_psShell_Inner,N_psShell_Inner_COMPLEX)
  equivalence (d, d_COMPLEX, psV_COMPLEX)
  equivalence (N_CrossCorr,N_psShell_Outer,N_psShell_Outer_COMPLEX, N_CrossCorr_COMPLEX) 
 

#ifdef FFT_TEST
 real, dimension(nc+2,nc,nc) :: Sync, Gauss, DeltaFunction, N_Hat, Hat
 real, dimension(nc) :: Xsi, HatAve, GaussAve, DeltaAve, FFT_Hat
#endif

  !real, dimension(0:hc,-hc+1:hc,-hc+1:hc) :: Hat
  integer, parameter :: MSL=150
  character (len=MSL) :: ofile,psFile
  character(len=MSL):: MAX_Inner_char, MAX_Outer_char, Thickness_Inner_char,Thickness_Outer_char,Run,Version

  integer nploc(nn),fstat
  real pi,omegak,MAX_Inner, MAX_Outer, Thickness_Inner, Thickness_Outer, dummy1, dummy2

  !! variables in CubePM
  integer cubepm_nts,cubepm_cur_checkpoint,cubepm_cur_projection,cubepm_cur_halofind,i1,j1,k1, node_coords(3)
  real  cubepm_a,cubepm_t,cubepm_tau,cubepm_dt_f_acc,cubepm_dt_pp_acc,cubepm_dt_c_acc, cubepm_mass_p

#ifdef HALO
  real, dimension(np_max) :: cubepm_v_disp, cubepm_radius_scale, halo_mass, halo_mass_pp, halo_mass1
  real, dimension(3,np_max) :: cubepm_halo_pos, peak, cubepm_l, var
#endif


  real, dimension(hc) :: PSFromFile


  common /rarr/ d,xv,ps 
  common /iarr/ ll,htoc,hoc

  pi=acos(-1.)
  omegak = 1-omegam-omegav

  call getarg(1,Run)
  call getarg(2,Version)



#ifndef RANDOM
#ifndef LegendreTest

#ifdef subsample
  call readdm_subsample 
#else
  call readdm
#endif

  call cic
#ifdef dump_density
  !write(*,*)'TMP::WRINTING DENSITY TO FILE'
  !open(unit = 33,file = '/mnt/scratch-3week/jharno/Densities/density.dat')
  ofile = '/mnt/scratch-3week/jharno/Densities/density-'//trim(Run)//'-'//trim(Version)//'.dat'
  open (unit=33,file=ofile,status='replace',iostat=fstat,form='binary')
  !open (unit=33,file=ofile,status='replace',iostat=fstat)

  if(fstat /= 0) then 
        write(*,*) 'could not open the file properly' 
        write(*,*) 'ofile =' , ofile
        write(*,*) 'max(d) =' , maxval(d) 
  else 
        write(*,*) 'WRITING DENSITY TO FILE AND STOP'
        write(33) d
        !do k1 = 1,nc
        !   do j1 = 1,nc
        !      do i1 = 1,nc
        !         !write(*,*)i1,j1,k1,d(i1,j1,k1)
        !         write(33, '(1f18.8)') d(i1,j1,k1)
        !      enddo
        !   enddo
        !enddo
        write(*,*) 'Wrote density file'
                 
  endif
     !write(*,*) 'could not open file properly'
  !endif
  close(33)
  !write(*,*)'wrote density file'
  stop
#endif
#endif
#endif

#ifdef GETPID
  call PID_analysis
#endif
#ifdef HALOPID
  call PID_analysis
#endif


!#ifndef LegendreTest
  call powerspectrum
  call writeps
!#endif



#ifdef shotnoise
  call random_number(xv(1:3,1:np))
  xv(1:3,1:np)=xv(1:3,1:np)*nc
  call cic
  call powerspectrum
  call writeps_rdm
  stop
#endif


#ifdef AngularCovariance


    ! Get P(k)
    ! Now, after calling powerspectrum, d is already complex, and equivalenced to psV_COMPLEX
    ! The following two ways are equivalent.

     psV_COMPLEX = lbox**3*abs(psV_COMPLEX/nc**3)**2

     !do i1 = 1, nc+2,2
     !   d(i1,:,:)=((d(i1,:,:)/nc**3)**2+(d(i1+1,:,:)/nc**3)**2)*lbox**3
     !   d(i1+1,:,:)=0
     !enddo

     !write(*,*) 'psV_COMPLEX = ', psV_COMPLEX(2,1,1),  psV_COMPLEX(1,2,1),  psV_COMPLEX(1,1,2)



     !Choose file to subtract the mean
 
#ifdef poisson  
     !psFile = '/home/jharno/AngularCovariance/ngppow_z05_N10-poisson-8000.txt'
     psFile = '/mnt/scratch-3week/jharno/AngularCovariance/data/PowerSpectrum/z_0.5/ngppow_z05_N10-poisson-8000000.txt' 
#else
#ifdef HALO
     !psFile = '/home/jharno/AngularCovariance/ngppow_z05_N10-halo.txt' 
     psFile = '/mnt/scratch-3week/jharno/AngularCovariance/data/PowerSpectrum/z_0.5/ngppow_z05_N10-halo.txt'
#else
#ifdef subsample
     !psFile = '/home/jharno/AngularCovariance/ngppow_z05_N10-8000.txt'
     psFile = '/mnt/scratch-3week/jharno/AngularCovariance/data/PowerSpectrum/z_0.5/ngppow_z05_N10-8000.txt'
#else
#ifdef LegendreTest
     psFile = '/mnt/scratch-3week/jharno/AngularCovariance/data/PowerSpectrum/z_0.5/ngppow-Legendre.dat'
#else
     !psFile = '/home/jharno/AngularCovariance/ngppow_z05_N10-wayne.txt'
     !psFile = '/home/jharno/AngularCovariance/ngppow_z05_N300-wayne.txt' 
     !psFile = '/mnt/scratch-3week/jharno/AngularCovariance/data/PowerSpectrum/z_0.5/ngppow_z05_N10.txt'
     !psFile = '/home/jharno/AngularCovariance/ngppow_z05_N200.txt'
     psFile = './ngppow_z05_N200.txt'
     !psFile = '/home/jharno/AngularCovariance/meanP_Z0-NGP.txt'
#endif
#endif
#endif
#endif


  !else
  !   write(*,*) 'wrong PSFromFile!'
  !endif


#ifndef GetPowerShells
  write(*,*) 'opening', psFile
  open(unit=24,file=psFile)
  do i1=1,hc/2

#ifdef LegendreTest
     read(24,*) dummy1,PSFromFile(i1), dummy2 
#else
     read(24,*) PSFromFile(i1)
#endif

#ifdef debug_PSFromFile
     write(*,*) i1, PSFromFile(i1)
#endif 

  enddo
  close(24)
#endif

#ifdef GAUSS
  do i1 = 1,hc
     PSFromFile(i1) = ps_gauss(i1)/((2*pi*real(i1)/lbox)**3)*2*pi**2
  enddo
  write(*,*) 'Assigned ps_gauss to PSFromFile'
#endif

#ifdef OverwritePower
  do i1 = 1,hc
     PSFromFile(i1) = ps(i1)/((2*pi*real(i1)/lbox)**3)*2*pi**2
  enddo
  write(*,*) 'Assigned ps to PSFromFile'
  write(*,*) 'PSFromFile(i1) =', PSFromFile(1)
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !call getarg(3,MAX_Inner_char)
    !call getarg(4,MAX_Outer_char)
    !call getarg(3,Thickness_Inner_char)
    !call getarg(4,Thickness_Outer_char)
    Thickness_Inner_char = '1'
    Thickness_Outer_char = '1'

    !write(*,*)'Inputs:',trim(MAX_Inner_char),trim(MAX_Outer_char),trim(Thickness_Inner_char),trim(Thickness_Outer_char)

    !do i1 = 1,128 
    !do i1 = 1,10
    do i1 =  31,32

       write(MAX_Outer_char,'(i4)') i1
       MAX_Outer_char=adjustl(MAX_Outer_char)

       !do j1 = i1,1,-1
       !do j1 = i1, 1,-1
       do j1 = i1,31,-1

          write(MAX_Inner_char,'(i4)') j1
          MAX_Inner_char=adjustl(MAX_Inner_char)


          read(MAX_Inner_char,*) MAX_Inner
          read(MAX_Outer_char,*) MAX_Outer
          read(Thickness_Inner_char,*) Thickness_Inner
          read(Thickness_Outer_char,*) Thickness_Outer

    !write(*,*)'MAX_Inner=', MAX_Inner
    !write(*,*)'MAX_Outer=', MAX_Outer
    !write(*,*)'Thickness_Inner=', Thickness_Inner
    !write(*,*)'Thickness_Inner=', Thickness_Inner

          call CrossCorrelationOfShells
          call writeCov
       enddo
    enddo
#endif

contains

  subroutine readdm
    implicit none
    character*100 fn, fn2, fn3

    character*8 t1,t2
    integer i,j, Nmax, nh_local, nh_total
    integer(kind=8) ip
    real HubbleScale
    real Conversion
    character (len=4) :: rank_s
    character (len=MSL) :: ofile,zstring,PIDofile
    character (len=MSL) :: ifile,PIDifile
#ifdef HALO
    real(4), dimension(20) :: halo_input_buffer
#endif

#ifdef MERGED_LIST
    real, dimension(17):: reading_buffer
#endif

#ifndef MERGED_LIST

    !! Read particle data file from decomposed volume
    ip=0
    nh_total=0       
    do i=1,nn

       nh_local = 0

       write(*,*) 'Reading Node ',i
       write(rank_s,'(i4)') i-1
       rank_s=adjustl(rank_s)
       write(zstring,'(f5.3)') z3dps       ! Need (f6.3) for Z > 10

#ifdef HALO
       ifile=trim(zstring)//'halo'//rank_s(1:len_trim(rank_s))//".dat"
#else 
       ifile=trim(zstring)//'xv'//rank_s(1:len_trim(rank_s))//".dat"
#endif

#ifdef IC
       ifile='xv'//rank_s(1:len_trim(rank_s))//".ic"
#endif

       !ofile='/mnt/scratch-3week/jharno/Lensing/LOS6/'//trim(ifile)
       !ofile='/mnt/scratch-2week/jharno/power_law/z_i=20/'//trim(ifile)
       !ofile='/mnt/scratch-2week/jharno/PID/'//trim(ifile)
       !ofile=proj_path//trim(Version)//'/out/RUN-'//trim(Run)//'/'//trim(ifile)
       !ofile='/mnt/scratch-3week/ngan/sims/out/'//trim(Version)//'/'//trim(ifile)
       !ofile='/mnt/scratch-3week/akhazr/output/box200/'//trim(Version)//'/'//trim(ifile)
       ofile='/mnt/scratch-3week/jharno/Simulations/200Mpc/Run-'//trim(Version)//'/'//trim(ifile)

!**********
       write(*,*) 'opening ',ofile
       open (unit=12,file=ofile,status='old',iostat=fstat,form='binary')
       if (fstat /= 0) then
          write(*,*) 'error opening catalog'
          write(*,*) 'rank=',nn, 'file:',ofile
          stop !call mpi_abort(mpi_comm_world,ierr,ierr)
       endif
       


#ifdef HALO
       read(12) nploc(i)
       write(*,*) 'nploc(i) =',nploc(i)
       do
          nh_total=nh_total+1
          nh_local=nh_local+1          
          !write(*,*)'Reading halo', nh_total

#ifdef HALOPID
          read(12,end=112,err=113) halo_input_buffer, halo_pid(:,nh_total)
#else
          read(12,end=112,err=113) halo_input_buffer
#endif
          !Assign variable to buffer
          !write(*,*) halo_input_buffer, halo_pid(:,nh_total)
          peak(:,nh_total)=halo_input_buffer(1:3)
          xv(:,nh_total)=halo_input_buffer(4:9)
          cubepm_l(:,nh_total)=halo_input_buffer(10:12)
          cubepm_v_disp(nh_total)=halo_input_buffer(13)
          cubepm_radius_scale(nh_total)=halo_input_buffer(14)
          halo_mass(nh_total)=halo_input_buffer(15)
          halo_mass_pp(nh_total)=halo_input_buffer(16)
          halo_mass1(nh_total)=halo_input_buffer(17)
          var(:,nh_total)=halo_input_buffer(18:20)

#ifdef debug_halo
          write(*,*)'peak=',peak(:,nh_total)
          write(*,*)'xv=',xv(:,nh_total)
          write(*,*)'cubepm_l=',cubepm_l(:,nh_total)
          write(*,*)'cubepm_v_disp=',cubepm_v_disp(nh_total)
          write(*,*)'cubepm_radius_scale=',cubepm_radius_scale(nh_total)
          write(*,*)'halo_mass=',halo_mass(nh_total)
          write(*,*)'halo_mass_pp=',halo_mass_pp(nh_total)
          write(*,*)'halo_mass1=',halo_mass1(nh_total)
          write(*,*)'var=',var(:,nh_total)
#ifdef HALOPID
          write(*,*)'PID=',halo_pid(:,nh_total)   
#endif
#endif

113       continue
       enddo
112    close(12)

       !remove the last increment. Reading empty halo
       nh_total = nh_total-1
       nh_local = nh_local-1
       !nploc(i) = nh_total
       write(*,*)'but found EOF. Working with ',nh_total
            

#else
#ifdef IC 
!#ifdef pmfast
       read(12) nploc(i)
#endif
#ifdef cubepm
       read(12) nploc(i), cubepm_a,cubepm_t,cubepm_tau,cubepm_nts,cubepm_dt_f_acc,cubepm_dt_pp_acc,cubepm_dt_c_acc,cubepm_cur_checkpoint, cubepm_cur_projection,cubepm_cur_halofind,cubepm_mass_p
#endif
       read(12) xv(1:6,ip+1:ip+nploc(i))
       close(12)
#endif

!***********

#ifdef pmfast       
       xv(3,ip+1:ip+nploc(i))=xv(3,ip+1:ip+nploc(i))+(i-1)*nc/nn
#endif
#ifdef cubepm
       do k1=1,nodes_dim
          do j1=1,nodes_dim
             do i1=1,nodes_dim
                if (i-1 == (i1-1)+(j1-1)*nodes_dim+(k1-1)*nodes_dim**2)  &
                    node_coords(:)=(/(i1-1),(j1-1),(k1-1)/)
             enddo
          enddo
       enddo
#ifndef HALO       
       xv(1,ip+1:ip+nploc(i))=modulo(xv(1,ip+1:ip+nploc(i))+node_coords(1)*ncc,rnc)
       xv(2,ip+1:ip+nploc(i))=modulo(xv(2,ip+1:ip+nploc(i))+node_coords(2)*ncc,rnc)
       xv(3,ip+1:ip+nploc(i))=modulo(xv(3,ip+1:ip+nploc(i))+node_coords(3)*ncc,rnc)
#endif
#endif
#ifdef Kaiser

    !Red Shift Distortion: x_z -> x_z +  v_z/H(Z)   
    !Converting seconds into simulation time units
    !cancels the H0...
    
    !xv(3,ip+1:ip+nploc(i))=xv(3,ip+1:ip+nploc(i)) + xv(6,ip+1:ip+nploc(i))*1.5*sqrt(omegam/cubepm_a)
    xv(3,ip+1:ip+nploc(i))=xv(3,ip+1:ip+nploc(i)) + xv(6,ip+1:ip+nploc(i))*1.5/sqrt(cubepm_a*(1+cubepm_a*omegak/omegam + omegav/omegam*cubepm_a**3))
    
    if(i==nn) then
       write(*,*) '**********************'
       write(*,*) 'Included Kaiser Effect'
       write(*,*) 'Omega_m =', omegam, 'a =', cubepm_a
       !write(*,*) '1/H(z) =', 1.5*sqrt(omegam/cubepm_a)
       write(*,*) '1/H(z) =', 1.5/sqrt(cubepm_a*(1+cubepm_a*omegak/omegam + omegav/omegam*cubepm_a**3))
       write(*,*) '**********************'
    endif
#endif

#ifdef GETPID
       PIDifile=trim(zstring)//'PID'//rank_s(1:len_trim(rank_s))//'.dat'
       !PIDofile=proj_path//trim(Version)//'/out/RUN-'//trim(Run)//'/'//trim(ifile)
       PIDofile='/cita/d/scratch-2week/jharno/PID/'//trim(PIDifile)

       write(*,*) 'opening ',PIDofile

       open (unit=19,file=PIDofile,form='binary')
       read(19) nploc(i), cubepm_a,cubepm_t,cubepm_tau,cubepm_nts,cubepm_dt_f_acc,cubepm_dt_pp_acc,cubepm_dt_c_acc,cubepm_cur_checkpoint, &
            cubepm_cur_projection,cubepm_cur_halofind,cubepm_mass_p       
       read(19) PID(ip+1:ip+nploc(i))
       
       close(19)

#endif

#ifndef HALO
       ip=ip+nploc(i)
       write(*,*) 'np cumulative = ', ip,', np local = ', nploc(i)
#else
       ip=ip+nh_local
       write(*,*) 'np cumulative = ', ip,'or', nh_total,', np local = ', nh_local       
#endif

    enddo

#else

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! If working with a MERGED LIST !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef HALO
    ip = 0

    write(zstring,'(f5.3)') z3dps
    ifile=trim(zstring)//'halo.dat' ! xv.dat
    ofile='/cita/d/scratch-2week/jharno/PID/'//trim(ifile)
    !ofile=proj_path//trim(Version)//'/out/RUN-'//trim(Run)//'/'//trim(ifile)

    write(*,*) 'opening ',ofile
    open (unit=12,file=ofile)

    do
       read(12,'(17f20.10)',end=112,err=113) reading_buffer
       ip = ip + 1           

       peak(1:3,ip) = reading_buffer(1:3)
       xv(1:6,ip) = reading_buffer(4:9)
       cubepm_l(1:3,ip) = reading_buffer(10:12)
       cubepm_v_disp(ip) = reading_buffer(13)
       cubepm_radius_scale(ip) = reading_buffer(14)
       halo_mass(ip) = reading_buffer(15)
       halo_mass_pp(ip) = reading_buffer(16)
       halo_mass1(ip) = reading_buffer(17)
113    continue
    enddo
112 close(12)

    np = ip
    
    write(*,*) 'n particle = ', np

    !test
    !fn3 = 'test_'//trim(zstring)//'halo.dat'
    !open(unit=22,file=fn3)
    !do i =1,np
    !   write(22,*) xv(1:3,i),cubepm_radius_scale(i),halo_mass(i) 
    !enddo
    !close(22)
    !write(*,*)'wrote', fn3


#ifdef Kaiser
    !Red Shift Distortion: x_z -> x_z +  v_z/H(Z)   
    !Converting seconds into simulation time units
    !cancels the H0...
    
    !xv(3,:)=xv(3,:) + xv(6,:)*1.5*sqrt(omegam/cubepm_a)
    xv(3,:)=xv(3,:) + xv(6,:)*1.5/sqrt(cubepm_a*(1+cubepm_a*omegak/omegam + omegav/omegam*cubepm_a**3))
    
    if(i==nn) then
       write(*,*) 'Included Kaiser Effect'
    endif
#endif



#else
    write(*,*) '***************************************************'
    write(*,*) '*** MERGE LIST NOT IMPLEMENTED FOR PARTICLES YET***'
    write(*,*) '***************************************************'
#endif

#ifdef debug 
    write(*,*) 'np = ', np
    write(*,*) 'xv', maxval(xv(1:3,:)),minval(xv(1:3,:)),maxval(xv(4:6,:)),minval(xv(4:6,:))
#endif


#endif


    

#ifdef HALO
    np = ip

    if(np==0)then
       write(*,*) 'No haloes found!'
       return
    endif

    !**********************************
    ! Reassigning the mass of the halos
    mp = real(nc**3)/real(np)
    !**********************************

    write(*,*) 'halo mass = ', mp
    write(*,*) 'total mass = ', mp*np
    write(*,*) 'total mass pp= ', 8*hc**3
    write(*,*) 'average mass per grid cell= ', mp*np/nc**3

#endif


    write(*,*) '*************'
    write(*,*) '*Done readdm*'
    write(*,*) '*************'
    return
  end subroutine readdm

  !***********************
  subroutine readdm_subsample
    implicit none
    character*100 fn, fn2, fn3

    character*8 t1,t2
    integer i,j, Nmax, nh_local, nh_total!, rdm_array(80000)
    integer(kind=8) ip, counter
    real HubbleScale, rdm, threshold
    real Conversion
    character (len=4) :: rank_s
    character (len=MSL) :: ofile,zstring
    character (len=MSL) :: ifile

    write(*,*) 'Read particle data file from decomposed volume'

    ip=0
    nh_total=0
    do i=1,nn

       nh_local = 0

       write(*,*) 'Reading Node ',i
       write(rank_s,'(i4)') i-1
       rank_s=adjustl(rank_s)
       write(zstring,'(f5.3)') z3dps       ! Need (f6.3) for Z > 10

#ifdef HALO
       !ifile=trim(zstring)//'halo'//rank_s(1:len_trim(rank_s))//".dat"
       write(*,*) 'Subsampling of halos not implemented yet'
       stop
#else 
       ifile=trim(zstring)//'xv'//rank_s(1:len_trim(rank_s))//".dat"
#endif

#ifdef IC
       ifile='xv'//rank_s(1:len_trim(rank_s))//".ic"
#endif

       !ofile='/mnt/scratch-3week/jharno/Lensing/LOS6/'//trim(ifile)
       !ofile='/mnt/scratch-2week/jharno/power_law/z_i=20/'//trim(ifile)
       !ofile='/mnt/scratch-2week/jharno/PID/'//trim(ifile)
       !ofile=proj_path//trim(Version)//'/out/RUN-'//trim(Run)//'/'//trim(ifile)
       !ofile='/mnt/scratch-3week/akhazr/output/box200/'//trim(Version)//'/'//trim(ifile)
       ofile='/mnt/scratch-3week/jharno/Simulations/200Mpc/Run-'//trim(Version)//'/'//trim(ifile)

!**********
       write(*,*) 'opening ',ofile
       open (unit=12,file=ofile,status='old',iostat=fstat,form='binary')
       if (fstat /= 0) then
          write(*,*) 'error opening catalog'
          write(*,*) 'rank=',nn, 'file:',ofile
          stop !call mpi_abort(mpi_comm_world,ierr,ierr)
       endif
       
#ifdef cubepm
       read(12) nploc(i), cubepm_a,cubepm_t,cubepm_tau,cubepm_nts,cubepm_dt_f_acc,cubepm_dt_pp_acc,cubepm_dt_c_acc,cubepm_cur_checkpoint, cubepm_cur_projection,cubepm_cur_halofind,cubepm_mass_p
#endif
       read(12) xv(1:6,ip+1:ip+nploc(i))
       close(12)


!***********

#ifdef pmfast       
       xv(3,ip+1:ip+nploc(i))=xv(3,ip+1:ip+nploc(i))+(i-1)*nc/nn
#endif
#ifdef cubepm
       do k1=1,nodes_dim
          do j1=1,nodes_dim
             do i1=1,nodes_dim
                if (i-1 == (i1-1)+(j1-1)*nodes_dim+(k1-1)*nodes_dim**2)  &
                    node_coords(:)=(/(i1-1),(j1-1),(k1-1)/)
             enddo
          enddo
       enddo
#ifndef HALO       
       xv(1,ip+1:ip+nploc(i))=modulo(xv(1,ip+1:ip+nploc(i))+node_coords(1)*ncc,rnc)
       xv(2,ip+1:ip+nploc(i))=modulo(xv(2,ip+1:ip+nploc(i))+node_coords(2)*ncc,rnc)
       xv(3,ip+1:ip+nploc(i))=modulo(xv(3,ip+1:ip+nploc(i))+node_coords(3)*ncc,rnc)
#endif
#endif

#ifndef HALO
       ip=ip+nploc(i)
       write(*,*) 'np cumulative = ', ip,', np local = ', nploc(i)
#else
       ip=ip+nh_local
       write(*,*) 'np cumulative = ', ip,'or', nh_total,', np local = ', nh_local
#endif

    enddo

    !write(*,*)'Opening random_array file'
    !open(unit=21, file='./random_array_80000', form='formatted')
    !   read(21,fmt='(i8)' ) rdm_array
    !close(21)

    write(*,*) 'picking xv from a random array'
    counter = 0
    threshold=8000.0/ip
    write(*,*) 'threshold = ', threshold
    call random_seed
    do i = 1,ip
      !call random_seed
      call random_number(rdm)
      if(rdm<(threshold)) then
        counter = counter + 1
        xv(:,counter)= xv(:,i)     
        !write(*,*)'Got one!', rdm, counter, ip, xv(:,counter)
      endif
 
      !write(*,*) i,rdm_array(i)
      !write(*,*) xv(1:3, i), xv(1:3, rdm_array(i))
      !xv(:,i) = xv(:,rdm_array(i))
      !write(*,*) xv(1:3, i), xv(1:3, rdm_array(i))
      !pause
    enddo

! I don't know why, but I need this line...
#ifdef subsample
    ! np = ip
    np = counter ! 8000
    !write(*,*) 'np = ', np
    if(np==0)then
       write(*,*) 'No particles found!'
       return
    endif


    !**************************************
    ! Reassigning the mass of the particles
    mp = real(nc**3)/real(np)
    !*************************************
#else
    write(*,*) 'Something is terribly wrong with this code'
#endif

    write(*,*) 'pp mass = ', mp
    write(*,*) 'total mass = ', mp*np
    write(*,*) 'total mass pp= ', np*hc**3
    write(*,*) 'average mass per grid cell= ', mp*np/nc**3
 
    return
  end subroutine readdm_subsample

  !************************
  subroutine writeps
    implicit none
    integer k
    character*250 fn1

#ifdef NGP_PS

    fn1=PS_dir//'z_'//Z//'/ngppow_'//RedShift//'_Run'//trim(Run)//'-'//trim(Version)//'.dat'   
#ifdef subsample 
    fn1=PS_dir//'z_'//Z//'/ngppow_'//RedShift//'_Run'//trim(Run)//'-'//trim(Version)//'-8000.dat'    
#endif
    !fn1=PS_dir//'z_'//Z//'/ngppow_'//RedShift//'_Run'//trim(Run)//'-'//trim(Version)//'.dat'
#ifdef TEST
    !fn1=Test_dir//'ngppow_'//RedShift//'_Run'//trim(Run)//'-'//trim(Version)//'.dat'
    fn1 = proj_path//trim(Version)//'/out/RUN-'//trim(Run)//'/'//Z//'ngpps.dat'
#endif

#ifdef HALO
    fn1=PS_dir//'z_'//Z//'/ngppow_'//RedShift//'_Run'//trim(Run)//'-'//trim(Version)//'-halo.dat'
#endif


#ifdef Kaiser
    fn1=PS_dir//'z_'//Z//'/npgpow_'//RedShift//'_Run'//trim(Run)//'-'//trim(Version)//'-RSD.dat'
#endif

#ifdef RANDOM
#ifdef poisson
    fn1=PS_dir//'z_'//Z//'/ngppow_'//RedShift//'_Run'//trim(Run)//'-'//trim(Version)//'-poisson-8000000.dat'
#else
    fn1=PS_dir//'z_'//Z//'/ngppow_'//RedShift//'_Run'//trim(Run)//'-'//trim(Version)//'-RDM.dat'
#endif  

#endif

#ifdef GAUSS
    fn1=PS_dir//'z_'//Z//'/ngppow_'//RedShift//'_Run'//trim(Run)//'-'//trim(Version)//'-GAUSS.dat'
#ifdef TEST
    fn1=Test_dir//'ngppow_'//RedShift//'_Run'//trim(Run)//'-'//trim(Version)//'-GAUSS.dat'
#endif
#endif

#ifdef power_law
    fn1=PS_dir//'z_'//Z//'/ngppow_'//RedShift//'_Run'//trim(Run)//'-'//trim(Version)//'-power_law_z_i=20.dat'
#endif

#ifdef LegendreTest
      fn1=PS_dir//'z_'//Z//'/ngppow-Legendre.dat'     
#endif

#else

    fn1=PS_dir//'z_'//Z//'/cicpow_'//RedShift//'_Run'//trim(Run)//'-'//trim(Version)//'.dat'
#ifdef TEST
    fn1=Test_dir//'cicpow_'//RedShift//'_Run'//trim(Run)//'-'//trim(Version)//'.dat'
    !fn1 = proj_path//trim(Version)//'/out/RUN-'//trim(Run)//'/'//Z//'00cicps.dat'
#endif

#ifdef HALO
    fn1=PS_dir//'cicpow_'//RedShift//'_Run'//trim(Run)//'-'//trim(Version)//'-halo.dat'
#endif


#ifdef Kaiser
    fn1=PS_dir//'z_'//Z//'/cicpow_'//RedShift//'_Run'//trim(Run)//'-'//trim(Version)//'-RSD.dat'
#endif

#ifdef RANDOM
    fn1=PS_dir//'z_'//Z//'/cicpow_'//RedShift//'_Run'//trim(Run)//'-'//trim(Version)//'-RDM.dat'
#endif

#ifdef GAUSS
    fn1=PS_dir//'z_'//Z//'/cicpow_'//RedShift//'_Run'//trim(Run)//'-'//trim(Version)//'-GAUSS.dat'
#ifdef TEST
    fn1=Test_dir//'cicpow_'//RedShift//'_Run'//trim(Run)//'-'//trim(Version)//'-GAUSS.dat'
#endif
#endif

#ifdef power_law
    fn1=PS_dir//'z_'//Z//'/cicpow_'//RedShift//'_Run'//trim(Run)//'-'//trim(Version)//'-power_law_z_i=20.dat'
#endif

#endif

    !! Output power spectrum
    !! First column is physical k
    !! Second column is \Delta^2

    open(11,file=fn1)
    do k=1,hc !hc+1
       write(11,*) 2*pi/lbox*(k), ps(k) ,err(k)
    enddo
    close(11)

    write(*,*) 'Wrote ', fn1

    write(*,*) '**************'
    write(*,*) '*done writeps*'
    write(*,*) '**************'

    return
  end subroutine writeps

  subroutine writeps_rdm
    implicit none
    integer k
    character*150 fn1

#ifdef HALO
    fn1=Test_dir//'halopow_'//RedShift//'_Run'//trim(Run)//'-'//trim(Version)//'-shotnoise.dat'
#else    
    fn1=PS_dir//'z_'//Z//'/cicpow_'//RedShift//'_Run'//trim(Run)//'-'//trim(Version)//'-shotnoise.dat'
#endif

    !! Output power spectrum
    !! First column is physical k
    !! Second column is \Delta^2

    open(11,file=fn1)
    do k=1,hc !hc+1
       write(11,*) 2*pi/lbox*(k), ps(k) 
    enddo
    close(11)

    write(*,*) 'Wrote ', fn1

    write(*,*) '******************'
    write(*,*) '*done writeps_rdm*'
    write(*,*) '******************'

    return
  end subroutine writeps_rdm

  !************************

  subroutine writeCov
    implicit none
    integer k,kr2
    character*350, fn2

    fn2=CrossCorr_dir//'z_'//Z//'/CrossCorr_'//RedShift//'_'//trim(MAX_Outer_char)//'_'//trim(MAX_Inner_char)//'_Run'//trim(Run)//'-'//trim(Version)//'.dat'

#ifdef subsample
    fn2=CrossCorr_dir//'z_'//Z//'/CrossCorr_'//RedShift//'_'//trim(MAX_Outer_char)//'_'//trim(MAX_Inner_char)//'_Run'//trim(Run)//'-'//trim(Version)//'-8000.dat'
#endif
#ifdef HALO
    fn2=CrossCorr_dir//'z_'//Z//'/CrossCorr_'//RedShift//'_'//trim(MAX_Outer_char)//'_'//trim(MAX_Inner_char)//'_Run'//trim(Run)//'-'//trim(Version)//'-halos.dat'
#endif

#ifdef MNT_SCRATCH
    fn2='/mnt/node_scratch/jharno/data/CrossCorr_'//RedShift//'_'//trim(MAX_Outer_char)//'_'//trim(MAX_Inner_char)//'_Run'//trim(Run)//'-'//trim(Version)//'.dat'
#endif

#ifdef TEST
    fn2=Test_dir//'CrossCorr_'//RedShift//'_'//trim(MAX_Outer_char)//'_'//trim(MAX_Inner_char)//'_Run'//trim(Run)//'-'//trim(Version)//'.dat'
#endif

#ifdef Kaiser
    fn2=CrossCorr_dir//'z_'//Z//'/CrossCorr_'//RedShift//'_'//trim(MAX_Outer_char)//'_'//trim(MAX_Inner_char)//'_Run'//trim(Run)//'-'//trim(Version)//'-RSD.dat'
#endif

#ifdef RANDOM
#ifdef poisson
    fn2=CrossCorr_dir//'z_'//Z//'/CrossCorr_'//RedShift//'_'//trim(MAX_Outer_char)//'_'//trim(MAX_Inner_char)//'_Run'//trim(Run)//'-'//trim(Version)//'-poisson-8000000.dat'
#else
    fn2=CrossCorr_dir//'z_'//Z//'/CrossCorr_'//RedShift//'_'//trim(MAX_Outer_char)//'_'//trim(MAX_Inner_char)//'_Run'//trim(Run)//'-'//trim(Version)//'-RDM.dat'
#endif
#endif

#ifdef LegendreTest
    fn2=CrossCorr_dir//'z_'//Z//'/CrossCorr_'//RedShift//'_'//trim(MAX_Outer_char)//'_'//trim(MAX_Inner_char)//'_Run'//trim(Run)//'-'//trim(Version)//'-Legendre.dat'
#endif

#ifdef power_law
    fn2=CrossCorr_dir//'z_'//Z//'/CrossCorr_'//RedShift//'_'//trim(MAX_Outer_char)//'_'//trim(MAX_Inner_char)//'_Run'//trim(Run)//'-'//trim(Version)//'-power_law_z_i=20.dat'
#endif

#ifdef GAUSS
    fn2=CrossCorr_dir//'z_'//Z//'/CrossCorr_'//RedShift//'_'//trim(MAX_Outer_char)//'_'//trim(MAX_Inner_char)//'_Run'//trim(Run)//'-'//trim(Version)//'-GAUSS.dat'
#ifdef TEST
    fn2=Test_dir//'CrossCorr_'//RedShift//'_'//trim(MAX_Outer_char)//'_'//trim(MAX_Inner_char)//'_Run'//trim(Run)//'-'//trim(Version)//'-GAUSS.dat'
#endif
#ifdef MNT_SCRATCH
    fn2='/mnt/node_scratch/jharno/data/CrossCorr_'//RedShift//'_'//trim(MAX_Outer_char)//'_'//trim(MAX_Inner_char)//'_Run'//trim(Run)//'-'//trim(Version)//'-GAUSS.dat'
#endif
#endif




#ifdef debug   
    write(*,*) 'TEST',MAX_Inner,MAX_Outer,2*MAX_Inner*MAX_Outer,2.0*MAX_Inner*MAX_Outer
    write(*,*) 'TEST',((MAX_Inner**2+MAX_Outer**2-(10-1.)**2)/(2*MAX_Inner*MAX_Outer))
#endif

    open(22,file=fn2)
    open(23, file='kr2weights_'//trim(MAX_Outer_char)//'_'//trim(MAX_Inner_char)//'.dat')
#ifdef kr2Bining
    
    do kr2=int(MAX_Outer-MAX_Inner)**2+1,int(MAX_Outer+MAX_Inner)**2+1 
       if(N_CrossCorrAve(kr2).gt.1e-13) write(22,*) (kr2-1), N_CrossCorrAve(kr2),CrossCorrAve(kr2),CrossCorrAve(kr2)/N_CrossCorrAve(kr2)   
       if(N_CrossCorrAve(kr2).gt.1e-13) write(23,*) kr2weight(kr2)
    enddo
    close(23)

#else

    do k=int(MAX_Outer-MAX_Inner)+1,int(MAX_Outer+MAX_Inner)+1!1,hc+1
       write(22,*) (k-1),acos((MAX_Inner**2+MAX_Outer**2-(k-1)**2)/(2*MAX_Inner*MAX_Outer)),N_CrossCorrAve(k),CrossCorrAve(k),CrossCorrAve(k)/N_CrossCorrAve(k) 
    enddo

#endif

    close(22)
    write(*,*) 'Wrote ', fn2

    


#ifdef FFT_TEST

    integer r
    character*100 fn8, fn9

    fn8=dir_work//'FFT_XSI.dat'
    fn9=dir_work//'FFT_HAT.dat'


    open(20,file=fn8)
    do r=2,nc! r=1 Asks for FFT_Hat(0) which is a NaN
       write(20,*) r-1, Xsi(r)   
    enddo
    close(20)

    open(21,file=fn9)
    do k=2,hc+1 
       write(21,*) k-1, HatAve(k), FFT_Hat(k),  GaussAve(k), DeltaAve(k) 
    enddo
    close(21)

    write(*,*) 'Wrote ', fn8
    write(*,*) 'Wrote ', fn9


#endif


    write(*,*) '**************'
    write(*,*) '*done writeCov*'
    write(*,*) '**************'

    return
  end subroutine writeCov



  !************************************
  !*** cic computes the density 'd' ***
  !************************************

  subroutine cic
    implicit none
    integer, parameter :: kpt=nc/nt

#ifdef HALO
    integer :: npt
#else
#ifdef subsample
    integer :: npt
#else
    integer, parameter :: npt=np/nt
#endif
#endif
    
    integer it,i,j,k,ip,OutBound
    real toe

#ifdef HALO
    npt=np/nt
    write(*,*) 'Performing cic on', npt, 'haloes'
#endif
#ifdef subsample
    npt=np/nt
    write(*,*) 'Performing cic on', npt, 'haloes'
#endif

    !! Construct chaining lists in parallel
    !$omp parallel do default(shared) private(it,ip,j,k)
    OutBound = 0
    do it=1,nt
       htoc(:,:,:,it)=0
       do ip=1+(it-1)*npt,min(np,it*npt)
          j=floor(xv(2,ip))+1
          k=floor(xv(3,ip))+1
          if((j > nc)) then
#ifdef debug_cic
             write (*,*) '#### PROBLEM!!! (j = floor(xv(2,ip))+1) =',j
             write (*,*) 'Enforcing Periodic BC Manually'
#endif
             j = j-nc
             OutBound = OutBound +1
             !pause
          endif
          if((k > nc)) then
#ifdef debug_cic
             write (*,*) '#### PROBLEM!!! (k = floor(xv(2,ip))+1) =',k
             write (*,*) 'Enforcing Periodic BC Manually'
#endif
             k = k-nc
             OutBound = OutBound +1
             !pause
          endif


          if (htoc(1,j,k,it) .eq. 0) then
             ll(ip)=0
             htoc(:,j,k,it)=ip
          else
             ll(htoc(2,j,k,it))=ip
             ll(ip)=0
             htoc(2,j,k,it)=ip
          endif
       enddo
    enddo
    !$omp end parallel do

    write(*,*) 'Enforced BC with ',OutBound, 'particles'

    !! Merge chaining lists
    !$omp parallel do default(shared) private(it,j,k,toe)
    do k=1,nc
       do j=1,nc
          hoc(j,k)=0
          do it=1,nt
             if (hoc(j,k) .eq. 0) then
                hoc(j,k)=htoc(1,j,k,it)
                toe=htoc(2,j,k,it)
             else
                if (htoc(1,j,k,it) .ne. 0) then
                   ll(toe)=htoc(1,j,k,it)
                   toe=htoc(2,j,k,it)
                endif
             endif
          enddo
       enddo
    enddo
    !$omp end parallel do


    !! Initialize density field
    !$omp parallel do default(shared) private(it,k)
    do it=1,nt
       do k=1+(it-1)*kpt,min(nc,it*kpt)
          d(:,:,k)=0
       enddo
    enddo
    !$omp end parallel do


    !! Add particle density to density field
!    do ko=1,4
    !$omp parallel do default(shared) private(j,k,ip)
       do k=1,nc
          do j=1,nc
             ip=hoc(j,k)
             call cicmass(ip)
          enddo
       enddo
!   enddo
    !$omp end parallel do


  write(*,*) '**********'
  write(*,*) '*done cic*'
  write(*,*) '**********'


    return
  end subroutine cic

  subroutine cicmass(ip)
    implicit none
    real, parameter :: ncr=nc

    integer ip

    integer i1,i2,j1,j2,k1,k2
    real x,y,z,dx1,dx2,dy1,dy2,dz1,dz2

    do
       if (ip .eq. 0) exit

       x=modulo(xv(1,ip)-0.5+ncr,ncr)
       y=modulo(xv(2,ip)-0.5+ncr,ncr)
       z=modulo(xv(3,ip)-0.5+ncr,ncr)

       i1=floor(x)+1
       i2=mod(i1,nc)+1
       dx1=i1-x
       dx2=1-dx1
       j1=floor(y)+1
       j2=mod(j1,nc)+1
       dy1=j1-y
       dy2=1-dy1
       k1=floor(z)+1
       k2=mod(k1,nc)+1
       dz1=k1-z
       dz2=1-dz1

       dz1=mp*dz1
       dz2=mp*dz2
       d(i1,j1,k1)=d(i1,j1,k1)+dx1*dy1*dz1
       d(i2,j1,k1)=d(i2,j1,k1)+dx2*dy1*dz1
       d(i1,j2,k1)=d(i1,j2,k1)+dx1*dy2*dz1
       d(i2,j2,k1)=d(i2,j2,k1)+dx2*dy2*dz1
       d(i1,j1,k2)=d(i1,j1,k2)+dx1*dy1*dz2
       d(i2,j1,k2)=d(i2,j1,k2)+dx2*dy1*dz2
       d(i1,j2,k2)=d(i1,j2,k2)+dx1*dy2*dz2
       d(i2,j2,k2)=d(i2,j2,k2)+dx2*dy2*dz2
          
       ip=ll(ip)
    enddo
   
    return
  end subroutine cicmass


!!--------------------------------------------------------------!!

  subroutine powerspectrum
    implicit none

    !real pst(2,nc,nt)
    integer i,j,k,n_counts, kz
    
#ifdef GAUSS
    real, dimension(nc) :: kgauss
    character(len=MSL) :: psGaussFile
#endif 

#ifdef RANDOM
    integer k_seed, clock
    integer, dimension(2) :: old, seed
#endif

!************

#ifdef LegendreTest

    !$omp parallel do default(shared) private(i,j,k,kz)
    do k=1,nc

       if(k .lt. hc+2) then
          kz = k-1 ! ->[0,hc]
       else
          kz = k-1-nc ! ->[-hc+1,-1]
       endif      
       do j=1,nc          
          do i=1,hc+1
             psV_COMPLEX(i,j,k) = kz**2
          enddo
       enddo
       !write(*,*) 'psV_COMPLEX(4,4,k) =', psV_COMPLEX(4,4,k) 
    enddo
    !$omp end parallel do

    write(*,*) '*******************************'
    write(*,*) 'Overwriting with Legendre Field'
    write(*,*) '*******************************'

    write(*,*) 'Check that these numbers agree with P(k) = k_z^2:'

    write(*,*) 'psV(1,1,1)=', psV_COMPLEX(1,1,1)
    write(*,*) 'psV(2,2,2)=', psV_COMPLEX(2,2,2)
    write(*,*) 'psV(128,128,128)=', psV_COMPLEX(128,128,128)

    !call AngleAverageKSpace(d,ps,err,nc)
    call Get_KAverage(d,ps,nc)

  return

#endif

!************
#ifdef RANDOM   
    write(*,*) '******************************'
    write(*,*) 'Generating Flat Random Density'
    write(*,*) '******************************'
    
    call random_seed
    call random_seed(size=k_seed)
    !write(*,*) 'size k_seed = ', k_seed    
    call system_clock(count=clock)
    seed = clock + 37 * (/ (i - 1, i = 1, 2) /)
    call random_seed(put = seed)
    call random_seed(get=old(1:k_seed))
    !write(*,*) ' Old starting value = ', old
    !write(*,*) ' Seed value =', seed
    call random_number(d) 
    !d = d-0.5

#ifdef poisson

    n_counts = 0
    write(*,*) 'Converting to poisson density'
 
   !$omp parallel do default(shared) private(i,j,k)
    do i = 1,nc
       do j = 1,nc
          do k = 1,nc
             !if(d(i,j,k) > 0.49993) then
             if(d(i,j,k) > 1.0 -(8000000.0/nc/nc/nc)) then
                 d(i,j,k) = 1.0
                !write(*,*) 'Got one'
                n_counts = n_counts+1
             else
                d(i,j,k) = 0.0
             endif
          enddo
       enddo
    enddo
    !$omp end parallel do

    write(*,*) 'Got ', n_counts, ' counts'

    write(*,*) 'Converting to densities contrast'
    d = d*nc*nc*nc/n_counts - 1.0


    !stop
#else
    !random field btw 0 and 1 becomes overdensities btw -0.5 and 0.5
    d = d-0.5
#endif
#else
    
    !! Convert density to overdensity
    d = d-1        
#endif

#ifdef debug_density
    write(*,*) 'd before fft', maxval(d),minval(d) 
#endif 

!************

    !call ps3_r2c(d,ps,nc)
#ifdef NGP_PS
    call ps3_r2c_err_NGP(d,ps,err,nc) 
#else    
    call ps3_r2c_err(d,ps,err,nc)
#endif

#ifdef debug_density
    write(*,*) 'd after fft', maxval(d),minval(d) 
    write(*,*) 'Delta^2(k) after fft', maxval(ps),minval(ps) 
#endif 

#ifdef GAUSS

    write(*,*) '********************************'
    write(*,*) 'Generating Gaussian Random Field'
    write(*,*) '********************************'

    do i = 1,nc
       kgauss(i) = i*2*pi/lbox
    enddo

#ifdef debug_gauss
    write(*,*)'Filled kgauss'
#endif

    call GaussRandomField_3d_r2c(d, lbox, nc, kgauss, ps, nc)

    call ps3_r2c_err_NGP(d,ps_gauss,err,nc)

    !Get Average PowerSpectrum from 43 simulations from File
    !psGaussFile = dir_work//'/PowerSpectrum/z_'//Z//'/psGaussAve_Z'//Z//'.txt'
    !
    !write(*,*) 'writing', psGaussFile
    !open(unit=25,file=psGaussFile)
    !***** To write/read the power spectrum to/from a file
    !do i=1,hc 
    !   write(25,*) kgauss(i), ps_gauss(i)
    !enddo
    !close(25)


#endif



    write(*,*) '********************'
    write(*,*) '*done powerspectrum*'
    write(*,*) '********************'



    return
  end subroutine powerspectrum

!!--------------------------------------------------------------!!

!!--------------------------------------------------------------!!

  subroutine CrossCorrelationOfShells
    implicit none


    integer i,j,k, kk, kx, ky, kz, N_offshell, N_onshell, N_onshell_kx
    !real MAX_Inner, Max_Outer, MIN_Inner, MIN_Outer, kr, Thickness_Inner, Thickness_Outer, sigma, r0, k0
    real kr, MIN_Inner, MIN_Outer, sigma, r0, k0, sigmaP_shell !,MAX_Inner, MAX_Outer, Thickness_Inner, Thickness_Outer
    complex(8) psShellCumul_COMPLEX, psShellMean_COMPLEX, N_psShellCumul_COMPLEX, N_psShellMean_COMPLEX

    !real, dimension(nc+2,nc,nc) ::  psShell_Inner, psShell_Outer, N_psShell_Inner, N_psShell_Outer
    !complex, dimension(nc/2 + 1,nc,nc) ::  psShell_Inner_COMPLEX, psShell_Outer_COMPLEX, N_psShell_Inner_COMPLEX,N_psShell_Outer_COMPLEX
    !complex(8), dimension(nc/2 + 1,nc,nc) ::  psShell_COMPLEX, N_psShell_COMPLEX

    !equivalence (psShell_Inner,psShell_Inner_COMPLEX)
    !equivalence (psShell_Outer,psShell_Outer_COMPLEX)
    !equivalence (N_psShell_Inner,N_psShell_Inner_COMPLEX)
    !equivalence (N_psShell_Outer,N_psShell_Outer_COMPLEX)

    real, dimension(hc) :: KFromFile
!    real, dimension(hc) :: PSFromFile
    real, dimension(hc) :: psAveFromShells
    real, dimension(hc/2) :: Num_OnShell,psShellCumul

!    character(len=MSL) :: psFile

   
    !***** Write/read the power spectra to/from a file

    !Get Average PowerSpectrum from 299 simulations from File
    !Either interpolated or straight from shell
    !*** Interpolated ***
    !psFile = dir_work//'/PowerSpectrum/z_'//Z//'/psAve_Z'//Z//'.txt'
    !*** Straight from shells ***
    !psFile = dir_work//'/PowerSpectrum/z_'//Z//'/ps128Shell_N50.dat'
     !psFile = '../data/rho_Small/meanP_Z0.txt'
!#ifdef GAUSS
!    psFile = dir_work//'/PowerSpectrum/z_'//Z//'/ps128Shell_N50-GAUSS.dat'
!#endif
!    if(z3dps.eq.0)then
!       !psFile = '/cita/h/home-1/jharno/Projects/bao_fisher/results/meanP_Z0-NGP.txt'
!       psFile = '/home/jharno/AngularCovariance/meanP_Z0-NGP.txt'
!    else
!       write(*,*) 'wrong PSFromFile!'
!    endif
! Actually, the cells on the shell is now done exactly the same as in NGP P(k),
! so we can precompute the average P(k) over all simulations and substract it directly.


!#ifndef GetPowerShells
!    write(*,*) 'opening', psFile
!    open(unit=24,file=psFile)
!    do i=1,hc/2 
!       read(24,*) PSFromFile(i)
!#ifdef debug
!       write(*,*) i, PSFromFile(i)
!#endif 
!    enddo
!    close(24)
!#endif

    ! Initialize variables
    N_offshell = 0 
    N_onshell = 0
    N_onshell_kx = 0
    N_psShellCumul_COMPLEX = 0.
    psShellCumul_COMPLEX = 0.
    sigmaP_shell = 0.

    !density d already in fourier space after "power-spectrum" subroutine
    !and equivalenced to psV

    !psV_COMPLEX = lbox**3*abs(d_COMPLEX/nc**3)**2
!    psV_COMPLEX = lbox**3*abs(psV_COMPLEX/nc**3)**2
    
#ifdef debug_AngCov
    write(*,*) 'Got Vectorial Power Spectrum, max Re = ',maxval(real(psV_COMPLEX)), 'min Re= ', minval(real(psV_COMPLEX))
    write(*,*) 'Got Vectorial Power Spectrum, max Im = ',maxval(aimag(psV_COMPLEX)), 'min Im= ', minval(aimag(psV_COMPLEX))
#endif

#ifdef LegendreTest


    !$omp parallel do default(shared) private(i,j,k,kz)
    do k=1,nc

       if(k .lt. hc+2) then
          kz = k-1 ! ->[0,hc]
       else
          kz = k-1-nc ! ->[-hc+1,-1]
       endif
        
       do j=1,nc
          
          do i=1,hc+1

             psV_COMPLEX(i,j,k) = kz**2

          enddo
       enddo

       !write(*,*) 'psV_COMPLEX(4,4,k) =', psV_COMPLEX(4,4,k) 

    enddo
    !$omp end parallel do

    write(*,*) '*******************************'
    write(*,*) 'Overwriting with Legendre Field'
    write(*,*) '*******************************'

    write(*,*) 'Check that these numbers agree with P(k) = k_z^2:'

    write(*,*) 'psV(1,1,1)=', psV_COMPLEX(1,1,1)
    write(*,*) 'psV(2,2,2)=', psV_COMPLEX(2,2,2)
    write(*,*) 'psV(128,128,128)=', psV_COMPLEX(128,128,128)

    

#endif



    !***************

    ! Define a shell, in fraction of hc (radius of biggest sphere in the k-box)
    ! and set the power spectrum to 0 everywhere else. 
    ! The Shell must be smaller than hc


    write(*,*) '******************************'
    write(*,*) 'Starting Loop over Inner Shell'
    write(*,*) '******************************'



    MIN_Inner = MAX_Inner -  Thickness_Inner


    ! Loop over all k-cells, get magnitude kr  

    !$omp  parallel do default(shared) private(i,j,k,kx,ky,kz,kr)
    do k=1,nc
       if(k .lt. hc+2) then
          kz = k-1 ! ->[0,hc]
       else
          kz = k-1-nc ! ->[-hc+1,-1]
       endif
       do j=1,nc
          if(j .lt. hc+2) then
             ky = j-1 !-> [0,hc]
          else
             ky = j-1-nc ! -> [-hc+1,-1]
          endif

          !do i=1,nc+2,2
          do i=1,hc+1 ! -> [1,hc+1]            
             kx = i-1 ! -> [0,hc]
             !kx = (i-1)/2

             !Spherical Coordinates
             kr= ( (kx)**2 + (ky)**2 + (kz)**2 )**(0.5)


             ! I should treat the x=0 slab separately, to match the NGP power spectrum
             ! But I will leave it to zero for now
             if(kx.eq.0 .and. ky <=0 .and. kz <=0) cycle;
             if(kx.eq.0 .and. ky >0 .and. kz <0) cycle;
             if(kr.eq.0) cycle;
           

             if (kr <= MIN_Inner .or. kr > MAX_Inner) then
                N_offshell = N_offshell +1
                psShell_Inner_COMPLEX(i,j,k) = 0 ! Set power spectrum to zero
                N_psShell_Inner_COMPLEX(i,j,k) = 0
             else
                N_onshell = N_onshell +1 ! Count # cells on shell
                psShell_Inner_COMPLEX(i,j,k) = psV_COMPLEX(i,j,k) - PSFromFile(MAX_Inner)! Assign power spectrum
                psShellCumul_COMPLEX = psShellCumul_COMPLEX + psShell_Inner_COMPLEX(i,j,k)
                N_psShell_Inner_COMPLEX(i,j,k) = 1
                N_psShellCumul_COMPLEX = N_psShellCumul_COMPLEX + N_psShell_Inner_COMPLEX(i,j,k)
                sigmaP_shell =  sigmaP_shell + (psV_COMPLEX(i,j,k) - PSFromFile(MAX_Inner))**2

                if(kx == 0) N_onshell_kx = N_onshell_kx+1
                !if(kr <= 1)write(*,*)'k=',kx,ky,kz,psV_COMPLEX(i,j,k),PSFromFile(MAX_Inner)


             endif


!**************************************************
! To get the power of each shells
! 1- Loop over cells in the shells (done above here)
! 2- Count and Accumulate the power for them
! 3- Divide ps/N to get mean, write to file...
#ifdef GetPowerShells
             do kk = 1,nc/4
                if ( (kk-1 < kr) .and. (kr <= kk) ) then
                   Num_OnShell(kk) = Num_OnShell(kk)+1
                   psShellCumul(kk) = psShellCumul(kk) + real(psV_COMPLEX(i,j,k))
                endif
             enddo
#endif
!***************************************************              
             
 
          enddo
       enddo
    enddo
    !$omp  end parallel do



!***************************************************
#ifdef GetPowerShells
    psFile = dir_work//'/PowerSpectrum/z_'//Z//'/AllShells_Run'//trim(Run)//'-'//trim(Version)//'.dat'
#ifdef GAUSS
    psFile = dir_work//'/PowerSpectrum/z_'//Z//'/AllShells_Run'//trim(Run)//'-'//trim(Version)//'-GAUSS.dat'
#endif

    !write(*,*) 'Writing', psFile ,'then closing'
    open(unit=30,file=psFile)

    do kk = 1,nc/4
       !write(*,*)'Cumul Shell Power =',psShellCumul(kk)
       !write(*,*)'Number of cells on Shell =',Num_OnShell(kk)
       !pause
       !write(*,*) Num_OnShell(kk), psShellCumul(kk), psShellCumul(kk)/Num_OnShell(kk)
       write(30,*) psShellCumul(kk)/Num_OnShell(kk)
    enddo

    close(30)

    write(*,*)'Wrote', psFile
    !return
    stop

#endif
!*************************************************


    ! **********************************

    if(N_onshell .eq. 0)then
       write(*,*)'N_onshell = 0, problem coming...'
    endif

    ! Compute Mean powerspectrum on the Shell
    psShellMean_COMPLEX = psShellCumul_COMPLEX/N_onshell
    N_psShellMean_COMPLEX = N_psShellCumul_COMPLEX/N_onshell
    sigmaP_shell = sqrt(sigmaP_shell/(N_onshell-1))

    ! **********************************


    write(*,*) 'MAX = ', MAX_Inner, 'MIN = ',MIN_Inner

#ifdef debug_shell

    write(*,*) 'N k-cells offshell : ', N_offshell
    write(*,*) 'N k-cells onshell : ', N_onshell
    write(*,*) 'N k-cells onshell_kx : ', N_onshell_kx
    write(*,*) '*******************************'
    !write(*,*) 'Total Power on Shell = ', psShellCumul_COMPLEX
    !write(*,*) 'Mean Power on Shell = ', psShellMean_COMPLEX
    write(*,*) 'Mean Power on Shell = ', psShellMean_COMPLEX + PSFromFile(MAX_Inner)
    write(*,*) 'After substraction  = ', psShellMean_COMPLEX
    write(*,*) 'std on shell        = ', sigmaP_shell
    write(*,*) 'Mean Norm on Shell  = ', N_psShellMean_COMPLEX ! Make sure this equals to 1!!!!
    write(*,*) '*******************************'
    write(*,*) '<P>(k-shell) = ',PSFromFile(MAX_Inner)
    write(*,*) '*******************************'

#endif



    ! **********************************

    !Substract mean from power spectrum on shell:

    !write(*,*) 'sum psShell (before mean substraction) =' , sum(psShell_Inner_COMPLEX)


    ! Loop again over all k-cells, get magnitude kr  

    !$omp  parallel do default(shared) private(i,j,k,kx,ky,kz,kr)
    !do k=1,nc
    !   if(k .lt. hc+2) then
    !      kz = k-1 ! ->[0,hc]
    !   else
    !      kz = k-1-nc ! ->[-hc+1,-1]
    !   endif
    !   do j=1,nc
    !      if(j .lt. hc+2) then
    !         ky = j-1 !-> [0,hc]
    !      else
    !         ky = j-1-nc ! -> [-hc+1,-1]
    !      endif
    !
    !      !do kx=1,hc+1 ! -> [1,hc+1]
    !      do i=1,hc+1 ! 
    !         
    !         kx = i-1 ! -> [0,hc]
    !         
    !         !kr=( (kx)**2 + (ky-hc)**2 + (kz-hc)**2)**(0.5)
    !         kr=( (kx)**2 + (ky)**2 + (kz)**2 )**(0.5)
    !         if(kx.eq.0 .and. ky <=0 .and. kz <=0) cycle;
    !         if(kx.eq.0 .and. ky >0 .and. kz <0) cycle;
    !         if(kr.eq.0) cycle;
    !
    !         if (kr <= MIN_Inner .or. kr > MAX_Inner) then
    !            psShell_Inner_COMPLEX(i,j,k) = psShell_Inner_COMPLEX(i,j,k)
    !            N_psShell_Inner_COMPLEX(i,j,k) = N_psShell_Inner_COMPLEX(i,j,k)
    !         else
    !            !Un-comment the next line to substract exact mean from the current shell
    !            psShell_Inner_COMPLEX(i,j,k) = psShell_Inner_COMPLEX(i,j,k) - psShellMean_COMPLEX                
    !
    !            !Un-comment the next line to to substract an interpolated average power spectrum over 50 simulations
    !            !psShell_Inner_COMPLEX(i,j,k) = psShell_Inner_COMPLEX(i,j,k) - 0.5*(PSFromFile(MAX_Inner)+PSFromFile(MAX_Inner-1))         
    !
    !            !Un-comment the next line to to substract the mean of the shell power spectrum over 50 simulations
    !            !psShell_Inner_COMPLEX(i,j,k) = psShell_Inner_COMPLEX(i,j,k) - PSFromFile(MAX_Inner)        
    !
    !            !Un-comment the next line to verify code:  Norm should be exactly == 0
    !            !N_psShell_Inner_COMPLEX(i,j,k) = N_psShell_Inner_COMPLEX(i,j,k) - N_psShellMean_COMPLEX 
    !
    !         endif
    !
    !      enddo
    !   enddo
    !enddo
    !$omp end parallel do


!#ifdef debug
!    write(*,*) 'sum psShell (after mean substraction) =' , sum(psShell_Inner_COMPLEX)
!#endif

    !psShell_Inner_COMPLEX = psShell_Inner_COMPLEX
    !N_psShell_Inner_COMPLEX = N_psShell_Inner_COMPLEX

    !****************************


    ! Re-Initialize variables
    N_offshell = 0 
    N_onshell = 0
    N_onshell_kx = 0
    N_psShellCumul_COMPLEX = 0.
    psShellCumul_COMPLEX = 0.



    write(*,*) '******************************'
    write(*,*) 'Starting Loop over Outer Shell'
    write(*,*) '******************************'



    MIN_Outer = MAX_Outer -  Thickness_Outer



    ! Loop over all k-cells, get magnitude kr  

    !$omp  parallel do default(shared) private(i,j,k,kx,ky,kz,kr)
    do k=1,nc
       if(k .lt. hc+2) then
          kz = k-1 ! ->[0,hc]
       else
          kz = k-1-nc ! ->[-hc+1,-1]
       endif
       do j=1,nc
          if(j .lt. hc+2) then
             ky = j-1 !-> [0,hc]
          else
             ky = j-1-nc ! -> [-hc+1,-1]
          endif

          !do i=1,nc+2,2! -> [1,hc+1]
          do i=1,hc+1! -> [1,hc+1]
             
             kx = i-1 ! -> [0,hc]
             !kx = (i-1)/2 ! -> [0,hc]

             !Spherical Coordinates
             kr= ( (kx)**2 + (ky)**2 + (kz)**2 )**(0.5)
             if(kx.eq.0 .and. ky <=0 .and. kz <=0) cycle;
             if(kx.eq.0 .and. ky >0 .and. kz <0) cycle;
             if(kr.eq.0) cycle;

             

             if (kr <= MIN_Outer .or. kr > MAX_Outer) then
                N_offshell = N_offshell +1
                psShell_Outer_COMPLEX(i,j,k) = 0 ! Set power spectrum to zero
                N_psShell_Outer_COMPLEX(i,j,k) = 0
             else
                N_onshell = N_onshell +1 ! Count # cells on shell
                psShell_Outer_COMPLEX(i,j,k) = psV_COMPLEX(i,j,k) -PSFromFile(MAX_Outer)! Assign power spectrum
                psShellCumul_COMPLEX = psShellCumul_COMPLEX + psShell_Outer_COMPLEX(i,j,k)
                N_psShell_Outer_COMPLEX(i,j,k) = 1
                N_psShellCumul_COMPLEX = N_psShellCumul_COMPLEX + N_psShell_Outer_COMPLEX(i,j,k)

                if(kx == 0) N_onshell_kx = N_onshell_kx+1

             endif

 
          enddo
       enddo
    enddo
    !$omp end parallel do

    ! **********************************

    ! Compute Mean powerspectrum on the Shell
    psShellMean_COMPLEX = psShellCumul_COMPLEX/N_onshell
    N_psShellMean_COMPLEX = N_psShellCumul_COMPLEX/N_onshell

    ! **********************************


    write(*,*) 'MAX = ', MAX_Outer, 'MIN = ',MIN_Outer

#ifdef debug_shell

    write(*,*) 'N k-cells offshell : ', N_offshell
    write(*,*) 'N k-cells onshell : ', N_onshell
    write(*,*) 'N k-cells onshell_kx : ', N_onshell_kx
    write(*,*) '*******************************'
    !write(*,*) 'Total Power on Shell = ', psShellCumul_COMPLEX
    !write(*,*) 'Mean Power on Shell = ', psShellMean_COMPLEX
    write(*,*) 'Mean Power on Shell = ', psShellMean_COMPLEX + PSFromFile(MAX_Outer)
    write(*,*) 'After substraction = ', psShellMean_COMPLEX
    write(*,*) 'Mean Norm on Shell = ', N_psShellMean_COMPLEX ! Make sure this equals to 1!!!!
    write(*,*) '*******************************'
    write(*,*) '<P>(k-shell) = ',PSFromFile(MAX_Outer)
    write(*,*) '*******************************'

#endif

    ! **********************************

    !Substract mean from power spectrum on shell:

!    write(*,*) 'sum psShell (before mean substraction) =' , sum(psShell_Outer_COMPLEX)


    ! Loop again over all k-cells, get magnitude kr  

    !$omp  parallel do default(shared) private(i,j,k,kx,ky,kz,kr)
!    do k=1,nc
!       if(k .lt. hc+2) then
!          kz = k-1 ! ->[0,hc]
!       else
!          kz = k-1-nc ! ->[-hc+1,-1]
!       endif
!       do j=1,nc
!          if(j .lt. hc+2) then
!             ky = j-1 !-> [0,hc]
!          else
!             ky = j-1-nc ! -> [-hc+1,-1]
!          endif
!
!          !do kx=1,hc+1 ! -> [1,hc+1]
!          do i=1,hc+1
!             
!             kx = i-1 ! -> [0,hc]
!             
!             !kr=( (kx)**2 + (ky-hc)**2 + (kz-hc)**2)**(0.5)
!             kr=( (kx)**2 + (ky)**2 + (kz)**2 )**(0.5)
!             if(kx.eq.0 .and. ky <=0 .and. kz <=0) cycle;
!             if(kx.eq.0 .and. ky >0 .and. kz <0) cycle;
!             if(kr.eq.0) cycle;
!
!             if (kr <= MIN_Outer .or. kr > MAX_Outer) then
!                psShell_Outer_COMPLEX(i,j,k) = psShell_Outer_COMPLEX(i,j,k)
!                N_psShell_Outer_COMPLEX(i,j,k) = N_psShell_Outer_COMPLEX(i,j,k)
!             else
!                !Un-comment the next line to substract exact mean from the current shell
!                psShell_Outer_COMPLEX(i,j,k) = psShell_Outer_COMPLEX(i,j,k) - psShellMean_COMPLEX                
!
!                !Un-comment the next line to to substract an interpolated average power spectrum over 50 simulations
!                !psShell_Outer_COMPLEX(i,j,k) = psShell_Outer_COMPLEX(i,j,k) - 0.5*(PSFromFile(MAX_Outer)+PSFromFile(MAX_Outer-1))         
!
!                !Un-comment the next line to to substract the mean of the shell power spectrum over 50 simulations
!                !psShell_Outer_COMPLEX(i,j,k) = psShell_Outer_COMPLEX(i,j,k) - PSFromFile(MAX_Outer)         
!
!                !Un-comment the next line to verify code:  Norm should be exactly == 0
!                !N_psShell_Outer_COMPLEX(i,j,k) = N_psShell_Outer_COMPLEX(i,j,k) - N_psShellMean_COMPLEX 
!
!             endif
!
!          enddo
!       enddo
!    enddo
    !$omp end parallel do

!#ifdef debug
!    write(*,*) 'sum psShell (after mean substraction) =' , sum(psShell_Outer_COMPLEX)
!#endif



    !******************************************************
    write(*,*) '*** Cross-Correlation of Normalization ***'
    !******************************************************

    call GetCrossCorr(N_psShell_Inner, N_psShell_Outer, N_CrossCorr, nc)

    write(*,*) 'Sum of Norm = ', sum(N_CrossCorr_COMPLEX)

#ifdef debug_AngCov
    write(*,*) 'Vectorial Cross-Correlation Re: ', maxval(real(N_CrossCorr_COMPLEX)),' > N_CrossCorr > ',minval(real(N_CrossCorr_COMPLEX))
    write(*,*) 'Vectorial Cross-Correlation Im: ', maxval(aimag(N_CrossCorr_COMPLEX)),' > N_CrossCorr > ',minval(aimag(N_CrossCorr_COMPLEX))
#endif


#ifdef kr2Bining
    call kr2BinAverage(N_CrossCorr,N_CrossCorrAve,nc,1.,1.)
#else
    call Get_KAverage(N_CrossCorr,N_CrossCorrAve,nc,1.,1.)
#endif

#ifdef debug_AngCov
    write(*,*) 'Averaged Cross-Correlation', maxval(N_CrossCorrAve),' > CrossCorrAve > ',minval(N_CrossCorrAve) 
#endif



    !**********************************************
    write(*,*) '*** Cross-Correlation of Shell ***'
    !**********************************************

#ifdef debug_AngCov
    write(*,*) 'sum psShell_Inner (after mean substraction) =' , sum(psShell_Inner_COMPLEX)
    write(*,*) 'sum psShell_Outer (after mean substraction) =' , sum(psShell_Outer_COMPLEX)
#endif

    call GetCrossCorr(psShell_Inner, psShell_Outer, CrossCorr, nc)

    write(*,*) 'Sum of Cross Correlation = ', sum(CrossCorr_COMPLEX)

#ifdef debug_AngCov
    write(*,*) 'Vectorial Cross-Correlation', maxval(CrossCorr),' > CrossCorr > ',minval(CrossCorr)
    write(*,*) 'Vectorial Cross-Correlation Re: ', maxval(real(CrossCorr_COMPLEX)),' > CrossCorr > ',minval(real(CrossCorr_COMPLEX))
    write(*,*) 'Vectorial Cross-Correlation Im: ', maxval(aimag(CrossCorr_COMPLEX)),' > CrossCorr > ',minval(aimag(CrossCorr_COMPLEX))

    write(*,*) 'TEST : SUM over CrossCorrelation Re = ', sum(real(CrossCorr_COMPLEX))
    write(*,*) 'TEST : SUM over CrossCorrelation Im = ', sum(aimag(CrossCorr_COMPLEX))
#endif


!! The following prints the Values of CrossCorr above a certain threshold, along with the coordinates
!! Typically, these appear 
!
!    do k=1,nc
!       if(k .lt. hc+2) then
!          kz = k-1 ! ->[0,hc]
!       else
!          kz = k-1-nc ! ->[-hc+1,-1]
!       endif
!       do j=1,nc
!          if(j .lt. hc+2) then
!             ky = j-1 !-> [0,hc]
!          else
!             ky = j-1-nc ! -> [-hc+1,-1]
!          endif
!          
!          !do kx=1,hc+1 ! -> [1,hc+1]
!          do i=1,hc+1
!             
!             kx = i-1 ! -> [0,hc]
!
!             kr=( (kx)**2 + (ky)**2 + (kz)**2 )**(0.5)
!             
!             if(real(CrossCorr_COMPLEX(i,j,k)) > 3) then
!                write(*,*) 'k = (',kx,ky,kz,'), kr = ',kr  
!                write(*,*) 'CrossCorr(k) = ', CrossCorr_COMPLEX(i,j,k)
!                pause
!             endif
!
!          enddo
!       enddo
!    enddo


#ifdef kr2Bining
    call kr2BinAverage(CrossCorr,CrossCorrAve,nc,1.,1.)
#else
    call Get_KAverage(CrossCorr,CrossCorrAve,nc,1.,1.)
#endif




    write(*,*) '***********************'
    write(*,*) '*done CrossCorrelation*'
    write(*,*) '***********************'



    return
  end subroutine CrossCorrelationOfShells

!*********************************

!! 1- FFT-Inverse a Power Spectrum to.get the Correlation Function
!! 2- AbsoluteValue-Square the Correlation Function
!! 3- FFT to get a Cross-Correlation(delta_k)

subroutine GetCrossCorr(ps1,ps2,Cross,nc)
!subroutine GetCrossCorr(ps1_COMPLEX,ps2_COMPLEX,Cross_COMPLEX,nc)
  implicit none

  integer nc
  real, dimension(nc+2,nc,nc) :: ps1, ps2, Cross  

#ifdef debug_xi
  write(*,*) 'P1:sum over density in Fourier Space =' , sum(ps1),'nc = ', nc
#endif

  call sfft3_r2c(ps1,nc,-1) 

#ifdef debug_xi
  write(*,*) 'Xi_1(0,0,0) = ', ps1(1,1,1)   ! NOTE: Xi(0,0,0) ~ 2*sum/nc**3 (this double counts the x=0 slab)
  !write(*,*) 'Xi_1(1,0,0) = ', ps1(2,1,1)   ! NOTE: Xi(0,0,0) ~ 2*sum/nc**3 (this double counts the x=0 slab)
  write(*,*) 'P2:sum over density in Fourier Space =' , sum(ps2),'nc = ', nc
#endif

  call sfft3_r2c(ps2,nc,-1) 
#ifdef debug_xi
  write(*,*) 'Xi_2(0,0,0) = ', ps2(1,1,1)   ! NOTE: Xi(0,0,0) ~ 2*sum/nc**3 (this double counts the x=0 slab)
  !write(*,*) 'Xi_2(1,0,0) = ', ps2(2,1,1)   ! NOTE: Xi(0,0,0) ~ 2*sum/nc**3 (this double counts the x=0 slab)
#endif


  !! Now ps1 and ps2 are a "real" function, going from 1 to nc 
  !! in all 3 dimensions. 

  Cross = ps1*ps2


#ifdef debug_xi

  write(*,*) '***************'
  write(*,*) 'Cross before FFT : '
  write(*,*) '***************'
  write(*,*) 'Xi_1*Xi_2(0,0,0) = ', Cross(1,1,1)
  write(*,*) 'Xi_1*Xi_2(1,0,0) = ', Cross(2,1,1)
  write(*,*) 'Xi_1*Xi_2(0,1,0) = ', Cross(1,2,1)
  write(*,*) 'Xi_1*Xi_2(0,0,1) = ', Cross(1,1,2)
  write(*,*) 'Xi_1*Xi_2(1,1,0) = ', Cross(2,2,1)
  write(*,*) 'Xi_1*Xi_2(1,0,1) = ', Cross(2,1,2)
  write(*,*) 'Xi_1*Xi_2(0,1,1) = ', Cross(1,2,2)
  write(*,*) 'Xi_1*Xi_2(1,1,1) = ', Cross(2,2,2)
  write(*,*) 'Xi_1*Xi_2(2,0,0) = ', Cross(3,1,1)
  write(*,*) 'Xi_1*Xi_2(0,2,0) = ', Cross(1,3,1)
  write(*,*) 'Xi_1*Xi_2(0,0,2) = ', Cross(1,1,3)
  write(*,*) 'Xi_1*Xi_2(2,1,0) = ', Cross(3,2,1)
  write(*,*) 'Xi_1*Xi_2(2,0,1) = ', Cross(3,1,2)
  write(*,*) 'Xi_1*Xi_2(1,2,0) = ', Cross(2,3,1)
  write(*,*) 'Xi_1*Xi_2(0,2,1) = ', Cross(1,3,2)
  write(*,*) 'Xi_1*Xi_2(0,1,2) = ', Cross(1,2,3)
  write(*,*) 'Xi_1*Xi_2(1,0,2) = ', Cross(2,1,3)
  write(*,*) '***************'
  write(*,*) maxval(Cross),' > Cross > ' ,minval(Cross)
  write(*,*) '***************'

#endif
 
  call sfft3_r2c(Cross,nc,1) 

  !Now Cross is in k-space. This subroutine defined Cross as a real, 
  !but even x-array elements represent imaginary parts


#ifdef debug_xi


  !Recall that in k-space, even-kx represent imaginary parts...
  write(*,*) '***************'
  write(*,*) 'Cross after FFT : '
  write(*,*) '***************'
  write(*,*) 'rho(0,0,0) = ', Cross(1,1,1)
  write(*,*) 'rho(1,0,0) = ', Cross(3,1,1)
  write(*,*) 'rho(0,1,0) = ', Cross(1,2,1)
  write(*,*) 'rho(0,0,1) = ', Cross(1,1,2)
  write(*,*) 'rho(1,1,0) = ', Cross(3,2,1)
  write(*,*) 'rho(1,0,1) = ', Cross(3,1,2)
  write(*,*) 'rho(0,1,1) = ', Cross(1,2,2)
  write(*,*) 'rho(1,1,1) = ', Cross(3,2,2)
  write(*,*) 'rho(2,0,0) = ', Cross(5,1,1)
  write(*,*) 'rho(0,2,0) = ', Cross(1,3,1)
  write(*,*) 'rho(0,0,2) = ', Cross(1,1,3)
  write(*,*) 'rho(2,1,0) = ', Cross(5,2,1)
  write(*,*) 'rho(2,0,1) = ', Cross(5,1,2)
  write(*,*) 'rho(1,2,0) = ', Cross(3,3,1)
  write(*,*) 'rho(0,2,1) = ', Cross(1,3,2)
  write(*,*) 'rho(0,1,2) = ', Cross(1,2,3)
  write(*,*) 'rho(1,0,2) = ', Cross(3,1,3)
  write(*,*) '***************'
  write(*,*) maxval(Cross),' > Cross > ' ,minval(Cross)
  write(*,*) '***************'

#endif


  return
end subroutine GetCrossCorr

!*******************************


!======================================================
 Subroutine Get_KAverage(map,Ave,n) !,normp,normk)
!======================================================
! Calculates the angular average in Fourier space


   Integer                            :: n
   !Real                               :: normp,normk
   Real, Dimension(1:n+2,1:n,1:n)  :: map
   Character(Len=120)              :: file

   !Real*4, Dimension(1:n+2,1:n,1:n):: fft
   Real, Dimension(3,1:n)        :: pst
   Real, Dimension(n)        :: Ave

   Real  :: w1,w2,kz,kx,ky,kr,pow
   Integer :: i,j,k,hn,k1,k2

   If (mod(n,2) .NE. 0) Then
      Write(*,*) ' n must be even'
      Stop
   End If


   !! == Dump power spectra
   !Open(unit=40,file=file,status='replace')
   hn  = n/2
   pst = 0.0

   !$omp  parallel do default(shared) private(i,j,k,kx,ky,kz,kr,k1,k2,pow,w1,w2)
   Do k = 1,n
      If (k .Lt. hn+2) Then
         kz = k-1
      Else
         kz = k-1-n
      Endif
      Do j = 1,n
         If (j .Lt. hn+2) Then
            ky = j-1
         Else
            ky = j-1-n
         Endif
         Do i = 1,n+2,2
            kx = (i-1)/2
            kr = Sqrt(kx**2+ky**2+kz**2)
            if(kx.eq.0 .and. ky <=0 .and. kz <=0) cycle;!write(*,*)'NO!!!';continue
            if(kx.eq.0 .and. ky >0 .and. kz <0) cycle;!write(*,*)'NO!!!';continue
 
            If (kr .Ne. 0.) Then
               k1  = Ceiling(kr)
               k2  = k1+1
               w1  = 1. ! k1-kr
               w2  = 0. ! 1-w1
               pow = map(i,j,k)
               pst(1,k1)=pst(1,k1)+w1*pow
               pst(2,k1)=pst(2,k1)+w1*pow**2
               pst(3,k1)=pst(3,k1)+w1 ! Count the number of elements
               pst(1,k2)=pst(1,k2)+w2*pow
               pst(2,k2)=pst(2,k2)+w2*pow**2
               pst(3,k2)=pst(3,k2)+w2 ! Count the number of elements
               if (kr<2) write(*,*) 'kr =',kr, 'k1 =',k1,'k2=',k2
               if (kr<2) write(*,*) 'w1 =',w1,'kw=',w2, 'pst(k1) =',pst(1,k1),'pst(k2) =',pst(1,k2) 
            Else
               !write(*,*) 'kr =',kr, 'pst(kr)=', pst(1,1) 
            Endif
         Enddo
      Enddo
   End Do
   !$omp end parallel do

   Ave(1) = map(1,1,1)
   Do k = 1,n
      If (pst(3,k) .Eq. 0) Then
         Ave(k) = 0
      Else 
         Ave(k) = pst(1,k)/pst(3,k)
      Endif
   Enddo
   !write(*,*) 'Ave(1) =', Ave(1), 'Ave(2) =', Ave(2)
   !Close(40)
   !Write(*,*) '3D power spectra written in ',Trim(file)
   Write(*,*) '3D Average Done '

   return

End Subroutine Get_KAverage


!*********************************************
!======================================================
 Subroutine kr2BinAverage(map,Ave,n,normp,normk)
!======================================================
! Calculates the Power Spectrum and Dumps it into a file

   Real                              :: normp,normk ! useless for now...

   Integer                           :: n
   Real, Dimension(1:n+2,1:n,1:n)    :: map
   Character(Len=120)                :: file

   Real, Dimension(2,3*(n/2)**2 + 1) :: pst
   Real, Dimension(3*(n/2)**2 + 1)   :: Ave ! kr max = 3*(256**2), and must count kr = 0 
   Real                              :: pow

   Integer                           :: i,j,k,hn,kz,kx,ky,kr2

   If (mod(n,2) .NE. 0) Then
      Write(*,*) ' n must be even'
      Stop
   End If


   !! == Dump power spectra
   !Open(unit=40,file=file,status='replace')
   hn  = n/2
   pst = 0.0

   !$omp  parallel do default(shared) private(i,j,k,kx,ky,kz,kr2,pow)
   Do k = 1,n
      If (k .Lt. hn+2) Then
         kz = k-1
      Else
         kz = k-1-n
      Endif
      Do j = 1,n
         If (j .Lt. hn+2) Then
            ky = j-1
         Else
            ky = j-1-n
         Endif
         Do i = 1,n+2,2
            kx = (i-1)/2
            kr2 = kx**2+ky**2+kz**2 ! = [0,1,2,3, ... 3*(n/2)**2+1]                      

            if(kx.eq.0 .and. ky <=0 .and. kz <=0) cycle;
            if(kx.eq.0 .and. ky >0 .and. kz <0) cycle;
            If (kr2 .Ne. 0.) Then             

               pow = map(i,j,k)
               pst(1,kr2+1) = pst(1,kr2+1) + pow
               pst(2,kr2+1) = pst(2,kr2+1) + 1

#ifdef kr2debug
               if (kr2.eq.2) then 
                  write(*,*) 'kr2 =',kr2, 'kx =',kx,'ky =',ky,'kz =',kz
                  write(*,*) 'rho =',map(i,j,k)
                  write(*,*) 'Cumul_rho =',pst(1,kr2+1)
                  write(*,*) 'Cumul_weight =',pst(2,kr2+1)
               endif

               if (kr2.eq.3) then 
                  write(*,*) 'kr2 =',kr2, 'kx =',kx,'ky =',ky,'kz =',kz
                  write(*,*) 'rho =',map(i,j,k)
                  write(*,*) 'Cumul_rho =',pst(1,kr2+1)
                  write(*,*) 'Cumul_weight =',pst(2,kr2+1)
               endif

               if (kr2.eq.4) then
                  write(*,*) 'kr2 =',kr2, 'kx =',kx,'ky =',ky,'kz =',kz
                  write(*,*) 'rho =',map(i,j,k)
                  write(*,*) 'Cumul_rho =',pst(1,kr2+1)
                  write(*,*) 'Cumul_weight =',pst(2,kr2+1)
               endif

               !if (kr2<2) write(*,*) 'w1 =',w1,'kw=',w2, 'pst(k1) =',pst(1,k1),'pst(k2) =',pst(1,k2) 

#endif

            Else
               !write(*,*) 'kr =',kr, 'pst(kr)=', pst(1,1) 
            Endif
         Enddo
      Enddo
   End Do
   !$omp end parallel do
   Ave(1) = map(1,1,1)
   kr2weight(1) = 1
   Do kr2 = 2,3*(n/2)**2+1
      If (pst(2,kr2) .Eq. 0) Then
         Ave(kr2) = 0
         kr2weight(kr2) = 0
      Else 
         !write(*,*) 'kr2 =', kr2, 'pst(1) =', pst(1,kr2), 'weigth = ', pst(2,kr2) 
         !pause
         Ave(kr2) = pst(1,kr2)/pst(2,kr2)
         kr2weight(kr2) = pst(2,kr2)
      Endif
   Enddo


#ifdef kr2debug
   write(*,*) 'rhoAve(0) =', Ave(1)
   write(*,*) 'rhoAve(1) =', Ave(2)
   write(*,*) 'rhoAve(2) =', Ave(3)
   write(*,*) 'rhoAve(3) =', Ave(4)
   write(*,*) 'rhoAve(4) =', Ave(5)
   write(*,*) 'rhoAve(5) =', Ave(6)
   write(*,*) 'rhoAve(6) =', Ave(7)
   write(*,*) 'rhoAve(7) =', Ave(8)
   write(*,*) 'rhoAve(8) =', Ave(9)
   write(*,*) 'rhoAve(9) =', Ave(10)
#endif

   !Close(40)
   !Write(*,*) '3D power spectra written in ',Trim(file)
   Write(*,*) '3D Average Done '

   return

End Subroutine kr2BinAverage


!*********************************************


!*********************************

! Usefull after a Fourier Transform the density, this routine rearranges indices to run from [1,hc+1],[1,nc],[1,nc]
! Corresponding to physical wave numbers of [1,hc+1],[-hc+1,hc] and [-hc+1,hc]

! d must be in 3-D Fourrier space
 
! NOTE : the density will have Real/Imaginary components arranged as follow:
! Re(rho) from kx -> [1,hc+1], then Im(rho) from kx -> [hc+2,nc+2], 

! Ex: After Rearrange(d,nc), a slice through the kx axis (ky = kz = 0) is given by:
! do k = 1:hc
!   d(k,hc,hc)
! enddo

subroutine Rearrange(d,nc)
  implicit none
  
  integer nc,hc,FB
  integer, parameter :: nt=4
  
  real, dimension(nc+2,nc,nc) :: d
  real, dimension(nc+2,nc,nc, nt) :: rhot
  
  integer it,i,j,k,kpt
  real kx,ky,kz,pi
  


  pi=acos(-1.)
  hc=nc/2
  kpt=nc/nt

  !$omp parallel do default(shared) &
  !$omp& private(it,i,j,k,kr,kx,ky,kz,k1,k2,w1,w2,pow)
  do it=1,nt
     rhot(:,:,:,it)=0
     !weightt(:,it)=0
     do k=1+(it-1)*kpt,min(nc,it*kpt)
        if (k .lt. hc+2) then
            kz=k-1 ! kz -> [0,hc]
        else
            kz=k-1-nc ! kz -> [-hc+1, -1] 
        endif

        ! Now, kz runs in the range [-hc+1,hc]
        ! I will shift values of kz so that they run from [1,nc]
        kz = kz + hc
        
        do j=1,nc
           if (j .lt. hc+2) then
              ky=j-1 ! ky -> [0,hc]
           else
              ky=j-1-nc ! ky -> [-hc+1, -1] 
           endif

           ! Now, ky runs in the range [-hc+1,hc]
           ! I will shift values of ky so that they run from [1,nc]
           ky = ky + hc


          ! kx is particular: (i odd = real, i even = imaginary)
          ! so I will bring Re(rho) from kx -> [1,hc+1], then Im(rho) from kx -> [hc+2,nc+2]


           !WEIRD!!! Kx should run from 0 to hc

           do i=1,nc+1,2 
              kx=(i+1)/2 

              ! kx runs from [1,hc+1]. 

              ! Filling Re(rho) from [1,hc+1]
              rhot(kx,ky,kz,it)=d(i,j,k)

              ! Filling Im(rho) from [hc+2,nc+2]
              rhot(hc+kx+1,ky,kz,it)=d(i+1,j,k)
                     
              !write(*,*) 'Test : kx = ', kx, 'ky = ', ky, 'kz = ',kz
              !pause


           enddo
        enddo
     enddo
  enddo
  !$omp end parallel do
  

  !! Merge density from threads

  d=0
  do it=1,nt  
     d=d+rhot(:,:,:,it)
  enddo

  !write(*,*) 'Density : ', maxval(d), ' > d > ',minval(d)
  
 
  return
end subroutine Rearrange

!*****************************************



!! 1- FFT-Inverse a Power Spectrum to.get the Correlation Function
!! 2- AbsoluteValue-Square the Correlation Function
!! 3- FFT to get a Covariance(delta_k)

subroutine GetCov(PS,Cov,nc)
  implicit none

  integer nc
  real, dimension(nc+2,nc,nc) :: PS, Cov  
  integer i,j,k

  write(*,*) 'PS before backward FFT : ', maxval(PS),' > PS > ' ,minval(PS)

  !do i = 1,nc
     !write(*,*) 'PPShell(', i, ', hc, hc) = ', d(i,hc,hc)  
  !enddo

  call sfft3_r2c(PS,nc,-1) 


  !! Now PS is a "mostly-real" function, going from 1 to nc 
  !! in all 3 dimensions.Let's absolute-Square the power spectrum

  Cov = abs(PS)**2

#ifdef debug_AngCov
  write(*,*) 'PS after backward FFT : ', maxval(PS),' > PS > ' ,minval(PS)

  write(*,*) 'Xi^2(0,0,0) = ', Cov(1,1,1)
  write(*,*) 'Xi^2(1,0,0) = ', Cov(2,1,1)
  write(*,*) 'Xi^2(0,1,0) = ', Cov(1,2,1)
  write(*,*) 'Xi^2(0,0,1) = ', Cov(1,1,2)

  write(*,*) 'Cov before FFT : ', maxval(Cov),' > Cov > ' ,minval(Cov)
#endif

  call sfft3_r2c(Cov,nc,1) 

#ifdef debug
  write(*,*) 'Cov after FFT : ', maxval(Cov),' > Cov > ' ,minval(Cov)
#endif

  !Now Cov is in k-space, orderer as the fftw wishes it.

  return
end subroutine GetCov

!*******************************
subroutine PID_analysis
  implicit none

  integer i,j,k
  integer(kind=8) PID1,PID2

#ifdef GET_PID
  do i=1,100
     write(*,*) PID(i)
  enddo
#endif

#ifdef HALOPID
  do i = 1,np
    !write(*,*) halo_pid(:,i) 
    do j = 1,N_p
       do k = j+1,N_p
          PID1 = halo_pid(j,i)
          PID2 = halo_pid(k,i)
          if((PID1 .eq. PID2) .and. PID1 .ne. 0 ) then
             write(*,*)'Halo PIDs :',halo_pid(:,i)
             write(*,*)'Found identical PIDs:' , PID1,PID2
          endif
       enddo
    enddo
  enddo

#endif

  write(*,*) '*******************'
  write(*,*) '*done PID_analysis*'
  write(*,*) '*******************'

  return
end subroutine PID_analysis
!*******************************
!*******************************
subroutine ps3_r2c_err_NGP(d,ps,err,nc)
  implicit none

  integer nc,hc
  integer, parameter :: nt=4

  real, dimension(nc) :: ps,weight,err
  real, dimension(nc+2,nc,nc) :: d
  real, dimension(nc,nt) :: pst,weightt,errt

  integer it,i,j,k,kpt
  real kr,kx,ky,kz,k1,k2,w1,w2,pow,pi

  character(len=MSL) :: ofile

  pi=acos(-1.)
  hc=nc/2
  kpt=nc/nt

  call sfft3_r2c(d,nc,1)

  !$omp parallel do default(shared) &
  !$omp& private(it,i,j,k,kr,kx,ky,kz,k1,k2,w1,w2,pow)
  do it=1,nt
     errt(:,it)=0
     pst(:,it)=0
      weightt(:,it)=0
     do k=1+(it-1)*kpt,min(nc,it*kpt)
        if (k .lt. hc+2) then
            kz=k-1
        else
            kz=k-1-nc
        endif
        do j=1,nc
           if (j .lt. hc+2) then
               ky=j-1
           else
               ky=j-1-nc
           endif
           do i=1,nc+2,2
              kx=(i-1)/2
              kr=sqrt(kx**2+ky**2+kz**2)
              if(kx.eq.0 .and. ky <=0 .and. kz <=0) cycle;!write(*,*)'NO!!!';continue
              if(kx.eq.0 .and. ky >0 .and. kz <0) cycle;!write(*,*)'NO!!!';continue
                if (kr .ne. 0) then
                  k1=ceiling(kr)
                  k2=k1+1
                  w1=1!k1-kr
                  w2=0!1-w1
                  pow=(d(i,j,k)/nc**3)**2+(d(i+1,j,k)/nc**3)**2
                  weightt(k1,it)=weightt(k1,it)+w1
                  pst(k1,it)=pst(k1,it)+w1*pow
                  errt(k1,it)=errt(k1,it)+w1*pow**2
                  weightt(k2,it)=weightt(k2,it)+w2
                  pst(k2,it)=pst(k2,it)+w2*pow
                  errt(k2,it)=errt(k2,it)+w2*pow**2

                  !if(kr <= 1)write(*,*)'k=',kx,ky,kz,pow*lbox**3

              endif
           enddo
        enddo
     enddo
  enddo
  !$omp end parallel do

  !! Merge power spectrum from threads
  ps=0
  weight=0
  do it=1,nt
     ps=ps+pst(:,it)
     err=err+errt(:,it)
     weight=weight+weightt(:,it)
  enddo

  !! Divide by weights
  !! Convert P(k) to \Delta^2(k)
  !! Store k in ps(1,i)

  !ps(1) = d(1,1,1)
  !err(1) = 0
  !write(*,*)'k',1, 'ps',ps(1),'weight',1

  !ofile = 'source/Weight128-NGP.dat'
  !open(unit=13,file=ofile)

  do k=1,nc
     if (weight(k) .ne. 0) then        
        ps(k) = ps(k)/weight(k)
        err(k)= 4*pi*(k**3)*(err(k)/weight(k) - ps(k)**2)
        ps(k) = 4*pi*(k**3)*ps(k)
       ! write(*,*)'k',k, 'ps',ps(k),'weight',weight(k)
     endif
     !write(13,*) weight(k)
  enddo
  
  !close(13)

  return
end subroutine ps3_r2c_err_NGP
!*******************************
subroutine AngleAverageKSpace(d,ps,err,nc)
  implicit none

  integer nc,hc
  integer, parameter :: nt=4

  real, dimension(nc) :: ps,weight,err
  real, dimension(nc+2,nc,nc) :: d
  real, dimension(nc,nt) :: pst,weightt,errt

  integer it,i,j,k,kpt
  real kr,kx,ky,kz,k1,k2,w1,w2,pow,pi

  character(len=MSL) :: ofile

  pi=acos(-1.)
  hc=nc/2
  kpt=nc/nt

  !call sfft3_r2c(d,nc,1)

  !$omp parallel do default(shared) &
  !$omp& private(it,i,j,k,kr,kx,ky,kz,k1,k2,w1,w2,pow)
  do it=1,nt
     errt(:,it)=0
     pst(:,it)=0
      weightt(:,it)=0
     do k=1+(it-1)*kpt,min(nc,it*kpt)
        if (k .lt. hc+2) then
            kz=k-1
        else
            kz=k-1-nc
        endif
        do j=1,nc
           if (j .lt. hc+2) then
               ky=j-1
           else
               ky=j-1-nc
           endif
           do i=1,nc+2,2
              kx=(i-1)/2
              kr=sqrt(kx**2+ky**2+kz**2)
              if(kx.eq.0 .and. ky <=0 .and. kz <=0) cycle;!write(*,*)'NO!!!';continue
              if(kx.eq.0 .and. ky >0 .and. kz <0) cycle;!write(*,*)'NO!!!';continue
                if (kr .ne. 0) then
                  k1=ceiling(kr)
                  k2=k1+1
                  w1=1!k1-kr
                  w2=0!1-w1
                  pow=(d(i,j,k)/nc**3)**2+(d(i+1,j,k)/nc**3)**2
                  weightt(k1,it)=weightt(k1,it)+w1
                  pst(k1,it)=pst(k1,it)+w1*pow
                  errt(k1,it)=errt(k1,it)+w1*pow**2
                  weightt(k2,it)=weightt(k2,it)+w2
                  pst(k2,it)=pst(k2,it)+w2*pow
                  errt(k2,it)=errt(k2,it)+w2*pow**2
              endif
           enddo
        enddo
     enddo
  enddo
  !$omp end parallel do

  !! Merge power spectrum from threads
  ps=0
  weight=0
  do it=1,nt
     ps=ps+pst(:,it)
     err=err+errt(:,it)
     weight=weight+weightt(:,it)
  enddo

  !! Divide by weights
  !! Convert P(k) to \Delta^2(k)
  !! Store k in ps(1,i)

  !ps(1) = d(1,1,1)
  !err(1) = 0
  !write(*,*)'k',1, 'ps',ps(1),'weight',1

  !ofile = 'source/Weight128-NGP.dat'
  !open(unit=13,file=ofile)

  do k=1,nc
     if (weight(k) .ne. 0) then        
        ps(k) = ps(k)/weight(k)
        err(k)= 4*pi*(k**3)*(err(k)/weight(k) - ps(k)**2)
        ps(k) = 4*pi*(k**3)*ps(k)
       ! write(*,*)'k',k, 'ps',ps(k),'weight',weight(k)
     endif
     !write(13,*) weight(k)
  enddo
  
  !close(13)

  return
end subroutine AngleAverageKSpace
!*******************************



end program cicpow
