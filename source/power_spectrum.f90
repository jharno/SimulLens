!! Modified by Joachim Harnois-Deraps
!! Power-spectrum modified on Hy Trac's + Ting Ting Lu's code.
!! Size: xv = 6*(nc/2)**3*(4 bytes )
!!       ll = (nc/2)**3*(4 bytes )
!!        d = nc**3*(4 bytes )
!! Total    = (7.5)*(nc)**3 bytes 
!! Ex.: nc = 2048, Total = 65 Gb      
program cicpow
  implicit none
  include 'source/BAO.fh'

  !! np should be set to nc (1:1) or hc (1:2)
  integer, parameter :: nt=8

#ifndef HALO
  integer(kind=8), parameter :: np=hc**3
  real, parameter :: mp=int(nc,kind=8)**3/np !nc**3/np
  !real, dimension(6,np) :: xv
  !integer, dimension(np) :: ll

  real, allocatable, dimension(:,:) :: xv
  real, allocatable, dimension(:)   :: ll

#ifdef GETPID
  integer(kind=8), dimension(np) :: PID
#endif

#else
  integer(kind=8), parameter :: np_max = hc**3 
  integer(kind=8) :: np  
  real :: mp
  !real, dimension(6,np_max) :: xv
  !integer, dimension(np_max) :: ll

  real, allocatable, dimension(:,:) :: xv
  real, allocatable, dimension(:)   :: ll

#ifdef HALOPID  
  integer(kind=8), dimension(10,10000) :: halo_pid
#endif
#endif



  !! cubepm
!  integer, parameter :: nodes_dim=2
  real, parameter  :: ncc=nc/nodes_dim !! number of cells / cubic 
  real, parameter  :: rnc = nc

  !! Dark matter arrays
  integer, dimension(nc,nc) :: hoc
  !integer, dimension(2,nc,nc,nt) :: htoc
  integer, allocatable, dimension(:,:,:,:) :: htoc



  !! Power spectrum arrays
  real, dimension(nc) :: ps,err
  real, dimension(nc) :: PS_Ave!,  CrossCorrAve, N_CrossCorrAve
  !real, dimension(3*(nc/2)**2+1) :: CrossCorrAve, N_CrossCorrAve
  !real, dimension(nc+2,nc,nc) :: d 
  real, allocatable, dimension(:,:,:) :: d

 

 !equivalence (d, d_COMPLEX, psV_COMPLEX, N_CrossCorr,N_psShell_Outer,N_psShell_Outer_COMPLEX, N_CrossCorr_COMPLEX) 
  
 
#ifdef FFT_TEST
 real, dimension(nc+2,nc,nc) :: Sync, Gauss, DeltaFunction, N_Hat, Hat
 real, dimension(nc) :: Xsi, HatAve, GaussAve, DeltaAve, FFT_Hat
#endif

  !real, dimension(0:hc,-hc+1:hc,-hc+1:hc) :: Hat
  integer, parameter :: MSL=100
  character (len=MSL) :: ofile,LOS,zstring

  integer nploc(nn),fstat
  real pi,omegak

  !! variables in CubePM
  integer cubepm_nts,cubepm_cur_checkpoint,cubepm_cur_projection,cubepm_cur_halofind,i1,j1,k1, node_coords(3)
  real  cubepm_a,cubepm_t,cubepm_tau,cubepm_dt_f_acc,cubepm_dt_pp_acc,cubepm_dt_c_acc, cubepm_mass_p

#ifdef HALO
  !real, dimension(np_max) :: halo_mass!, cubepm_radius_scale, halo_mass_pp, halo_mass1!,cubepm_v_disp
  real, allocatable, dimension(:) :: halo_mass
  !real, dimension(3,np_max) :: cubepm_halo_pos, peak, cubepm_l, var,cubepm_v_disp
#endif

  character(*), parameter :: output_form=    'unformatted' !'binary'
  character(*), parameter :: output_access=  'stream' !'sequential'
  character(*), parameter :: input_form=     'binary'

  !character(*), parameter :: output_form=    'unformatted'
  !character(*), parameter :: output_access=  'stream' 
  !character(*), parameter :: input_form=     'unformatted'

#ifdef HALO
  allocate(xv(6,np_max))
  allocate(ll(np_max))
  allocate(halo_mass(np_max))
#else
  allocate(xv(6,np))
  allocate(ll(np))
#endif

  allocate(htoc(2,nc,nc,nt))
  allocate(d(nc+2,nc,nc))


  !common /rarr/ d,xv,ps 
  !common /iarr/ ll,htoc,hoc

  pi=acos(-1.)
  omegak = 1-omegam-omegav

  !$ call omp_set_num_threads(nt)
    

#ifdef shotnoise

#ifdef HALO
  np = 20000
  write(*,*) 'Using', np , 'halos for shotnoise'
  mp = real(nc**3)/real(np)
  
#endif

  call random_number(xv(1:3,1:np))

  xv(1:3,1:np)=xv(1:3,1:np)*nc

  call cic
  call powerspectrum
  call writeps_rdm
  stop
#endif


#ifndef RANDOM
#ifndef LegendreTest

  !call getarg(1,LOS)
  call getarg(1,zstring)
  call readdm 
  call cic
#ifdef debug
       do i1 =1,10
          write(*,*) 'd(',i1,',1,1)=',d(i1,1,1)
       enddo
#endif


#ifdef dump_density
  !write(*,*)'TMP::WRINTING DENSITY TO FILE'
  !open(unit = 33,file = '/mnt/scratch-3week/jharno/Densities/density.dat')
  ofile = '/mnt/scratch-3week/jharno/Densities/density-'//Run//'-'//Version//'.dat'
  !open (unit=33,file=ofile,status='replace',iostat=fstat,form='binary')
  open (unit=33,file=ofile,status='replace',iostat=fstat)

  if(fstat /= 0) then 
        write(*,*) 'could not open the file properly' 
        write(*,*) 'ofile =' , ofile
        write(*,*) 'max(d) =' , maxval(d) 
  else 
        write(*,*) 'WRITING DENSITY TO FILE AND STOP'
        !write(33) d
        do k1 = 1,nc
           do j1 = 1,nc
              do i1 = 1,nc
                 !write(*,*)i1,j1,k1,d(i1,j1,k1)
                 write(33, '(1f18.8)') d(i1,j1,k1)
              enddo
           enddo
        enddo
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

#ifndef LegendreTest
  call powerspectrum
  call writeps
#endif

deallocate(xv)
deallocate(ll)
deallocate(d)
deallocate(htoc)

#ifdef HALO
   deallocate(halo_mass)
#endif

contains

  subroutine readdm
    implicit none
    character*100 fn, fn2, fn3

    character*8 t1,t2
    integer i,j,ii, Nmax, nh_local, nh_total, high_mass_counter
    integer low_mass_counter, blocksize,num_writes,nplow,nphigh
    integer(kind=8) ip
    real HubbleScale
    real Conversion
    character (len=4) :: rank_s
    character (len=MSL) :: ofile,PIDofile!, zstring
    character (len=MSL) :: ifile,PIDifile
#ifdef HALO
    real(4), dimension(28) :: halo_input_buffer
#endif

#ifdef MERGED_LIST
    real, dimension(28):: reading_buffer
#endif

#ifndef MERGED_LIST

    !! Read particle data file from decomposed volume
    ip=0
    nh_total=0     
    high_mass_counter=0
    low_mass_counter=0
    do i=1,nn

       nh_local = 0      

       write(*,*) 'Reading Node ',i
       write(rank_s,'(i4)') i-1
       rank_s=adjustl(rank_s)
       !write(zstring,'(f5.3)') z3dps       ! Need (f6.3) for Z > 10
       !write(zstring,'(f5.3)') z_str
        

#ifdef HALO
       ifile=trim(zstring)//'halo'//rank_s(1:len_trim(rank_s))//".dat"
#else 
       ifile=trim(zstring)//'xv'//rank_s(1:len_trim(rank_s))//".dat"
#endif

#ifdef IC
       ifile='xv'//rank_s(1:len_trim(rank_s))//".ic"
#endif

       ofile=proj_path//'/'//trim(ifile)
       !ofile=proj_path//trim(LOS)//'/'//trim(ifile)
       !ofile='/scratch/jharno/NewLensing_2048/LOS'//trim(LOS)//'/'//trim(ifile)
       !ofile='/mnt/scratch-2week/jharno/power_law/z_i=20/'//trim(ifile)
       !ofile='/mnt/scratch-2week/jharno/MoveBack/without/'//trim(ifile)
       !ofile='/mnt/scratch-3week/akhazr/output/box200/'//trim(LOS)//'/'//trim(ifile)
       !ofile='/mnt/scratch-3week/akhazr/scinet1024/Run-'//trim(LOS)//'/'//trim(ifile)
       !ofile=proj_path//Version//'/out/RUN-'//Run//'/'//trim(ifile)


!**********
       write(*,*) 'opening ',ofile
       open (unit=12,file=ofile,status='old',iostat=fstat,form=input_form)
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
          !peak(:,nh_total)=halo_input_buffer(1:3)
          xv(:,nh_total)=halo_input_buffer(4:9)
          !cubepm_l(:,nh_total)=halo_input_buffer(10:12)
          !!cubepm_v_disp(nh_total)=halo_input_buffer(13)
          !cubepm_v_disp(:,nh_total)=halo_input_buffer(13:15)
          !cubepm_radius_scale(nh_total)=halo_input_buffer(16)
          halo_mass(nh_total)=halo_input_buffer(17)
          !halo_mass_pp(nh_total)=halo_input_buffer(18)
          !halo_mass1(nh_total)=halo_input_buffer(19)
          !var(:,nh_total)=halo_input_buffer(20:22)
          !write(*,*) halo_mass(nh_total)

          !if(halo_mass(nh_total)>2000) high_mass_counter=high_mass_counter+1
          !if(halo_mass(nh_total)<550) low_mass_counter=low_mass_counter+1

          !!!!!!!!!!!!!!!!!!!!
          ! High mass selector
          !!!!!!!!!!!!!!!!!!!!

          !if(halo_mass(nh_total)>2000) then 
          !   high_mass_counter=high_mass_counter+1
          !else
          !   nh_total=nh_total-1
          !   nh_local=nh_local-1
          !   continue
          !endif

          !!!!!!!!!!!!!!!!!!!
          ! Low mass selector
          !!!!!!!!!!!!!!!!!!!

          if(halo_mass(nh_total)<550) then 
             low_mass_counter=low_mass_counter+1
          !   !write(*,*) 'halo number      = ', nh_total
          !   !write(*,*) 'low_mass_counter = ', low_mass_counter
          !   !write(*,*) 'halo mass        = ', halo_mass(nh_total)
          !   !write(*,*) 'halo xv          = ', xv(:,nh_total)
          else
             nh_total=nh_total-1
             nh_local=nh_local-1
             continue
          endif


#ifdef debug_halo
          !write(*,*)'peak=',peak(:,nh_total)
          write(*,*)'xv=',xv(:,nh_total)
          !write(*,*)'cubepm_l=',cubepm_l(:,nh_total)
          !write(*,*)'cubepm_v_disp=',cubepm_v_disp(nh_total)
          !write(*,*)'cubepm_radius_scale=',cubepm_radius_scale(nh_total)
          !write(*,*)'halo_mass=',halo_mass(nh_total)
          !write(*,*)'halo_mass_pp=',halo_mass_pp(nh_total)
          !write(*,*)'halo_mass1=',halo_mass1(nh_total)
          !write(*,*)'var=',var(:,nh_total)
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


       blocksize=(32*1024*1024)/24
       num_writes=nploc(i)/blocksize+1

       do ii=1,num_writes
          nplow=(ii-1)*blocksize+1
          nphigh=min(ii*blocksize,nploc(i))
          !!      print *,rank,nplow,nphigh,np_local
          do j=nplow,nphigh
             read(12) xv(:,ip+j-1)
          enddo
       enddo
       !read(12) xv(1:6,ip+1:ip+nploc(i))
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

!#ifdef debug
       do i1=1,30
          write(*,*) 'x=',xv(:3,i1)
       enddo
!#endif

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
       !PIDofile=proj_path//Version//'/out/RUN-'//Run//'/'//trim(ifile)
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

    write(*,*) 'n high_mass halos = ', high_mass_counter
    write(*,*) 'n low_mass halos = ', low_mass_counter
 

#else

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! If working with a MERGED LIST !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef HALO
    ip = 0

    !write(zstring,'(f5.3)') z3dps
    ifile=trim(zstring)//'halo.dat' ! xv.dat
    !ofile='/cita/d/scratch-2week/jharno/PID/'//trim(ifile)
    ofile=proj_path//Version//'/out/RUN-'//Run//'/'//trim(ifile)

    write(*,*) 'opening ',ofile
    open (unit=12,file=ofile)

    do
       read(12,'(20f20.10)',end=112,err=113) reading_buffer
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
   
    !write(*,*) 'xv halos=' , xv(1:3,1:np+1)

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

    write(*,*) 'testing:', real(nc,kind=8)**3, np,  real(nc,kind=8)**3/real(np)

    !**********************************
    ! Reassigning the mass of the halos
    mp = real(nc,kind=8)**3/real(np)
    !**********************************

    write(*,*) 'halo mass = ', mp
    write(*,*) 'total mass = ', mp*np
    write(*,*) 'total mass pp= ', 8*hc**3
    write(*,*) 'average mass per grid cell= ', mp*np/real(nc,kind=8)**3

#endif

    write(*,*) '*************'
    write(*,*) '*Done readdm*'
    write(*,*) '*************'
    return
  end subroutine readdm

  !************************
  subroutine writeps
    implicit none
    integer k
    character*150 fn1!, zstring

    !write(zstring,'(f5.3)') z       ! Need (f6.3) for Z > 10

#ifdef NGP_PS

    !fn1=PS_dir//'z_'//Z//'/ngppow_'//RedShift//'_Run'//Run//'-'//Version//'.dat'
    fn1=proj_path//trim(LOS)//'/'//trim(zstring)//'ngpps.dat'
    !fn1='/cita/d/scratch-2week/jharno/MoveBack/without/'//trim(zstring)//'ngpps.dat'
    !fn1 = Test_dir//'z_0.5/'//trim(zstring)//'ngpps_Run'//trim(LOS)//'.dat'

#ifdef TEST
    !fn1=Test_dir//'ngppow_'//RedShift//'_Run'//Run//'-'//Version//'-with.dat'
    !fn1 = proj_path//Version//'/out/RUN-'//Run//'/'//Z//'ngpps.dat'
    fn1=Test_dir//'ngppow_'//RedShift//'_Run'//Run//'-'//Version//'-with.dat'
#endif

#ifdef HALO
    !fn1=PS_dir//'z_'//Z//'/ngppow_'//RedShift//'_Run'//Run//'-'//Version//'-halo.dat'
    !fn1=Test_dir//'z_0.5/'//trim(zstring)//'ngpps_Run'//trim(LOS)//'-light-halo.dat'
    fn1 = PS_dir//trim(zstring)//'halops.dat'
#endif


#ifdef Kaiser
    fn1=PS_dir//'z_'//Z//'/npgps_'//RedShift//'_Run'//Run//'-'//Version//'-RSD.dat'
#endif

#ifdef RANDOM
#ifdef poisson
    fn1=PS_dir//'z_'//Z//'/npgps_'//RedShift//'_Run'//Run//'-'//Version//'-poisson-80000.dat'
#else
    fn1=PS_dir//'z_'//Z//'/npgps_'//RedShift//'_Run'//Run//'-'//Version//'-RDM.dat'
#endif
#endif

#ifdef GAUSS
    fn1=PS_dir//'z_'//Z//'/ngppow_'//RedShift//'_Run'//Run//'-'//Version//'-GAUSS.dat'
#ifdef TEST
    fn1=Test_dir//'ngppow_'//RedShift//'_Run'//Run//'-'//Version//'-GAUSS.dat'
#endif
#endif

#ifdef power_law
    fn1=PS_dir//'z_'//Z//'/ngppow_'//RedShift//'_Run'//Run//'-'//Version//'-power_law_z_i=20.dat'
#endif

#else

    fn1=PS_dir//'z_'//Z//'/cicpow_'//RedShift//'_Run'//Run//'-'//Version//'.dat'
#ifdef TEST
    fn1=Test_dir//'cicpow_'//RedShift//'_Run'//Run//'-'//Version//'.dat'
    !fn1 = proj_path//Version//'/out/RUN-'//Run//'/'//Z//'00cicps.dat'
#endif

#ifdef HALO
    fn1=PS_dir//'cicpow_'//RedShift//'_Run'//Run//'-'//Version//'-halo.dat'
#endif


#ifdef Kaiser
    fn1=PS_dir//'z_'//Z//'/cicpow_'//RedShift//'_Run'//Run//'-'//Version//'-RSD.dat'
#endif

#ifdef RANDOM
    fn1=PS_dir//'z_'//Z//'/cicpow_'//RedShift//'_Run'//Run//'-'//Version//'-RDM.dat'
#endif

#ifdef GAUSS
    fn1=PS_dir//'z_'//Z//'/cicpow_'//RedShift//'_Run'//Run//'-'//Version//'-GAUSS.dat'
#ifdef TEST
    fn1=Test_dir//'cicpow_'//RedShift//'_Run'//Run//'-'//Version//'-GAUSS.dat'
#endif
#endif

#ifdef power_law
    fn1=PS_dir//'z_'//Z//'/cicpow_'//RedShift//'_Run'//Run//'-'//Version//'-power_law_z_i=20.dat'
#endif

#endif

    !! Output power spectrum
    !! First column is physical k
    !! Second column is \Delta^2

    open(11,file=fn1)
    do k=1,hc !hc+1
#ifdef NGP_PS
       write(11,*) 2*pi/lbox*(k), ps(k) ,err(k)
#else
       write(11,*) 2*pi/lbox*(k-1), ps(k) ,err(k)
#endif
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

#ifdef NGP_PS

#ifdef HALO
    fn1=PS_dir//'ngppow-20000halos-200Mpc-1024-shotnoise.dat'
#else    
    fn1=PS_dir//'ngppow-shotnoise-600Mpc-512.dat'
#endif

#else

#ifdef HALO
    fn1=PS_dir//'cicpow-halo-shotnoise.dat'
#else    
    fn1=PS_dir//'cicpow-shotnoise.dat'
#endif

#endif

    !! Output power spectrum
    !! First column is physical k
    !! Second column is \Delta^2

    open(11,file=fn1)
    do k=1,hc !hc+1
#ifdef NGP_PS
       write(11,*) 2*pi/lbox*(k), ps(k) ,err(k)
#else
       write(11,*) 2*pi/lbox*(k-1), ps(k) ,err(k)
#endif
    enddo
    close(11)

    write(*,*) 'Wrote ', fn1

    write(*,*) '******************'
    write(*,*) '*done writeps_rdm*'
    write(*,*) '******************'

    return
  end subroutine writeps_rdm



  !************************************
  !*** cic computes the density 'd' ***
  !************************************

  subroutine cic
    implicit none
    integer, parameter :: kpt=nc/nt

#ifdef HALO
    integer :: npt
#else
    integer, parameter :: npt=np/nt
#endif
    
    integer it,i,j,k,ip,OutBound
    real toe

#ifdef HALO
    npt=np/nt
    write(*,*) 'Performing cic on', npt, 'haloes and ', nt ,' threads'
#endif

    !! Construct chaining lists in parallel
    !$omp parallel do default(shared) private(it,ip,j,k) int(nc,kind=8)
    OutBound = 0
    do it=1,nt
       htoc(:,:,:,it)=0
       do ip=1+(it-1)*npt,min(np,int(it*npt, kind=8))
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
          if((j <= 0)) then
#ifdef debug_cic
             write (*,*) '#### PROBLEM!!! (j = floor(xv(2,ip))+1) =',j
             write (*,*) 'Enforcing Periodic BC Manually'
#endif
             j = j+nc
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
          if((k <= 0)) then
#ifdef debug_cic
             write (*,*) '#### PROBLEM!!! (k = floor(xv(2,ip))+1) =',k
             write (*,*) 'Enforcing Periodic BC Manually'
#endif
             k = k+nc
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
    integer i,j,k, n_counts

#ifdef GAUSS
    real, dimension(nc) :: kgauss
    real, dimension(hc) :: ps_gauss
    character(len=MSL) :: psGaussFile
#endif 
#ifdef RANDOM
    integer k_seed, clock
    integer, dimension(2) :: old, seed
#endif

#ifdef debug
    write(*,*) 'd before fft', maxval(d),minval(d) 
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
             if(d(i,j,k) > 1.0 -(80000.0/nc/nc/nc)) then
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

!*************

    !call ps3_r2c(d,ps,nc)
#ifdef NGP_PS
    call ps3_r2c_err_NGP(d,ps,err,nc) 
#else    
    call ps3_r2c_err(d,ps,err,nc)
#endif

#ifdef debug
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

#ifdef debug
    write(*,*)'Filled kgauss'
#endif

    call GaussRandomField_3d_r2c(d, lbox, nc, kgauss, ps, nc)

    call ps3_r2c(d,ps_gauss,nc)

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



!======================================================
 Subroutine Get_KAverage(map,Ave,n,normp,normk)
!======================================================
! Calculates the Power Spectrum and Dumps it into a file


   Integer                            :: n
   Real                               :: normp,normk
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

            If (kr .Ne. 0.) Then
               k1  = Ceiling(kr)
               k2  = k1+1
               w1  = k1-kr
               w2  = 1-w1
               pow = map(i,j,k)
               pst(1,k1)=pst(1,k1)+w1*pow
               pst(2,k1)=pst(2,k1)+w1*pow**2
               pst(3,k1)=pst(3,k1)+w1 ! Count the number of elements
               pst(1,k2)=pst(1,k2)+w2*pow
               pst(2,k2)=pst(2,k2)+w2*pow**2
               pst(3,k2)=pst(3,k2)+w2 ! Count the number of elements
               !if (kr<2) write(*,*) 'kr =',kr, 'k1 =',k1,'k2=',k2
               !if (kr<2) write(*,*) 'w1 =',w1,'kw=',w2, 'pst(k1) =',pst(1,k1),'pst(k2) =',pst(1,k2) 
            Else
               !write(*,*) 'kr =',kr, 'pst(kr)=', pst(1,1) 
            Endif
         Enddo
      Enddo
   End Do
   !$omp end parallel do

   Ave(1) = map(1,1,1)
   Do k = 2,n
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

!*************************************

subroutine ps3_r2c_err_NGP(d,ps,err,nc)
  implicit none

  integer nc,hc
  integer, parameter :: nt=8

  real, dimension(nc) :: ps,weight,err, weight_x,k_mean
  real, dimension(nc+2,nc,nc) :: d
  real, dimension(nc,nt) :: pst,weightt,errt,weightt_x,k_meant

  integer it,i,j,k,kpt
  real kr,kx,ky,kz,k1,k2,w1,w2,pow,pi

  character(len=MSL) :: ofile

  pi=acos(-1.)
  hc=nc/2
  kpt=nc/nt
  write(*,*) 'Performing P(k) on ',nt,'threads'
  call sfft3_r2c(d,nc,1)

  !$omp parallel do default(shared) &
  !$omp& private(it,i,j,k,kr,kx,ky,kz,k1,k2,w1,w2,pow)
  do it=1,nt
     errt(:,it)=0
     pst(:,it)=0
      weightt(:,it)=0
      weightt_x(:,it)=0
      k_meant(:,it)=0
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
                  k_meant(k1,it)=k_meant(k1,it)+w1*kr
                  weightt(k2,it)=weightt(k2,it)+w2
                  pst(k2,it)=pst(k2,it)+w2*pow
                  k_meant(k2,it)=k_meant(k2,it)+w2*kr
                  errt(k2,it)=errt(k2,it)+w2*pow**2
                  if(kx .eq. 0)then
                     weightt_x(k1,it)=weightt_x(k1,it)+w1
                     weightt_x(k2,it)=weightt_x(k2,it)+w2
                  endif
              endif
   !           if(kr<2) write(*,*)kx,ky,kz,kr,pow

           enddo
        enddo
     enddo
  enddo
  !$omp end parallel do

  !! Merge power spectrum from threads
  ps=0
  weight=0
  weight_x=0
  k_mean=0
  do it=1,nt
     ps=ps+pst(:,it)
     err=err+errt(:,it)
     weight=weight+weightt(:,it)
     weight_x=weight_x+weightt_x(:,it)
     k_mean=k_mean + k_meant(:,it)
  enddo

  !! Divide by weights
  !! Convert P(k) to \Delta^2(k)
  !! Store k in ps(1,i)

  !ps(1) = d(1,1,1)
  !err(1) = 0
  !write(*,*)'k',1, 'ps',ps(1),'weight',1

  ofile = 'source/Weight512-NGP_full.dat'
  open(unit=13,file=ofile)

  do k=1,nc
     if (weight(k) .ne. 0) then        
        ps(k) = ps(k)/weight(k)
        err(k)= 4*pi*(k**3)*(err(k)/weight(k) - ps(k)**2)
        ps(k) = 4*pi*(k**3)*ps(k)
        k_mean(k) = k_mean(k)/weight(k)
       ! write(*,*)'k',k, 'ps',ps(k),'weight',weight(k)
  !      write(*,*) 'k =',k , 'k_mean = ', k_mean(k) , 'weight =',weight(k)
     endif
     write(13,*) weight(k),weight_x(k),2*weight(k) - weight_x(k) ,k_mean(k)
  enddo
  
  close(13)

  return
end subroutine ps3_r2c_err_NGP
!*******************************



end program cicpow
