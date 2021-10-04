!! Angular Covariance written by Joachim Harnois- Deraps
!! Power-spectrum modified on Hy Trac's + Ting Ting Lu's code.
! Compile on TCS with:
program GetHaloField
  implicit none
  include './par.fh'

  !! np should be set to nc (1:1) or hc (1:2)
  integer, parameter :: nt=1
  integer(kind=8), parameter :: np_max = 1000000!hc**3 
  integer(kind=8) :: np  
  real :: mp
  real, dimension(6,np_max) :: xv
  integer, dimension(np_max) :: ll
#ifdef HALOPID  
  integer(kind=8), dimension(10,10000) :: halo_pid
#endif

  !! cubepm
!  integer, parameter :: nodes_dim=2
  real, parameter  :: ncc=nc/nodes_dim !! number of cells / cubic 
  real, parameter  :: rnc = nc

  !! Dark matter arrays
  integer, dimension(nc,nc) :: hoc
  integer, dimension(2,nc,nc,nt) :: htoc


  real, dimension(nc+2,nc,nc) :: d
    
  !real, dimension(0:hc,-hc+1:hc,-hc+1:hc) :: Hat
  integer, parameter :: MSL=100

  integer nploc(nn),fstat
  real pi

  !! variables in CubePM
  integer cubepm_nts,cubepm_cur_checkpoint,cubepm_cur_projection,cubepm_cur_halofind,i1,j1,k1, node_coords(3)
  real  cubepm_a,cubepm_t,cubepm_tau,cubepm_dt_f_acc,cubepm_dt_pp_acc,cubepm_dt_c_acc, cubepm_mass_p

  real, dimension(np_max) :: cubepm_radius_scale, halo_mass, halo_mass_pp, halo_mass1
  real, dimension(3,np_max) :: cubepm_halo_pos, peak, cubepm_l, var,  cubepm_v_disp2
  real, dimension(6,np_max) :: I_xy

  character(len=MSL) LOS, zstring

  common /rarr/ d,xv 
  common /iarr/ ll,htoc,hoc

  pi=acos(-1.)

  write(*,*) '****************'
  write(*,*) '* Halo Program *'
  write(*,*) '****************'
  
  call getarg(1,LOS)
  call getarg(2,zstring)

  call readdm 
!  call cic

  call write_info

contains

  subroutine write_info
    implicit none
    character*100 root, fn, fn2, fn3
    !character (len=MSL) :: zstring
    real sbox, lbox
    integer i

    root = '/scratch/jharno/SimulLens_Random_mix/LOS'//trim(LOS)//'/zs_3_147_mix_231Mpc_1024/'
    !write(zstring,'(f5.3)') z3dps       ! Need (f6.3) for Z > 10

! ---- halo mass

    fn = trim(root)//trim(zstring)//'halomass.dat'
    write(*,*) 'Writing ',fn
    open(11,file=fn,recl=500)
    do i = 1, np
       write(11,*) halo_mass_pp(i)
    enddo
    close(11)


! ---- halo 

    fn = trim(root)//trim(zstring)//'halo_full_info.dat'
    write(*,*) 'Writing ',fn
    open(11,file=fn,recl=500)
    do i = 1, np
       write(11,*) xv(:,i), cubepm_l(:,i), cubepm_v_disp2(:,i), cubepm_radius_scale(i), halo_mass_pp(i), var(:,i), I_xy(:,i)
    enddo
    close(11)


! ---- power spectrum
!    write(*,*) 'Writing ',fn
!    open(11,file=fn,recl=500)
!    do k=2,hc+1
!       kr=2*pi*(k-1)/box
!       write(11,*) kr,pkdm(:,k-1)
!    enddo
!    close(11)




    return

end subroutine write_info

  subroutine readdm
    implicit none
    character*100 fn, fn2, fn3

    character*8 t1,t2
    integer i,j, Nmax, nh_total
    integer(kind=8) ip
    real HubbleScale
    real Conversion
    character (len=4) :: rank_s
    character (len=MSL) :: ofile,PIDofile!,zstring
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
    do i=1,nn

       write(*,*) 'Reading Node ',i
       write(rank_s,'(i4)') i-1
       rank_s=adjustl(rank_s)
       !write(zstring,'(f5.3)') z3dps       ! Need (f6.3) for Z > 10

#ifdef HALO
       ifile=trim(zstring)//'zoomed_halo'//rank_s(1:len_trim(rank_s))//".dat"
#else
       ifile=trim(zstring)//'xv'//rank_s(1:len_trim(rank_s))//".dat"
#endif

       ofile='/scratch/jharno/SimulLens_Random_mix/LOS'//trim(LOS)//'/zs_3_147_mix_231Mpc_1024/'//trim(ifile)
       !ofile='/cita/d/scratch-2week/jharno/PID/'//trim(ifile)
       !ofile=proj_path//Version//'/out/RUN-'//Run//'/'//trim(ifile)


!**********
       write(*,*) 'opening ',ofile
       open (unit=12,file=ofile,status='old',iostat=fstat,form='unformatted', access = 'stream')
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
          !write(*,*)'Reading halo', nh_total
#ifdef HALOPID
          read(12,end=112,err=113) halo_input_buffer, halo_pid(:,nh_total)
#else
          read(12,end=112,err=113) halo_input_buffer
#endif
          !write(*,*) halo_input_buffer, halo_pid(:,nh_total)
          peak(:,nh_total)=halo_input_buffer(1:3)
          xv(:,nh_total)=halo_input_buffer(4:9)
          cubepm_l(:,nh_total)=halo_input_buffer(10:12)
          cubepm_v_disp2(:,nh_total)=halo_input_buffer(13:15)
          cubepm_radius_scale(nh_total)=halo_input_buffer(16)
          halo_mass(nh_total)=halo_input_buffer(17)
          halo_mass_pp(nh_total)=halo_input_buffer(18)
          halo_mass1(nh_total)=halo_input_buffer(19)
          var(:,nh_total)=halo_input_buffer(20:22)
          I_xy(:,nh_total)=halo_input_buffer(23:28)

#ifdef verbose
          write(*,*)'peak=',peak(:,nh_total)
          write(*,*)'xv=',xv(:,nh_total)
          write(*,*)'cubepm_l=',cubepm_l(:,nh_total)
          write(*,*)'cubepm_v_disp2=',cubepm_v_disp2(:,nh_total)
          write(*,*)'cubepm_radius_scale=',cubepm_radius_scale(nh_total)
          write(*,*)'halo_mass=',halo_mass(nh_total)
          write(*,*)'halo_mass_pp=',halo_mass_pp(nh_total)
          write(*,*)'halo_mass1=',halo_mass1(nh_total)
          write(*,*)'var=',var(:,nh_total)
          write(*,*)'I_xy=',I_xy(:,nh_total)
#ifdef HALOPID
          write(*,*)'PID=',halo_pid(:,nh_total)   
#endif
#endif

!          stop

113       continue
       enddo
112    close(12)
       write(*,*) 'Cumulative Halos Read =', nh_total

#else
#ifdef pmfast
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
       xv(1,ip+1:ip+nploc(i))=modulo(xv(1,ip+1:ip+nploc(i))+node_coords(1)*ncc,rnc)
       xv(2,ip+1:ip+nploc(i))=modulo(xv(2,ip+1:ip+nploc(i))+node_coords(2)*ncc,rnc)
       xv(3,ip+1:ip+nploc(i))=modulo(xv(3,ip+1:ip+nploc(i))+node_coords(3)*ncc,rnc)
#endif
#ifdef Kaiser

    !Red Shift Distortion: x_z -> x_z +  v_z/H(Z)   
    !Converting seconds into simulation time units
    !cancels the H0...
    
    xv(3,ip+1:ip+nploc(i))=xv(3,ip+1:ip+nploc(i)) + xv(6,ip+1:ip+nploc(i))*1.5*sqrt(omegam/cubepm_a)
    
    if(i==nn) then
       write(*,*) '**********************'
       write(*,*) 'Included Kaiser Effect'
       write(*,*) 'Omega_m =', omegam, 'a =', cubepm_a
       write(*,*) '1/H(z) =', 1.5*sqrt(omegam/cubepm_a)
       write(*,*) '**********************'
    endif
#endif

#ifdef GETPID
       PIDifile=trim(zstring)//'PID'//rank_s(1:len_trim(rank_s))//'.dat'
       PIDofile=proj_path//Version//'/out/RUN-'//Run//'/'//trim(ifile)
       !PIDofile='/cita/d/scratch-2week/jharno/PID/'//trim(PIDifile)

       write(*,*) 'opening ',PIDofile

       open (unit=19,file=PIDofile,form='binary')
       read(19) nploc(i), cubepm_a,cubepm_t,cubepm_tau,cubepm_nts,cubepm_dt_f_acc,cubepm_dt_pp_acc,cubepm_dt_c_acc,cubepm_cur_checkpoint, &
            cubepm_cur_projection,cubepm_cur_halofind,cubepm_mass_p       
       read(19) PID(ip+1:ip+nploc(i))
       
       close(19)

#endif


       ip=ip+nploc(i)
       write(*,*) 'np cumulative = ', ip,', np local = ', nploc(i)
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
    !ofile=proj_path//Version//'/out/RUN-'//Run//'/'//trim(ifile)

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
    
    xv(3,:)=xv(3,:) + xv(6,:)*1.5*sqrt(omegam/cubepm_a)
    
    if(i==nn) then
       write(*,*) 'Included Kaiser Effect'
    endif
#endif



#else
    write(*,*) '***************************************************'
    write(*,*) '*** MERGE LIST NOT IMPLEMENTED FOR PARTICLES YET***'
    write(*,*) '***************************************************'
#endif
! END OF HALOS


#ifdef debug 
    write(*,*) 'np = ', np
    write(*,*) 'xv', maxval(xv(1:3,:)),minval(xv(1:3,:)),maxval(xv(4:6,:)),minval(xv(4:6,:))
#endif


#endif
! END OF MERGED LIST

    

#ifdef HALO
    np = nh_total

    if(np==0)then
       write(*,*) 'No halos found!'
       return
    else

      !**********************************
      ! Reassigning the mass of the halos
      mp = real(nc**3)/real(np)
      !**********************************

      write(*,*) 'n_halos = ', np
      write(*,*) 'halo mass = ', mp
      !write(*,*) 'total mass = ', mp*np
      !write(*,*) 'total mass pp= ', 8*hc**3
      write(*,*) 'average mass per grid cell= ', mp*np/nc**3
    endif

#endif


    write(*,*) '*************'
    write(*,*) '*Done readdm*'
    write(*,*) '*************'
    return
  end subroutine readdm


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
    write(*,*) 'Performing cic on', npt, 'haloes'
#endif

    OutBound = 0
    !! Construct chaining lists in parallel
    !$omp parallel do default(shared) private(it,ip,j,k)
    do it=1,nt
       htoc(:,:,:,it)=0
       do ip=1+(it-1)*npt,min(np,int(it*npt, kind=8))
          j=floor(xv(2,ip))+1
          k=floor(xv(3,ip))+1
          if((j > nc)) then
#ifdef debug
             write (*,*) '#### PROBLEM!!! (j = floor(xv(2,ip))+1) =',j
             write (*,*) 'Enforcing Periodic BC Manually'
#endif
             j = j-nc
             OutBound = OutBound +1
             !pause
          endif
          if((k > nc)) then
#ifdef debug
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


end program GetHaloField
