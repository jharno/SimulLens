!! Angular Covariance written by Joachim Harnois- Deraps
!! Power-spectrum modified on Hy Trac's + Ting Ting Lu's code.
program cicpow
  implicit none
#include 'par.fh'

  !! np should be set to nc (1:1) or hc (1:2)
  integer, parameter :: nt=1

#ifndef HALO
  integer(kind=8), parameter :: np=hc**3
  real, parameter :: mp=int(nc,kind=8)**3/np !nc**3/np
  !real, dimension(6,np) :: xv
  !integer, dimension(np) :: ll

#ifdef GETPID
  integer(kind=8), dimension(np) :: PID
#endif

#else
  integer(kind=8), parameter :: np_max = 1000000!hc**3 
  integer(kind=8) :: np  
  real :: mp
  real, dimension(6,np_max) :: xv
  integer, dimension(np_max) :: ll
#ifdef HALOPID  
  integer,parameter :: N_p = 10
  integer(kind=8), dimension(N_p,100000) :: halo_pid
#endif
#endif



  !! cubepm
!  integer, parameter :: nodes_dim=2
  real, parameter  :: ncc=nc/nodes_dim !! number of cells / cubic 
  real, parameter  :: rnc = nc

  !! Dark matter arrays
  !integer, dimension(nc,nc) :: hoc
  !integer, dimension(2,nc,nc,nt) :: htoc

  !! Power spectrum arrays
  !real, dimension(nc) :: ps,err
  !real, dimension(nc) :: PS_Ave!,  CrossCorrAve, N_CrossCorrAve
  !real, dimension(3*(nc/2)**2+1) :: CrossCorrAve, N_CrossCorrAve
  !real, dimension(nc+2,nc,nc) :: d, CrossCorr, N_CrossCorr
  !real, dimension(nc+2,nc,nc) ::  psShell_Inner, psShell_Outer, N_psShell_Inner, N_psShell_Outer
  !complex, dimension(nc/2 + 1,nc,nc) ::  d_COMPLEX, psV_COMPLEX, CrossCorr_COMPLEX, N_CrossCorr_COMPLEX 
  !complex, dimension(nc/2 + 1,nc,nc) ::  psShell_Inner_COMPLEX, psShell_Outer_COMPLEX, N_psShell_Inner_COMPLEX,N_psShell_Outer_COMPLEX
  
  !equivalence (psShell_Inner,psShell_Inner_COMPLEX,CrossCorr, CrossCorr_COMPLEX)
  !equivalence (psShell_Outer,psShell_Outer_COMPLEX)
  !equivalence (N_psShell_Inner,N_psShell_Inner_COMPLEX)
  !equivalence (d, d_COMPLEX, psV_COMPLEX, N_CrossCorr,N_psShell_Outer,N_psShell_Outer_COMPLEX, N_CrossCorr_COMPLEX) 
  
 
#ifdef FFT_TEST
 real, dimension(nc+2,nc,nc) :: Sync, Gauss, DeltaFunction, N_Hat, Hat
 real, dimension(nc) :: Xsi, HatAve, GaussAve, DeltaAve, FFT_Hat
#endif

  !real, dimension(0:hc,-hc+1:hc,-hc+1:hc) :: Hat
  integer, parameter :: MSL=100
  character (len=MSL) :: ofile

  integer nploc(nn),fstat
  real pi,omegak

  !! variables in CubePM
  integer cubepm_nts,cubepm_cur_checkpoint,cubepm_cur_projection,cubepm_cur_halofind,i1,j1,k1, node_coords(3)
  real  cubepm_a,cubepm_t,cubepm_tau,cubepm_dt_f_acc,cubepm_dt_pp_acc,cubepm_dt_c_acc, cubepm_mass_p

#ifdef HALO
  real, dimension(np_max) :: cubepm_v_disp, cubepm_radius_scale, halo_mass, halo_mass_pp, halo_mass1
  real, dimension(3,np_max) :: cubepm_halo_pos, peak, cubepm_l, var
#endif

  !common /rarr/ d,xv,ps 
  !common /iarr/ ll,htoc,hoc

  pi=acos(-1.)

  call readdm 
  call PID_analysis



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
       write(rank_s,'(i4)') i-1!+8 too allow for a shift in nodes sampling
       rank_s=adjustl(rank_s)
       if(z3dps>=10.0)then 
          write(zstring,'(f6.3)') z3dps       ! Need (f6.3) for Z > 10
       else
          write(zstring,'(f5.3)') z3dps
       endif
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
       ofile='/mnt/scratch-3week/jharno/PID2/'//trim(ifile)
       !ofile=proj_path//Version//'/out/RUN-'//Run//'/'//trim(ifile)


!**********
       

#ifdef HALO       
       write(*,*) 'Opening ',ofile
       open (unit=12,file=ofile,status='old',iostat=fstat,form='binary')
       if (fstat /= 0) then
          write(*,*) 'error opening catalog'
          write(*,*) 'rank=',nn, 'file:',ofile
          stop !call mpi_abort(mpi_comm_world,ierr,ierr)
       endif
#else       
       write(*,*) '**Not** opening ',ofile
#endif

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
       !read(12) nploc(i)
#endif
#ifdef cubepm
       !read(12) nploc(i), cubepm_a,cubepm_t,cubepm_tau,cubepm_nts,cubepm_dt_f_acc,cubepm_dt_pp_acc,cubepm_dt_c_acc,cubepm_cur_checkpoint, cubepm_cur_projection,cubepm_cur_halofind,cubepm_mass_p
#endif
       !read(12) xv(1:6,ip+1:ip+nploc(i))
       !close(12)
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
       !xv(1,ip+1:ip+nploc(i))=modulo(xv(1,ip+1:ip+nploc(i))+node_coords(1)*ncc,rnc)
       !xv(2,ip+1:ip+nploc(i))=modulo(xv(2,ip+1:ip+nploc(i))+node_coords(2)*ncc,rnc)
       !xv(3,ip+1:ip+nploc(i))=modulo(xv(3,ip+1:ip+nploc(i))+node_coords(3)*ncc,rnc)
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
#ifdef IC
       PIDifile='PID'//rank_s(1:len_trim(rank_s))//'.ic'
#else
       PIDifile=trim(zstring)//'PID'//rank_s(1:len_trim(rank_s))//'.dat'
#endif
       !PIDofile=proj_path//Version//'/out/RUN-'//Run//'/'//trim(ifile)
       PIDofile='/mnt/scratch-3week/jharno/PID2/'//trim(PIDifile)

       write(*,*) 'opening ',PIDofile

       open (unit=19,file=PIDofile,form='binary')
#ifdef IC
       read(19) nploc(i)
#else
       read(19) nploc(i), cubepm_a,cubepm_t,cubepm_tau,cubepm_nts,cubepm_dt_f_acc,cubepm_dt_pp_acc,cubepm_dt_c_acc,cubepm_cur_checkpoint, &
            cubepm_cur_projection,cubepm_cur_halofind,cubepm_mass_p       
#endif
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
    !mp = real(nc**3)/real(np)
    ! might have to cancel large numbers first
    mp = (real(nc)/(real(np**(1.0/3.0))))**3

    !**********************************
   
    write(*,*) 'halo mass = ', mp
    write(*,*) 'total mass = ', mp*np
    write(*,*) 'total mass pp= ', 8*hc**3
    !write(*,*) 'average mass per grid cell= ', mp*np/nc**3
    write(*,*) 'average mass per grid cell= ', mp*(real(np**(1.0/3.0))/real(nc))**3

#endif


    write(*,*) '*************'
    write(*,*) '*Done readdm*'
    write(*,*) '*************'
    return
  end subroutine readdm


!*******************************
subroutine PID_analysis
  implicit none

  integer(kind=8) i,j,k,l
  integer(kind=8) PID1,PID2,test_PID
  integer(kind=8), dimension(np):: isortPID
  integer(kind=8), dimension(10):: isort_test
  integer n_rep, n_zeros
  real(4),dimension(np)  ::       distance

#ifdef GETPID

  !isort_test(:)=(/ (i,i=1,10)  /)
  !write(*,*) 'Sorting test'
  !write(*,*) 'PID before', PID(10:1:-1)
  !write(*,*) 'isort array before',isort_test(:)

  !call indexedsort_PID(10,PID(10:1:-1),isort_test(:))

  !write(*,*) 'sorting PIDs'
  !write(*,*) 'PID after', PID(10:1:-1)
  !write(*,*) 'isort array after',isort_test(:)
  !pause

  !write(*,*)'np = ', np
  !write(*,*)'PID begin before:' ,PID(1:10)
  !write(*,*)'PID end   before:' ,PID(np-10:np)
  isortPID(:)=(/ (i,i=1,np) /)
  !write(*,*) 'isortPID(1:10) = ', isortPID(1:10)
  write(*,*) 'sorting PIDs'
  write(*,*) 'PID(np/8)=', PID(np/8:np/8+100)
  write(*,*) 'PID(np/4)=', PID(np/4:np/4+100)
  write(*,*) 'PID(np/2)=', PID(np/2:np/2+100)
  write(*,*) 'PID(3*np/4)=', PID(3*np/4:3*np/4+100)
  write(*,*) 'PID(np)=', PID(np-100:np)

  call indexedsort_PID(np,PID(:),isortPID(:))
  !PID(:) = PID(isortPID(:))
  !write(*,*) 'isortPID(1:10) = ', isortPID(1:10)
  write(*,*) '*************'
  write(*,*) 'sorted PIDs:'
  write(*,*) '*************'
  write(*,*) 'PID(np/8)=', PID(np/8:np/8+100)
  write(*,*) 'PID(np/4)=', PID(np/4:np/4+100)
  write(*,*) 'PID(np/2)=', PID(np/2:np/2+100)
  write(*,*) 'PID(3*np/4)=', PID(3*np/4:3*np/4+100)
  write(*,*) 'PID(np)=', PID(np-100:np)


  !stop

  n_rep = 0
  n_zeros = 0
  do i=1,np-1       
     PID1 = PID(i)
     PID2 = PID(i+1)
     if(PID1 .eq. PID2 .and. PID1 .ne. 0) then
        !write(*,*) 'Identical PID found:',i,PID1,PID2
        n_rep = n_rep + 1
     endif
     if(PID1 .eq. 0)then
        n_zeros = n_zeros + 1
     endif
  enddo
  write(*,*) 'N repetitions = ', n_rep
  write(*,*) 'N zeros = ', n_zeros
  write(*,*) 'N total =', np

#endif

#ifdef HALOPID

  !*** look for identical PIDs in the same halo
  n_rep = 0

  do i = 1,np
    !write(*,*) halo_pid(:,i) 
    do j = 1,N_p-1
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

  write(*,*) 'Done looking for repeated PIDs in same halos'

  !*** look for identical PIDs all halos

  do i = 1,np-1
    !write(*,*) halo_pid(:,i) 
    do j = 1,N_p

       do k = i+1,np
          do l = 1,N_p
             PID1 = halo_pid(j,i)
             PID2 = halo_pid(l,k)
             !write(*,*) i,j,k,l
             !pause
             if((PID1 .eq. PID2) .and. PID1 .ne. 0 ) then

                n_rep = n_rep+1

                write(*,*)'Halo PIDs 1:',halo_pid(:,i)
                write(*,*)'Halo PIDs 2:',halo_pid(:,k)
                write(*,*)'Found identical PIDs:' , PID1,PID2
                     

                distance(n_rep) =  sqrt((xv(1,i)-xv(1,k))**2 + (xv(2,i)-xv(2,k))**2 + (xv(3,i)-xv(3,k))**2)
                
             endif
          enddo
       enddo
    enddo
  enddo

  write(*,*) distance(1:n_rep)

#endif

  write(*,*) '*******************'
  write(*,*) '*done PID_analysis*'
  write(*,*) '*******************'

  return
end subroutine PID_analysis
!*******************************
!*******************************



end program cicpow
