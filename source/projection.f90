program projection
  implicit none
#include 'Lens.fh'

  integer, parameter :: nt=1

  integer(kind=8), parameter :: np=hc**3
  real, parameter :: mp=int(nc,kind=8)**3/np !nc**3/np
  !real, dimension(6,np) :: xv
  !integer, dimension(np) :: ll
  integer, parameter :: MSL=100
 
  real, parameter  :: ncc=nc/nodes_dim !! number of cells / cubic
  real, parameter  :: rnc = nc

  !! Dark matter arrays
  integer, dimension(nc,nc) :: hoc
! integer, dimension(2,nc,nc,nt) :: htoc
! real, dimension(nc,nc,nc):: d
  real, dimension(nc,nc)::proj 
  integer(1), dimension(nc,nc)::map

  real, allocatable, dimension(:,:)        :: xv
  integer, allocatable, dimension(:)       :: ll
  integer, allocatable, dimension(:,:,:,:) :: htoc
  real, allocatable, dimension(:,:,:)      :: d



  integer nploc(nn),fstat
  real pi,current_z

  !! variables in cubepm header
  integer cubepm_nts,cubepm_cur_checkpoint,cubepm_cur_projection,cubepm_cur_halofind,i1,j1,k1, node_coords(3)
  real  cubepm_a,cubepm_t,cubepm_tau,cubepm_dt_f_acc,cubepm_dt_pp_acc,cubepm_dt_c_acc, cubepm_mass_p

!  character(len = 100) :: ifile,ofile


  pi=acos(-1.)

  print *,'starting serial density projection program'
 
  allocate(xv(6,np))
  allocate(ll(np))
  allocate(htoc(2,nc,nc,nt))
  allocate(d(nc,nc,nc))
 
  call zslice
  !call readdm
  !call cic
  !call project

  deallocate(d)
  deallocate(htoc)
  deallocate(ll)
  deallocate(xv)

contains

  subroutine zslice
    implicit none
    integer i,j
    real, dimension(nslice) :: z_write

    open(11,file=fn_z)
    do i = 1,nslice
       read(11,*)z_write(i)
       current_z = z_write(i)
       call readdm
       call cic
#ifdef debug
       do j =1,10
          write(*,*) 'd(',j,',1,1)=',d(j,1,1)
       enddo
#endif 
       if(current_z .eq. 0.025) call project
    enddo
    close(11)

  end subroutine zslice

  subroutine readdm
    implicit none
    character*100 fn, fn2, fn3

    character*8 t1,t2
    integer i,j, ii,Nmax, nh_local, nh_total, blocksize, num_writes, nplow, nphigh
    integer(kind=8) ip
    real HubbleScale
    real Conversion
    character (len=4) :: rank_s
    character (len=MSL) :: ofile,zstring
    character (len=MSL) :: ifile

    !! Read particle data file from decomposed volume
    ip=0
    nh_total=0
    do i=1,nn

       nh_local = 0

       write(*,*) 'Reading Node ',i
       write(rank_s,'(i4)') i-1
       rank_s=adjustl(rank_s)
       write(zstring,'(f5.3)') current_z !z3dps       ! Need (f6.3) for Z > 10
       ifile=trim(zstring)//'xv'//rank_s(1:len_trim(rank_s))//".dat" 
  
       if(current_z <=1.0) ofile=proj_path//'/'//trim(ifile)
       if(current_z >1.0) ofile=proj_path2//'/'//trim(ifile)
       !ofile=proj_path//'/'//trim(ifile)
       write(*,*) 'opening ',ofile
       open (unit=12,file=ofile,status='old',iostat=fstat,form='binary')
       if (fstat /= 0) then
          write(*,*) 'error opening catalog'
          write(*,*) 'rank=',nn, 'file:',ofile
          stop !call mpi_abort(mpi_comm_world,ierr,ierr)
       endif
       read(12) nploc(i), cubepm_a,cubepm_t,cubepm_tau,cubepm_nts,cubepm_dt_f_acc,cubepm_dt_pp_acc,cubepm_dt_c_acc,cubepm_cur_checkpoint, cubepm_cur_projection,cubepm_cur_halofind,cubepm_mass_p
       !write(*,*) 'Scale factor =',cubepm_a,'redshift = ',1./cubepm_a - 1

       blocksize=(32*1024*1024)/24
       num_writes=nploc(i)/blocksize+1


       do ii=1,num_writes
         nplow=(ii-1)*blocksize+1
         nphigh=min(ii*blocksize,nploc(i))
!!       print *,rank,nplow,nphigh,np_local
         do j=nplow,nphigh
           read(12) xv(:,j)
         enddo
       enddo

       !read(12) xv(1:6,ip+1:ip+nploc(i))
       close(12)

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
#ifdef debug
       do i1=1,30
          write(*,*) 'x=',xv(:3,i1)
       enddo
#endif
       ip=ip+nploc(i)
       write(*,*) 'np cumulative = ', ip,', np local = ', nploc(i)

    enddo

    write(*,*) '*************'
    write(*,*) '*Done readdm*'
    write(*,*) '*************'
    return

end subroutine readdm

!**********************************


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

    !! Construct chaining lists in parallel
    !$omp parallel do default(shared) private(it,ip,j,k)
    OutBound = 0
    do it=1,nt
       htoc(:,:,:,it)=0
       do ip=1+(it-1)*npt,min(np,it*npt)
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

subroutine project
  implicit none

  integer i,j,k
  character (len=MSL) :: zstring
  character (len=MSL) :: ifile

  write(zstring,'(f5.3)') current_z!z3dps       ! Need (f6.3) for Z > 10
  ifile=trim(zstring)
  !ofile=proj_path//'/'//trim(ifile)
  !ofile=trim(ofile) 
  ! x-y
#ifdef debug
  write(*,*)'max d = ', maxval(d), 'min d = ',minval(d)
  do i = 1,10
     write(*,*) 'd(',i,',1,1)=',d(i,1,1)
  enddo
#endif

  proj=0.0
  do k=1,nc!/nodes_dim
     do j = 1,nc
        do i = 1,nc  
           proj(i,j)=proj(i,j)+d(i,j,k)
#ifdef debug
           write(*,*) i,j,k,d(i,j,k),proj(i,j)
           pause
#endif
        enddo
     enddo
  enddo
  open(10,file=proj_path//'/'//trim(ifile)//'proj_xy_test.dat',status='replace',form='binary')
  write(10) 1.0/(current_z+1.0)
  write(10)  proj
  close(10)

#ifdef debug
  write(*,*)'proj max value=',maxval(proj)
  write(*,*)'proj min value=',minval(proj)
#endif

  proj=proj-minval(proj)
  proj=proj/maxval(proj)
  map=127*sqrt(sqrt(proj))

#ifdef debug
  write(*,*)'map max value=',maxval(map)
  write(*,*)'map min value=',minval(map)
  do k = 1,40
    write(*,*) 'diagonal of Map = ',map(k,k)
  enddo
#endif

  open(10,file=proj_path//'/'//trim(ifile)//'proj_xy_test.pgm',status='replace')
  write(10,'(2hP5)')
  write(10,*) nc,nc
  write(10,*) 127
  close(10)
  open(10,file=proj_path//'/'//trim(ifile)//'proj_xy_test.pgm',access='append',form='binary')
  write(10) map
  close(10)
  write(*,*) 'Wrote' , proj_path,'/',trim(ifile),'proj_xy_test.pgm'

  ! x-z
  proj=0.0
  do k = 1,nc
     do j=1,nc!/nodes_dim
        do i = 1,nc
           proj(i,k)=proj(i,k)+d(i,j,k)
        enddo
     enddo
  enddo
  open(10,file=proj_path//'/'//trim(ifile)//'proj_xz_test.dat',status='replace',form='binary')
  write(10) 1.0/(current_z+1.0)
  write(10)  proj
  close(10)
  proj=proj-minval(proj)
  proj=proj/maxval(proj)
  map=127*sqrt(sqrt(proj))

  open(10,file=proj_path//'/'//trim(ifile)//'proj_xz_test.pgm',status='replace')
  write(10,'(2hP5)')
  write(10,*) nc,nc
  write(10,*) 127
  close(10)
  open(10,file=proj_path//'/'//trim(ifile)//'proj_xz_test.pgm',access='append',form='binary')
  write(10) map
  close(10)
  write(*,*) 'Wrote' , proj_path,'/',trim(ifile),'proj_xz.pgm'


  ! y-z

  proj=0.0
  do k=1,nc
     do j = 1,nc
        do i=1,nc!/nodes_dim
           proj(j,k)=proj(j,k)+d(i,j,k)
        enddo
     enddo
  enddo
  open(10,file=proj_path//'/'//trim(ifile)//'proj_yz_test.dat',status='replace',form='binary')
  write(10) 1.0/(current_z+1.0)
  write(10)  proj
  close(10)
  proj=proj-minval(proj)
  proj=proj/maxval(proj)
  map=127*sqrt(sqrt(proj))

  open(10,file=proj_path//'/'//trim(ifile)//'proj_yz_test.pgm',status='replace')
  write(10,'(2hP5)')
  write(10,*) nc,nc
  write(10,*) 127
  close(10)
  open(10,file=proj_path//'/'//trim(ifile)//'proj_yz_test.pgm',access='append',form='binary')
  write(10) map
  close(10)
  write(*,*) 'Wrote' , proj_path,'/',trim(ifile),'proj_yz.pgm'


end subroutine project

!****************************

end program projection


 

 
