
!! modified on Hy Trac's code.
program cicpow
  implicit none
#include 'par.fh'

  !! np should be set to nc (1:1) or hc (1:2)
  integer, parameter :: nt=1
  integer(kind=8), parameter :: np=hc**3
  real, parameter :: mp=int(nc,kind=8)**3/np !nc**3/np

  !! cubepm
!  integer, parameter :: nodes_dim=2
  real, parameter  :: ncc=nc/nodes_dim !! number of cells / cubic 
  real, parameter  :: rnc = nc

  !! Dark matter arrays
  real, dimension(6,np) :: xv
  integer, dimension(np) :: ll
  integer, dimension(nc,nc) :: hoc
  integer, dimension(2,nc,nc,nt) :: htoc

  !! Power spectrum arrays
  real, dimension(2,nc) :: ps
  real, dimension(nc+2,nc,nc) :: d 

  integer, parameter :: MSL=100
  
  integer nploc(nn)
  real pi

  !! variables in CubePM
  integer cubepm_nts,cubepm_cur_checkpoint,cubepm_cur_projection,cubepm_cur_halofind,i1,j1,k1, node_coords(3)
  real  cubepm_a,cubepm_t,cubepm_tau,cubepm_dt_f_acc,cubepm_dt_pp_acc,cubepm_dt_c_acc, cubepm_mass_p

  common /rarr/ d,xv,ps 
  common /iarr/ ll,htoc,hoc

  pi=acos(-1.)

  call readdm 
  call cic
  call powerspectrum
  call writeps

contains

  subroutine readdm
    implicit none
    character*100 fn

    character*8 t1,t2
    integer i,j 
    integer(kind=8) ip
    character (len=4) :: rank_s
    character (len=MSL) :: ofile,zstring
    character (len=MSL) :: ifile

    !! Read particle data file
    ip=0
    do i=1,nn
       write(*,*) i
       write(rank_s,'(i4)') i-1
       rank_s=adjustl(rank_s)
       write(zstring,'(f5.3)') z3dps
       ifile=trim(zstring)//'xv'//rank_s(1:len_trim(rank_s))//".dat"
       ofile=PROJ_PATH//trim(ifile)
       write(*,*) ofile
       open (unit=12,file=ofile,form='binary')
#ifdef pmfast
       read(12) nploc(i)
#endif
#ifdef cubepm
       read(12) nploc(i), cubepm_a,cubepm_t,cubepm_tau,cubepm_nts,cubepm_dt_f_acc,cubepm_dt_pp_acc,cubepm_dt_c_acc,cubepm_cur_checkpoint, cubepm_cur_projection,cubepm_cur_halofind,cubepm_mass_p
#endif
       read(12) xv(1:6,ip+1:ip+nploc(i))
       close(12)
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
       ip=ip+nploc(i)
       write(*,*) ip,nploc(i)
    enddo
 
    write(*,*) 'Read ',fn
    return
  end subroutine readdm

  subroutine writeps
    implicit none
    integer k
    character*100 fn

    !! Output power spectrum
    !! First column is k
    !! Second column is \Delta^2

    fn=dir_work//'3d_z_0.200'
    open(11,file=fn)
    do k=1,nc !hc+1
       write(11,*) 2*pi/lbox*k,ps(2,k) !, exp(-2*ps(1,k)**2/sig**2)/box**3*ps(1,k)**3/2/pi**2*AA**2
    enddo
    close(11)
   
    write(*,*) 'Called writeps'
    return
  end subroutine writeps


  subroutine cic
    implicit none
    integer, parameter :: kpt=nc/nt
    integer, parameter :: npt=np/nt

    integer it,i,j,k,ip
    real toe

    !! Construct chaining lists in parallel
    !$omp parallel do default(shared) private(it,ip,j,k)
    do it=1,nt
       htoc(:,:,:,it)=0
       do ip=1+(it-1)*npt,min(np,it*npt)
          j=floor(xv(2,ip))+1
          k=floor(xv(3,ip))+1
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

    write(*,*) 'Called cic'
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

    real pst(2,nc,nt)
    integer i,j,k

    !! Convert density to overdensity
    do k=1,nc
       do j=1,nc
          do i=1,nc
             d(i,j,k)=d(i,j,k)-1
          enddo
       enddo
    enddo

    call ps3_r2c(d,ps(2,:),nc)

    write(*,*) 'Called power spectrum'
    return
  end subroutine powerspectrum

!!--------------------------------------------------------------!!

end program cicpow
