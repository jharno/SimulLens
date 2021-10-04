! compile with: ifort ...
! Estimators.f90 written by Joachim Harnois-Deraps 08/04/2013
! Program to compute accurate weak lensing estimators from 2D midplanes output from N-body simulations.

program Estimators 
  use Lensing
  implicit none
  include 'Lens.fh'

interface 
   function chi (z,omegam,h)
      real:: chi
      real, intent(in):: z,omegam,h
   end function chi

end interface

  ! Which analyses step to do?
  
  ! Compute l2Cl on peridic slabs
  logical, parameter :: z_slice_2pt = .true.
  logical, parameter :: z_integrate_2pt = .false.
  logical, parameter :: GaussPower_2pt = .false.


  !*************
  !1- initialize variables

  ! Arrays on original grid
  real, dimension(nc,nc) :: input_map,x1,x2 
  complex, dimension(nc,nc):: map_cplx
  !real, dimension(2*nc,2*nc):: zoom_map

!#ifdef calshear
!  real, dimension(nc,nc) :: phi,tmp_map 
!  type(Spin2), dimension(nc,nc) :: shear  
!  type(vec2D), dimension(nc,nc) :: defl 
!#endif


!  ! Arrays on pixels
!  real, dimension(npc,npc,nslice) :: map_3D
!  real, dimension(npc,npc) :: cumul_kappa
!#ifdef calshear
!  real, dimension(npc,npc) :: cumul_gamma1,cumul_gamma2,cumul_deflx,cumul_defly,tmp_pix_map 
!  type(vec2D), dimension(npc,npc,nslice) :: newdefl
!  type(vec2D), dimension(npc,npc) :: CorrBornDefl
!  type(Spin2), dimension(npc,npc,nslice) :: newshear
!#endif

  ! Other stuff
  real, dimension(nc) :: power
  real, dimension(nslice) :: z_write,z_write_s
!  type(vec2D) shift
  integer cur_z,i,icount, file_status,fu
  !integer j,i1,j1,j2,i3,kx,ky,ir,index_nbody_run, newLOS,patch!,redshift, z
  real lense_weight,frac,ncr,angle, angle_slice, dang,pi,random_index_run,box,kernel ,ScaleFactor
  real(kind=8) rhomean
  character(len=180) :: fn,fp! ,fn1,my_status, test_str
  character(len=7) z_string, LOS_str!, index_str, newLOS_str

  !-------- Multi-plateform code:
  ! in/output_form =  TCS: 'unformatted' ; GPC: 'binary'
  ! output_access  =  TCS: 'stream'      ; GPC: 'sequential'

  ! GPC:
  character(*), parameter :: output_form=   'binary'
  character(*), parameter :: output_access= 'sequential' 
  character(*), parameter :: input_form=    'binary' 

  ! TCS:
  !character(*), parameter :: output_form=   'unformatted'
  !character(*), parameter :: output_access= 'stream' 
  !character(*), parameter :: input_form=    'unformatted' 

  box = lbox
  ncr = nc
  pi=acos(-1.)

  angle = sqrt(Area)/180.*pi ! in rad

  Write(*,*) '******************'
  write(*,*) 'Running Estimators'
  write(*,*) 'Maps size     = ', nc
  write(*,*) 'Field of view = ', (angle*180/pi)**2, 'sq.degrees'
  write(*,*) 'N slices      = ', nslice 

  !-----------
  !read z-lens
  open(11,file=fn_z)  
  do cur_z=1,nslice
     read(11,*) z_write(cur_z)
  enddo
  close(11)
  write(*,*)'z_lens = ', z_write

  !read z_sources
  open(11,file=fn_z_s)  
  do cur_z=1,nslice
     read(11,*) z_write_s(cur_z)
  enddo
  close(11)
  !write(*,*)'z_source = ', z_write_s
  !------------

  call GetArg(1,LOS_str)

  !**************
  if(z_slice_2pt)     call Get_z_slice_2pt
  if(z_integrate_2pt) call Get_z_integrate_2pt
  if(GaussPower_2pt)      call Get_GaussPower           
  !**************

contains 

!-----------------------------
subroutine Get_z_integrate_2pt
  implicit none

  real ::  dCl_1, dCl_hc
  real, dimension(hc) ::  ell_grid, d2ydx2
  real, dimension(npc/2) :: ell_pixel
  real, dimension(npc/2,nslice) :: Cl_dm, Cl_zs
  integer :: cur_z_s

  write(*,*) '**************************'
  write(*,*) 'Summing l2Cl on the pixels'
  write(*,*) '**************************'


  !1- Get ell_pixel
  ! ell + 1 = 2*pi/theta_rad
  do i=1,npc/2
     ell_pixel(i) = (2*pi/angle-1)*i
  enddo

  do cur_z=1,nslice

     !2- load l2Cl(z) (convert to Cl?)
     write(z_string,'(f7.3)') z_write(cur_z)
     z_string=adjustl(z_string)
     open(30,file=Lens_output_path//'LOS'//trim(LOS_str)//'_'//trim(z_string)//'l2Cl_slice.dat')
     do i=1,hc
        read(30,*) ell_grid(i),power(i)
        !power(i)=2*pi*power(i)/ell_grid(i)**2
     enddo
     close(30)     

     !----------------------------------------------
     !3- interpolate Cl(ell_grid) onto Cl(ell_pixel)
     dCl_1  = (power(2) - power(1))/(ell_grid(2) - ell_grid(1))
     dCl_hc = (power(hc) - power(hc-1))/(ell_grid(hc) - ell_grid(hc-1))
     call spline(ell_grid, power, hc, dCl_1,dCl_hc,d2ydx2)
     
     do i=1,npc/2
       call splint(ell_grid, power, d2ydx2, hc, ell_pixel(i), Cl_dm(i,cur_z))
     enddo      

  enddo


  Cl_zs = 0;
  do cur_z_s=1,nslice
     do cur_z=1,cur_z_s

        !4- Get lensing kernel at z_lens and weight
        kernel = (3./2)*omegam*(box/nc)/(3.E3)**2*(chi(z_write(cur_z),omegam,h)*h)*(1+z_write(cur_z))* &
                   (1-chi(z_write(cur_z),omegam,h)/chi(z_write_s(cur_z_s),omegam,h))
     
        Cl_zs(:,cur_z_s) = Cl_zs(:,cur_z_s) + Cl_dm(:,cur_z)*kernel**2

     enddo
  enddo 

  !5- Write [ell_pixel Cl_zs] to file
  open(30,file=Lens_output_path//'LOS'//trim(LOS_str)//'_l2Cl_Allsource.dat')
  do i=1,npc/2
     write(30,"(25e15.5)") ell_pixel(i),Cl_zs(i,:)
  enddo
  close(30)

  write(*,*) 'Wrote', Lens_output_path//'LOS'//trim(LOS_str)//'_l2Cl_Allsource.dat'

  return
end subroutine Get_z_integrate_2pt
!--------------------------------
subroutine Get_z_slice_2pt
  implicit none

  write(*,*) '**************************'
  write(*,*) 'Computing l2Cl on 2D slabs'
  write(*,*) '**************************'

  icount=1
  !icount=0
  do cur_z=1,nslice

     !1- load projections

     i=mod(icount,3)

     ! Get the Redshift
     write(z_string,'(f7.3)') z_write(cur_z)
     z_string=adjustl(z_string)
 
     write(*,*) '******************'
     write(*,*) 'Processing z = ',  z_string

     fn=proj_path//z_string(1:len_trim(z_string))
     fn=adjustl(fn)

     if(i.eq.0)fp=fn(1:len_trim(fn))//'proj_half_finer_xy.dat_LOS'//trim(LOS_str)
     if(i.eq.1)fp=fn(1:len_trim(fn))//'proj_half_finer_yz.dat_LOS'//trim(LOS_str)
     if(i.eq.2)fp=fn(1:len_trim(fn))//'proj_half_finer_xz.dat_LOS'//trim(LOS_str)

     fu=10+i
     open(unit=fu,file=fp, form=input_form,status='old',iostat=file_status)
     if(file_status.gt.0) then
        close(fu)
        write(*,*) 'Failed to open file',fp
        write(*,*) 'file_status=',file_status
     else
        write(*,*) 'Opened',fp
     endif

     read(fu) ScaleFactor
     write(*,*) 'Read z = ', 1.0/ScaleFactor - 1.0
     read(fu) input_map
     write(*,*) 'Read proj'
     close(fu)
     !write(*,*) 'Got projection!'

     !2- compute l2Cl(z_slice)    
  
     !Correct for wrong mean in finer projections
     input_map = input_map*64.
 
     ! turn into overdensity
     rhomean = sum(real(input_map,kind=8))/nc/nc
     write(*,*) 'Mean before subtraction = ', rhomean
     input_map = input_map - rhomean

     power = 0.0
     call ps2_deconvolve(input_map,input_map,power,nc)
     !call ps2(input_map,input_map,power,nc)
     write(*,*)'Done l2Cl'

     angle_slice =atan( box/2./chi(z_write(cur_z),omegam,h)/h )*2.0

     !3- write l2Cl(z) estimators to file for each z     
     open(30,file=Lens_output_path2//'LOS'//trim(LOS_str)//'_'//trim(z_string)//'l2Cl_slice_new.dat')
     do i=1,nc/2
        write(30,*) i*(2*pi/angle_slice -1 ),power(i) 
     enddo
     close(30)
     icount=icount+1
     write(*,*) 'wrote', Lens_output_path2//'LOS'//trim(LOS_str)//'_'//trim(z_string)//'l2Cl_slice_new.dat' 

  enddo
  return
end subroutine Get_z_slice_2pt
!-----------------------------

subroutine Get_GaussPower
  implicit none

  real, dimension(100000) :: kgauss
  real, dimension(100000) :: ps_gauss
  !character(len=MSL) :: psGaussFile

  !real, dimension(nc,nc) :: x1,x2 

  write(*,*) '**********************'
  write(*,*) ' Computing Gauss Power'
  write(*,*) '**********************'

  z_string = '0.528'
  cur_z = 13

    !open(30,file=Lens_output_path//'Cl_theory.dat')
    !open(30,file=Lens_output_path//'LOS'//trim(LOS_str)//'_'//trim(z_string)//'l2Cl_slice_deconvolved.dat')
    open(30,file='./elCl_par_WMAP9_SN_BAO_z0_z0.582_nokcut.dat')
    do i=1,100000
       read(30,*) kgauss(i),ps_gauss(i)
    enddo
    close(30)
    !do i = 1,nc
    !   kgauss(i) = i*2*pi/lbox
    !enddo

!#ifdef debug_gauss
!    write(*,*)'Filled kgauss:',kgauss,ps_gauss
!#endif

    call GaussRandomField_2d_r2c(input_map, lbox, nc, kgauss, ps_gauss, 100000,x1,x2)

    write(*,*) 'Got Random Field:', minval(input_map), maxval(input_map), sum(input_map)

    !call ps3_r2c(d,ps_gauss,nc)

    power = 0.0
    call ps2_deconvolve(input_map,input_map,power,nc)
    !call ps2(input_map,input_map,power,nc)
    write(*,*)'Done l2Cl'

    angle_slice =atan( box/2./chi(z_write(cur_z),omegam,h)/h )*2.0

    !3- write l2Cl(z) estimators to file for each z     
    open(30,file=Lens_output_path//'LOS'//trim(LOS_str)//'_'//trim(z_string)//'l2Cl_slice_Gauss.dat')
    do i=1,nc/2
       write(30,*) i*2*pi/angle_slice,power(i)
    enddo
    close(30)
    icount=icount+1
    write(*,*) 'wrote', Lens_output_path//'LOS'//trim(LOS_str)//'_'//trim(z_string)//'l2Cl_slice_Gauss.dat' 

  return
end subroutine Get_GaussPower

!*************


!endif

end program Estimators
