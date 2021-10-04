
program CMBLens
  use Lensing
  implicit none
#include 'Lens.fh'
 
  integer, parameter :: nk=npc/2 
  real, dimension(npc,npc) :: map,mapL,phi_fake
  real, dimension(nl) :: l,l2Cl,cl_E,cl_B,dcl,dp2
  real, dimension(nk) :: kg,delta2
  type(vec2D), dimension(npc,npc) :: defl,defl_fake
  type(Spin2), dimension(npc,npc) :: map_pol,mapL_pol
  integer i,j,is,i1,j1 
  real pi,chi,angle,dang,DA,rx,ry

  external chi
   
  pi=acos(-1.)
  DA=chi(zs,omegam,h)*h
  angle=box/DA
  dang=angle/npc

  open(10,file=dir_work//fn_CMBCl)
  do i=1,nl
     read(10,*) l(i),l2Cl(i),cl_E(i),cl_B(i)
!     read(10,*) l(i),cl_B(i),cl_E(i),l2Cl(i)
     !kg(i)=l(i)/DA 
     !l2Cl(i)=kg(i)**2
  enddo

  !! l2cl actually means cl_E
  l2cl=cl_E 

  dcl=0
  dcl(1)=(l2cl(2)-l2cl(1))/(l(2)-l(1))
  dcl(nl)=(l2cl(nl)-l2cl(nl-1))/(l(nl)-l(nl-1))
  call spline(l,l2cl,nl,dcl(1),dcl(nl),dp2)

  delta2=0
  do i=1, nk
     kg(i)=2*pi/box*i
     if(kg(i)*DA.ge.l(1).and.kg(i)*DA.le.l(nl))then
       call splint(l,l2cl,dp2,nl,kg(i)*DA,delta2(i))
     endif
  enddo

  call random_seed

!!! GaussRandomField..-ttlu,cmbmap-tjzhang

  call GaussRandomField_2d_c2c(map(:,:), box, npc, kg, delta2, nk)

  !! CMB_EE ----> Q,U, similar to kappa -----> gamma1,gamma2
  call kappa_to_defl(map(:,:),defl_fake,map_pol,phi_fake,npc) 

  call GenLens(dang,DA,defl)

  call lensmap(map,defl,dang,DA,mapL)
  call lensmap(map_pol%p,defl,dang,DA,mapL_pol%p)
  call lensmap(map_pol%c,defl,dang,DA,mapL_pol%c)

  call outputmap(map,mapL,map_pol,mapL_pol,DA)
!

end program CMBLens 
  subroutine GenLens(dang,DA,defl)
    use Lensing
    implicit none
#include 'Lens.fh'
  
    real, intent(in) :: dang,DA
    type(vec2D), intent(out), dimension(npc,npc) :: defl 
    real av1,av2,rx,aa,pi
    integer i,j

    open(10,file=dir_work//fn_defl,form='binary')
    read(10)defl%x,defl%y
    close(10) 
  
    return
  end subroutine GenLens


 subroutine lensmap(map,defl,dang,DA,nmap)
    use Lensing
   implicit none
#include 'Lens.fh'
   
   real, intent(in) :: dang, DA 
   type(vec2D), dimension(npc,npc), intent(in) :: defl 
   real, dimension(npc,npc), intent(in) :: map
   real, dimension(npc,npc), intent(out) :: nmap

   real, dimension(npc,npc) :: x,y
   integer i,j,k,hpc
   real ang_i,ang_j,rho

   hpc=npc/2
!! NOTE the bug: the position of rays on each plane have to be in the centre 
!! of grids. Remember to subtract 0.5 !!! (April 29,2007) 
   do j=1,npc
      do i=1,npc
         ang_i=(i-0.5-hpc)*dang+hpc*dang
         ang_j=(j-0.5-hpc)*dang+hpc*dang
         x(i,j)=DA*tan(ang_i-1*defl(i,j)%x)/box*npc
         y(i,j)=DA*tan(ang_j-1*defl(i,j)%y)/box*npc
         call cic2d(map(:,:),npc,x(i,j),y(i,j),rho)
         nmap(i,j)=rho         
      enddo
   enddo

   return
  end subroutine lensmap

  subroutine outputmap(map,mapL,map_pol,mapL_pol,DA)
    use Lensing
    implicit none
#include 'Lens.fh'

    real, dimension(npc,npc) :: map,mapL 
    type(vec2D), dimension(npc,npc) :: map_pol,mapL_pol
    real, dimension(npc) :: power
    integer i,j
    real pi,DA,angle,chi

    pi=acos(-1.)
    DA=chi(zs,omegam,h)*h
    angle=box/DA

    open(10,file=dir_work//fn_source,form='binary')
    write(10) map
    close(10)    
    open(20,file=dir_work//fn_image,form='binary')
    write(20) mapL
    close(20)    
    open(30,file=dir_work//fn_source_Q,form='binary')
    write(30) map_pol%p
    close(30)    
    open(40,file=dir_work//fn_source_U,form='binary')
    write(40) map_pol%c
    close(40)    
    open(30,file=dir_work//fn_image_Q,form='binary')
    write(30) mapL_pol%p
    close(30)    
    open(40,file=dir_work//fn_image_U,form='binary')
    write(40) mapL_pol%c
    close(40)    

    call shear_to_power(map_pol,power,npc)
    open(30,file=dir_work//'pol_power.dat')
    do i=1,npc/2
        write(30,*) i*2*pi/angle,power(i)
    enddo
    close(30)
  
    call shear_to_power(mapL_pol,power,npc)
    open(30,file=dir_work//'polL_power.dat')
    do i=1,npc/2
        write(30,*) i*2*pi/angle,power(i)
    enddo
    close(30)

    call ps2(map,map,power,npc)
    open(30,file=dir_work//'cmb-bfls.dat')
    do i=1,npc/2
        write(30,*) i*2*pi/angle,power(i)
    enddo
    close(30)

    call ps2(mapL,mapL,power,npc)
    open(30,file=dir_work//'cmb-afls.dat')
    do i=1,npc/2
       write(30,*) i*2*pi/angle,power(i)
    enddo
    close(30)

    return
  end subroutine outputmap 


