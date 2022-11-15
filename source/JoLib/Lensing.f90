module Lensing

  type :: vec2D
       real x
       real y
  end type vec2D
  
  type :: vec2D_r8
       real*8 x
       real*8 y
  end type vec2D_r8

  type :: Spin2 
       real p !! +
       real c !! x 
  end type Spin2  

contains

 subroutine kappa_to_shear(kappa,defl,shear,phi,fmap,nc)
    implicit none

    integer nc
    !real, intent(in), dimension(nc,nc) :: kappa
    real, dimension(nc,nc) :: kappa
    real,dimension(nc,nc) :: phi
    !type(vec2D_r8), dimension(nc,nc) :: defl
    type(vec2D), dimension(nc,nc) :: defl
    type(Spin2), dimension(nc,nc) :: shear 

    ! should perform these locally?
    complex, dimension(nc,nc) :: fmap
    integer kx,ky,i1,j1
    real pi

    pi=acos(-1.)


    !--------
    ! Get phi

    fmap=kappa(:,:)+(0.0,0.0)

    write(*,*)' kappa' , sum(kappa(:,:)/nc/nc), minval(kappa(:,:)),maxval(kappa(:,:)) 
    write(*,*)' fmap (real) before fft' , sum(real(fmap(:,:))/nc/nc), minval(real(fmap(:,:))),maxval(real(fmap(:,:))) 
    write(*,*)' fmap (imag) before fft' , sum(aimag(fmap(:,:))/nc/nc), minval(aimag(fmap(:,:))),maxval(aimag(fmap(:,:))) 

    call sfft2_c2c(fmap,nc,1)

    write(*,*)'fmap (real) after fft' , sum(real(fmap(:,:))/nc/nc), minval(real(fmap(:,:))),maxval(real(fmap(:,:))) 
    write(*,*)'fmap (imag) after fft' , sum(aimag(fmap(:,:))/nc/nc), minval(aimag(fmap(:,:))),maxval(aimag(fmap(:,:))) 

    do kx=0,nc-1
       do ky=0,nc-1
          if((kx.eq.0).and.(ky.eq.0)) then
             fmap(1,1) = (0.0,0.0)
          else
             fmap(kx+1,ky+1)=-1.0/2*fmap(kx+1,ky+1)*(1**2)&
               &/((sin(kx*pi*1./nc))**2+(sin(ky*pi*1./nc))**2)
             !fgama1(kx+1,ky+1)=fphi(kx+1,ky+1)*(-2.0)/(1**2)&
             !  &*((sin(kx*pi*1./nc))**2-(sin(ky*pi*1./nc))**2)
             !fgama2(kx+1,ky+1)=fphi(kx+1,ky+1)/(1**2)&
             !  &*(1+(cos(2*pi*kx*1./nc)+(0.0,1.0)*sin(2*pi*kx*1./nc))*&
             !  &(cos(2*pi*ky*1./nc)+(0.0,1.0)*sin(2*pi*ky*1./nc))-&
             !  &(cos(2*pi*kx*1./nc)+(0.0,1.0)*sin(2*pi*kx*1./nc))-&
             !  &(cos(2*pi*ky*1./nc)+(0.0,1.0)*sin(2*pi*ky*1./nc)))
          endif
       enddo
    enddo

    call sfft2_c2c(fmap,nc,-1)
    phi=fmap

    !---------
    ! Get defl

    do j1=1,nc
       do i1=1,nc
          defl(i1,j1)%x=(phi(modulo(i1,nc)+1,j1)-phi(modulo(i1-2,nc)+1,j1))/2
          defl(i1,j1)%y=(phi(i1,modulo(j1,nc)+1)-phi(i1,modulo(j1-2,nc)+1))/2
       enddo
    enddo
    write(*,*) maxval(defl(:,:)%x),minval(defl(:,:)%y)


    !-----------
    ! Get gamma1

    fmap=kappa(:,:)+(0.0,0.0)

    write(*,*)' kappa' , sum(kappa(:,:)/nc/nc), minval(kappa(:,:)),maxval(kappa(:,:)) 
    write(*,*)' fmap (real) before fft' , sum(real(fmap(:,:))/nc/nc), minval(real(fmap(:,:))),maxval(real(fmap(:,:))) 
    write(*,*)' fmap (imag) before fft' , sum(aimag(fmap(:,:))/nc/nc), minval(aimag(fmap(:,:))),maxval(aimag(fmap(:,:))) 

    call sfft2_c2c(fmap,nc,1)

    write(*,*)'fmap (real) after fft' , sum(real(fmap(:,:))/nc/nc), minval(real(fmap(:,:))),maxval(real(fmap(:,:))) 
    write(*,*)'fmap (imag) after fft' , sum(aimag(fmap(:,:))/nc/nc), minval(aimag(fmap(:,:))),maxval(aimag(fmap(:,:))) 

    do kx=0,nc-1
       do ky=0,nc-1
          if((kx.eq.0).and.(ky.eq.0)) then
             fmap(1,1) = (0.0,0.0)
          else
             fmap(kx+1,ky+1)=-1.0/2*fmap(kx+1,ky+1)*(1**2)&
               &/((sin(kx*pi*1./nc))**2+(sin(ky*pi*1./nc))**2)
             fmap(kx+1,ky+1)=fmap(kx+1,ky+1)*(-2.0)/(1**2)&
               &*((sin(kx*pi*1./nc))**2-(sin(ky*pi*1./nc))**2)
             !fgama2(kx+1,ky+1)=fphi(kx+1,ky+1)/(1**2)&
             !  &*(1+(cos(2*pi*kx*1./nc)+(0.0,1.0)*sin(2*pi*kx*1./nc))*&
             !  &(cos(2*pi*ky*1./nc)+(0.0,1.0)*sin(2*pi*ky*1./nc))-&
             !  &(cos(2*pi*kx*1./nc)+(0.0,1.0)*sin(2*pi*kx*1./nc))-&
             !  &(cos(2*pi*ky*1./nc)+(0.0,1.0)*sin(2*pi*ky*1./nc)))
          endif
       enddo
    enddo

    call sfft2_c2c(fmap,nc,-1)
    shear%p=fmap


    !-----------
    ! Get gamma2

    fmap=kappa(:,:)+(0.0,0.0)

    write(*,*)' kappa' , sum(kappa(:,:)/nc/nc), minval(kappa(:,:)),maxval(kappa(:,:)) 
    write(*,*)' fmap (real) before fft' , sum(real(fmap(:,:))/nc/nc), minval(real(fmap(:,:))),maxval(real(fmap(:,:))) 
    write(*,*)' fmap (imag) before fft' , sum(aimag(fmap(:,:))/nc/nc), minval(aimag(fmap(:,:))),maxval(aimag(fmap(:,:))) 

    call sfft2_c2c(fmap,nc,1)

    write(*,*)'fmap (real) after fft' , sum(real(fmap(:,:))/nc/nc), minval(real(fmap(:,:))),maxval(real(fmap(:,:))) 
    write(*,*)'fmap (imag) after fft' , sum(aimag(fmap(:,:))/nc/nc), minval(aimag(fmap(:,:))),maxval(aimag(fmap(:,:))) 

    do kx=0,nc-1
       do ky=0,nc-1
          if((kx.eq.0).and.(ky.eq.0)) then
             fmap(1,1) = (0.0,0.0)
          else
             fmap(kx+1,ky+1)=-1.0/2*fmap(kx+1,ky+1)*(1**2)&
               &/((sin(kx*pi*1./nc))**2+(sin(ky*pi*1./nc))**2)
             !fgama1(kx+1,ky+1)=fphi(kx+1,ky+1)*(-2.0)/(1**2)&
             !  &*((sin(kx*pi*1./nc))**2-(sin(ky*pi*1./nc))**2)
             fmap(kx+1,ky+1)=fmap(kx+1,ky+1)/(1**2)&
               &*(1+(cos(2*pi*kx*1./nc)+(0.0,1.0)*sin(2*pi*kx*1./nc))*&
               &(cos(2*pi*ky*1./nc)+(0.0,1.0)*sin(2*pi*ky*1./nc))-&
               &(cos(2*pi*kx*1./nc)+(0.0,1.0)*sin(2*pi*kx*1./nc))-&
               &(cos(2*pi*ky*1./nc)+(0.0,1.0)*sin(2*pi*ky*1./nc)))
          endif
       enddo
    enddo

    call sfft2_c2c(fmap,nc,-1)
    shear%c=fmap

    write(*,*)' Unzoomed shear p mean,min,max:' , sum(shear(:,:)%p/nc/nc), minval(shear(:,:)%p),maxval(shear(:,:)%p) 
    write(*,*)' Unzoomed shear c mean,min,max:' , sum(shear(:,:)%c/nc/nc), minval(shear(:,:)%c),maxval(shear(:,:)%c) 


    return
  end subroutine kappa_to_shear

 subroutine kappa_to_defl(kappa,defl,shear,phi,fmap,nc)
    implicit none

    integer nc,nt
    real, dimension(nc,nc) :: kappa
    real,dimension(nc,nc) :: phi
    type(vec2D), dimension(nc,nc) :: defl
    type(Spin2), dimension(nc,nc) :: shear 

    complex, dimension(nc,nc) :: fmap,fphi,fgama1,fgama2
    integer kx,ky,i1,j1
    real pi

    pi=acos(-1.)
    fmap=kappa(:,:)+(0.0,0.0)

    write(*,*)' kappa' , sum(kappa(:,:)/nc/nc), minval(kappa(:,:)),maxval(kappa(:,:)) 
    write(*,*)' fmap (real) before fft' , sum(real(fmap(:,:))/nc/nc), minval(real(fmap(:,:))),maxval(real(fmap(:,:))) 
    write(*,*)' fmap (imag) before fft' , sum(aimag(fmap(:,:))/nc/nc), minval(aimag(fmap(:,:))),maxval(aimag(fmap(:,:))) 

    call sfft2_c2c(fmap,nc,1)

    write(*,*)'fmap (real) after fft' , sum(real(fmap(:,:))/nc/nc), minval(real(fmap(:,:))),maxval(real(fmap(:,:))) 
    write(*,*)'fmap (imag) after fft' , sum(aimag(fmap(:,:))/nc/nc), minval(aimag(fmap(:,:))),maxval(aimag(fmap(:,:))) 
 
    nt = 24
    write(*,*) 'Looping with openmp, n_threads = ', nt
    !$ call omp_set_num_threads(nt)
    !$omp parallel do default(private) shared(fphi,fgama1,fgama2) 
    do kx=0,nc-1
       do ky=0,nc-1
          if((kx.eq.0).and.(ky.eq.0)) then
             fphi(1,1)=(0.0,0.0)
             fgama1(1,1)=(0.0,0.0)
             fgama2(1,1)=(0.0,0.0)
          else
             fphi(kx+1,ky+1)=-1.0/2*fmap(kx+1,ky+1)*(1**2)&
               &/((sin(kx*pi*1./nc))**2+(sin(ky*pi*1./nc))**2)
             fgama1(kx+1,ky+1)=fphi(kx+1,ky+1)*(-2.0)/(1**2)&
               &*((sin(kx*pi*1./nc))**2-(sin(ky*pi*1./nc))**2)
             fgama2(kx+1,ky+1)=fphi(kx+1,ky+1)/(1**2)&
               &*(1+(cos(2*pi*kx*1./nc)+(0.0,1.0)*sin(2*pi*kx*1./nc))*&
               &(cos(2*pi*ky*1./nc)+(0.0,1.0)*sin(2*pi*ky*1./nc))-&
               &(cos(2*pi*kx*1./nc)+(0.0,1.0)*sin(2*pi*kx*1./nc))-&
               &(cos(2*pi*ky*1./nc)+(0.0,1.0)*sin(2*pi*ky*1./nc)))
          endif
       enddo
    enddo
    !$omp end parallel do 
    write(*,*) 'Got Fourier maps'

    call sfft2_c2c(fphi,nc,-1)
    phi=fphi
    write(*,*) 'Got phi'
    call sfft2_c2c(fgama1,nc,-1)
    shear%p=fgama1
    write(*,*) 'Got gamma1'
    call sfft2_c2c(fgama2,nc,-1)
    shear%c=fgama2
    write(*,*) 'Got gamma2'

    write(*,*)' Unzoomed shear c mean,min,max:' , sum(shear(:,:)%c/nc/nc), minval(shear(:,:)%c),maxval(shear(:,:)%c) 
    write(*,*)' Unzoomed shear p mean,min,max:' , sum(shear(:,:)%p/nc/nc), minval(shear(:,:)%p),maxval(shear(:,:)%p) 

    !!!$omp parallel do default(private) shared(defl) 
    do j1=1,nc
       do i1=1,nc
          defl(i1,j1)%x=(phi(modulo(i1,nc)+1,j1)-phi(modulo(i1-2,nc)+1,j1))/2
          defl(i1,j1)%y=(phi(i1,modulo(j1,nc)+1)-phi(i1,modulo(j1-2,nc)+1))/2
       enddo
    enddo
    !!!$omp end parallel do
    write(*,*) 'Got defl'
    write(*,*)' Unzoomed defl x mean,min,max:' , sum(defl(:,:)%x/nc/nc), minval(defl(:,:)%x),maxval(defl(:,:)%x) 

    return
  end subroutine kappa_to_defl

   !---
  subroutine kappa_to_shear_nodefl(kappa,shear,phi,fmap,nc)
    implicit none

    integer nc,nt
    real, dimension(nc,nc) :: kappa
    real,dimension(nc,nc) :: phi
    type(Spin2), dimension(nc,nc) :: shear 

    complex, dimension(nc,nc) :: fmap,fphi,fgama1,fgama2
    integer kx,ky,i1,j1
    real pi

    pi=acos(-1.)
    fmap=kappa(:,:)+(0.0,0.0)

    write(*,*)' kappa' , sum(kappa(:,:)/nc/nc), minval(kappa(:,:)),maxval(kappa(:,:)) 
    write(*,*)' fmap (real) before fft' , sum(real(fmap(:,:))/nc/nc), minval(real(fmap(:,:))),maxval(real(fmap(:,:))) 
    write(*,*)' fmap (imag) before fft' , sum(aimag(fmap(:,:))/nc/nc), minval(aimag(fmap(:,:))),maxval(aimag(fmap(:,:))) 

    call sfft2_c2c(fmap,nc,1)

    write(*,*)'fmap (real) after fft' , sum(real(fmap(:,:))/nc/nc), minval(real(fmap(:,:))),maxval(real(fmap(:,:))) 
    write(*,*)'fmap (imag) after fft' , sum(aimag(fmap(:,:))/nc/nc), minval(aimag(fmap(:,:))),maxval(aimag(fmap(:,:))) 
 
    nt = 24
    write(*,*) 'Looping with openmp, n_threads = ', nt
    !$ call omp_set_num_threads(nt)
    !$omp parallel do default(private) shared(fphi,fgama1,fgama2) 
    do kx=0,nc-1
       do ky=0,nc-1
          if((kx.eq.0).and.(ky.eq.0)) then
             fphi(1,1)=(0.0,0.0)
             fgama1(1,1)=(0.0,0.0)
             fgama2(1,1)=(0.0,0.0)
          else
             fphi(kx+1,ky+1)=-1.0/2*fmap(kx+1,ky+1)*(1**2)&
               &/((sin(kx*pi*1./nc))**2+(sin(ky*pi*1./nc))**2)
             fgama1(kx+1,ky+1)=fphi(kx+1,ky+1)*(-2.0)/(1**2)&
               &*((sin(kx*pi*1./nc))**2-(sin(ky*pi*1./nc))**2)
             fgama2(kx+1,ky+1)=fphi(kx+1,ky+1)/(1**2)&
               &*(1+(cos(2*pi*kx*1./nc)+(0.0,1.0)*sin(2*pi*kx*1./nc))*&
               &(cos(2*pi*ky*1./nc)+(0.0,1.0)*sin(2*pi*ky*1./nc))-&
               &(cos(2*pi*kx*1./nc)+(0.0,1.0)*sin(2*pi*kx*1./nc))-&
               &(cos(2*pi*ky*1./nc)+(0.0,1.0)*sin(2*pi*ky*1./nc)))
          endif
       enddo
    enddo
    !$omp end parallel do 
    write(*,*) 'Got Fourier maps'

    !call sfft2_c2c(fphi,nc,-1)
    !phi=fphi
    !write(*,*) 'Got phi'
    call sfft2_c2c(fgama1,nc,-1)
    shear%p=fgama1
    write(*,*) 'Got gamma1'
    call sfft2_c2c(fgama2,nc,-1)
    shear%c=fgama2
    write(*,*) 'Got gamma2'

    write(*,*)' Unzoomed shear c mean,min,max:' , sum(shear(:,:)%c/nc/nc), minval(shear(:,:)%c),maxval(shear(:,:)%c) 
    write(*,*)' Unzoomed shear p mean,min,max:' , sum(shear(:,:)%p/nc/nc), minval(shear(:,:)%p),maxval(shear(:,:)%p) 

    return
  end subroutine kappa_to_shear_nodefl

  !---

  subroutine kappa2gamma(kappa,defl,shear,phi,nc)
    implicit none

    integer nc
    real,dimension(nc,nc) :: kappa,phi
    !type(vec2D_r8), dimension(nc,nc) :: defl
    type(vec2D), dimension(nc,nc) :: defl
    type(Spin2), dimension(nc,nc) :: shear

    complex, dimension(2*nc,2*nc) :: fmap,fphi,fgamma1,fgamma2
    integer nx,ny,i1,j1
    real pi,kx,ky

    pi=acos(-1.)

    !Get kappa map in Fourier space
    !Work in a space twice as big!
    !fmap = kappa(:,:) + (0.0,0.0)
    !call sfft2_c2c(fmap,nc,1)
 
    fmap = 0.0
    fmap(nc/2+1:3*nc/2, nc/2+1:3*nc/2)=kappa(1:nc,1:nc)
    call sfft2_c2c(fmap,2*nc,1)
 
    !Get Graviational potential and shear maps in Fourier space

    do nx=1,2*nc
       kx = real(nx-1)*2.0*pi/(2*nc)
       do ny=1,nc
          ky = real(ny-1)*2.0*pi/(2*nc)
          if((kx.eq.0).and.(ky.eq.0)) then
             fphi(1,1)=(0.0,0.0)
             fgamma1(1,1)=(0.0,0.0)
             fgamma2(1,1)=(0.0,0.0)
          else
             fphi(nx,ny)=-2.0*fmap(nx,ny)/(kx**2+ky**2)
             fgamma1(nx,ny)=-0.5*fphi(nx,ny)*(kx**2-ky**2)
             fgamma2(nx,ny)=-1.0*kx*ky*fphi(nx,ny)
          endif
       enddo
    enddo


    ! go back in real space
    call sfft2_c2c(fphi,2*nc,-1)
    phi=fphi(nc/2+1:3*nc/2,nc/2+1:3*nc/2)
    call sfft2_c2c(fgamma1,2*nc,-1)
    call sfft2_c2c(fgamma2,2*nc,-1)

    shear%p=fgamma1(nc/2+1:3*nc/2,nc/2+1:3*nc/2)
    shear%c=fgamma2(nc/2+1:3*nc/2,nc/2+1:3*nc/2)

    do j1=1,nc
       do i1=1,nc
          !defl(i1,j1)%x=phi(modulo(i1,nc)+1,j1)-phi(i1,j1)
          !defl(i1,j1)%y=phi(i1,modulo(j1,nc)+1)-phi(i1,j1)
          defl(i1,j1)%x=(phi(modulo(i1,nc)+1,j1)-phi(modulo(i1-2,nc)+1,j1))/2
          defl(i1,j1)%y=(phi(i1,modulo(j1,nc)+1)-phi(i1,modulo(j1-2,nc)+1))/2
       enddo
    enddo

  
  end subroutine kappa2gamma

  subroutine shear_to_power(shear,power,nc)
    implicit none

    integer nc,i,j,i1,j1,ir
    type(Spin2), dimension(nc,nc) :: shear
    real a(nc,nc),b(nc,nc),ps(nc),w,r,weight(nc),tmp,cx,cy
    real, dimension(nc) :: power
    complex fa(nc,nc),fb(nc,nc)
  
  
    fa=shear%p+(0.,0.)
    fb=shear%c+(0.,0.)
    call sfft2_c2c(fa,nc,1)
    call sfft2_c2c(fb,nc,1)
  
    fa=fa/nc**2
    fb=fb/nc**2
    weight=0
    ps=0
    do j=1,nc
       j1=j-1
       if(j1>nc/2) j1=j1-nc
       do i=1,nc
          i1=i-1
          if(i1>nc/2) i1=i1-nc
          r=(i1**2+j1**2)**0.5
          ir=r
          cx=0
          cy=0
          if(ir>=1) then
            cx=i1/r
            cy=j1/r
            !tmp=fa(i,j)*conjg(fa(i,j))*((cx-cy)*(cx+cy))**2&
            !    +fb(i,j)*conjg(fb(i,j))*(2*cx*cy)**2
            !tmp=fa(i,j)*conjg(fa(i,j))&
            !    +fb(i,j)*conjg(fb(i,j))
            tmp=(fa(i,j)*(cx**2-cy**2)+fb(i,j)*(2*cx*cy))&
                *conjg(fa(i,j)*(cx**2-cy**2)+fb(i,j)*(2*cx*cy))
            w=r-ir
            ps(ir)=ps(ir)+tmp*(1-w)
            ps(ir+1)=ps(ir+1)+tmp*w
            weight(ir)=weight(ir)+1-w
            weight(ir+1)=weight(ir+1)+w
          endif
       enddo
    enddo
    power=0
    do i=1,nc
       if(weight(i)>0) power(i)=ps(i)/weight(i)*2*3.1415927*i**2
       !write(*,*) i,power(i)
       !pause 
    end do


    return
  end subroutine shear_to_power

end module Lensing
