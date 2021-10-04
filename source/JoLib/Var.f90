
!! This subroutine is to measure the mean variance of a 3D box,
!! and the 3D field is smoothed by a function w
Subroutine Varm(den,nc,box,w,ps,mvar)
  implicit none

  integer nc
  real box, w, mvar
  real, dimension(nc) :: ps
  real, dimension(nc+2,nc,nc) :: den 
  
  external spline,splint 

  !integer, parameter :: nk=ceiling(nc/2*3**0.5)
  !real, dimension(0:nk) :: pk,kg,dpk,dp2
  real, dimension(0:ceiling(nc/2*3**0.5)) :: pk,kg,dpk,dp2

  integer nk,hc,i,j,k,kx,ky,kz,k1,k2,nk1,nk2,nr2c
  real pi,npk,kri,kr
  real(kind=8) mvar_D

  pi=acos(-1.)
  hc=nc/2
  nk=ceiling(nc/2*3**0.5)

  write(*,*) 'mean', 2*sum(real(den**2,kind=8)/nc/nc/nc)
  call ps3_r2c(den,ps,nc)
 
  pk=0
  pk(1:nk)=ps(1:nk) 
  do k=1,nk
     kg(k)=2*pi/box*k
     pk(k)=pk(k)*2*pi**2/kg(k)**3
  enddo
 
  dpk=0
  dpk(0)=(pk(1)-pk(0))/(kg(1)-kg(0))
  dpk(nk)=(pk(nk)-pk(nk-1))/(kg(nk)-kg(nk-1))
  call spline(kg,pk,nk+1,dpk(0),dpk(nk),dp2)

  mvar_D=0
  do k=1,nc
     kz=k-1
     if(kz>hc)kz=k-1-nc
     do j=1,nc
        ky=j-1
        if(ky>hc)ky=j-1-nc
        do i=1,nc+1,2
           kx=(i-1)/2
           kr=sqrt(1.*kx**2+ky**2+kz**2)*2*pi/box
           kri=sqrt(1.*kx**2+ky**2+kz**2)
           nk1=floor(kri)
           nk2=nk1+1
           if(nk1.ge.1.and.nk2.le.nk) then
              call splint(kg,pk,dp2,nk+1,kr,npk)
           else
              npk=pk(0)
           endif

           nr2c=1
           if(kx.ne.0.and.kx.ne.hc)nr2c=2
           mvar_D=mvar_D+real(nr2c*npk*w(kr)**2/box**3,kind=8)
        enddo
     enddo
  enddo
  mvar=mvar_D
  write(*,*) 'integral ps', 2*mvar

  return
end subroutine Varm
