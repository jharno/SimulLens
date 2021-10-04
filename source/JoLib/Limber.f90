
subroutine Limber(kr,delta2,nk,w,zn,nz,l,lcl,omegam,h0)
!! x_2d=/int d\chi w(\chi) x_3d
!! delta2 is the 3d dimensionless powerspectrum of x_3d
  implicit none
  integer nk,nz
  integer j
  real kr(nk,nz),delta2(nk,nz),zn(nz),w(nz),y(nz),l,lcl
  real chi,omegam,omegav,h0
  real k0,delta20,x,H
  real z,pi

  external chi, linear_interpolation

  pi=acos(-1.) 
  omegav=1-omegam

  lcl=0
  do j=1,nz
     z=zn(j)
     x=chi(z,omegam,h0)*h0
     k0=l/x
     if(k0>kr(nk,j))then
        delta20=0.0
     elseif(k0<kr(1,j))then
        delta20=0.0
     else
       call linear_interpolation(kr(:,j),delta2(:,j),Nk,k0,delta20)
     endif
     !! NOTE we replace H0/c by 1/(3.E3) in the unit of h Mpc^{-1}
     H=1/(3.E3)*((1+z)**2*(omegam*z+1)-omegav*z*(z+2))**(1.0/2)
     y(j)=w(j)**2*x*delta20/H
     if(j.gt.1) then
        lcl=lcl+(y(j)+y(j-1))*(zn(j-1)-zn(j))/2
     endif
  enddo
  lcl=pi*(l+1)/l**2*lcl

  return
end subroutine Limber

