
!! linear_interpolation is from Pengjie zhang

subroutine linear_interpolation(xa,ya,n,x,y)
  implicit none
  integer n,i,k,khi,klo
  real x,xa(n),ya(n),y,a,b,h

  if(x>xa(n)) then
    !write(*,*) ' beyond the upper range of dataset'
    y=xa(n)
  else if(x<xa(1)) then
    !write(*,*) ' beyond the lower range of dataset'
    y=xa(1)
  end if

  klo=1
  khi=n
  1 if(khi-klo .gt. 1) then
      k=(khi+klo)/2
      if(xa(k) .gt. x) then
  khi=k
      else
  klo=k
      end if
  go to 1
  end if
  h=xa(khi)-xa(klo)
  if (h .eq. 0.) pause 'bad xa input in splint'
  ! the xa's must be distinct
  a=(xa(khi)-x)/h
  b=(x-xa(klo))/h
  y=a*ya(klo)+b*ya(khi)
  
  return
end subroutine linear_interpolation


subroutine cic2d(den,n,xx,yy,rho)
  implicit none

  integer n
  real den(n,n)
  real xx,yy,rho

  integer i1,i2,j1,j2
  real x,y,dx1,dx2,dy1,dy2,rn

  rn=n
  x=modulo(real(xx)-0.5+rn,rn)
  y=modulo(real(yy)-0.5+rn,rn)
  i1=int(x)+1
  i2=mod(i1,n)+1
  dx1=i1-x
  dx2=1-dx1
  j1=int(y)+1
  j2=mod(j1,n)+1
  dy1=j1-y
  dy2=1-dy1

  rho=0
  rho=rho+den(i1,j1)*dx1*dy1
  rho=rho+den(i2,j1)*dx2*dy1
  rho=rho+den(i1,j2)*dx1*dy2
  rho=rho+den(i2,j2)*dx2*dy2

  return
end subroutine cic2d


subroutine tsc3d(den,nc,xx,yy,zz,rho)
  implicit none

  integer nc,i0,in1,i1,j0,jn1,j1,k0,kn1,k1
  real xx,yy,zz,rho,den(nc,nc,nc)
  real ncr,x,y,z,dxn1,dx0,dx1,dyn1,dy0,dy1,dzn1,dz0,dz1

  ncr=nc
  !i0=NINT(xv(1,ip))
  !j0=NINT(xv(2,ip))
  !k0=NINT(xv(3,ip))
  !x=xv(1,ip)-i0
  !y=xv(2,ip)-j0
  !z=xv(3,ip)-k0
  i0=modulo(NINT(xx-0.5),nc)+1
  j0=modulo(NINT(yy-0.5),nc)+1
  k0=modulo(NINT(zz-0.5),nc)+1
  x=xx-0.5-NINT(xx-0.5)
  y=yy-0.5-NINT(yy-0.5)
  z=zz-0.5-NINT(zz-0.5)

  i0=mod(i0-1+nc,nc)+1
  j0=mod(j0-1+nc,nc)+1
  k0=mod(k0-1+nc,nc)+1
  in1=mod(i0-1-1+nc,nc)+1
  i1=mod(i0+1-1+nc,nc)+1
  jn1=mod(j0-1-1+nc,nc)+1
  j1=mod(j0+1-1+nc,nc)+1
  kn1=mod(k0-1-1+nc,nc)+1
  k1=mod(k0+1-1+nc,nc)+1
  dxn1=1./2*(3./2-(x+1))**2
  dx0=3./4-x**2
  dx1=1./2*(3./2-(1-x))**2
  dyn1=1./2*(3./2-(y+1))**2
  dy0=3./4-y**2
  dy1=1./2*(3./2-(1-y))**2
  dzn1=1./2*(3./2-(z+1))**2
  dz0=3./4-z**2
  dz1=1./2*(3./2-(1-z))**2

  rho=0
  rho=rho+den(i0,j0,k0)*dx0*dy0*dz0
  rho=rho+den(in1,j0,k0)*dxn1*dy0*dz0
  rho=rho+den(i1,j0,k0)*dx1*dy0*dz0
  rho=rho+den(i0,j0,kn1)*dx0*dy0*dzn1
  rho=rho+den(in1,j0,kn1)*dxn1*dy0*dzn1
  rho=rho+den(i1,j0,kn1)*dx1*dy0*dzn1
  rho=rho+den(i0,j1,k0)*dx0*dy1*dz0
  rho=rho+den(in1,j1,k0)*dxn1*dy1*dz0
  rho=rho+den(i1,j1,k0)*dx1*dy1*dz0
  rho=rho+den(i0,j1,kn1)*dx0*dy1*dzn1
  rho=rho+den(in1,j1,kn1)*dxn1*dy1*dzn1
  rho=rho+den(i1,j1,kn1)*dx1*dy1*dzn1
  rho=rho+den(i0,j1,k1)*dx0*dy1*dz1
  rho=rho+den(in1,j1,k1)*dxn1*dy1*dz1
  rho=rho+den(i1,j1,k1)*dx1*dy1*dz1
  rho=rho+den(i0,jn1,k0)*dx0*dyn1*dz0
  rho=rho+den(in1,jn1,k0)*dxn1*dyn1*dz0
  rho=rho+den(i1,jn1,k0)*dx1*dyn1*dz0
  rho=rho+den(i0,jn1,k1)*dx0*dyn1*dz1
  rho=rho+den(in1,jn1,k1)*dxn1*dyn1*dz1
  rho=rho+den(i1,jn1,k1)*dx1*dyn1*dz1
  rho=rho+den(i0,jn1,kn1)*dx0*dyn1*dzn1
  rho=rho+den(in1,jn1,kn1)*dxn1*dyn1*dzn1
  rho=rho+den(i1,jn1,kn1)*dx1*dyn1*dzn1
  rho=rho+den(i0,j0,k1)*dx0*dy0*dz1
  rho=rho+den(in1,j0,k1)*dxn1*dy0*dz1
  rho=rho+den(i1,j0,k1)*dx1*dy0*dz1

  return
end subroutine tsc3d

subroutine tsc2d(den,nc,xx,yy,rho)
  implicit none

  integer nc,i0,in1,i1,j0,jn1,j1
  real xx,yy,rho,den(nc,nc)
  real ncr,x,y,dxn1,dx0,dx1,dyn1,dy0,dy1

  ncr=nc
  !i0=NINT(xv(1,ip))
  !j0=NINT(xv(2,ip))
  !x=xv(1,ip)-i0
  !y=xv(2,ip)-j0
  i0=modulo(NINT(xx-0.5),nc)+1
  j0=modulo(NINT(yy-0.5),nc)+1
  x=xx-0.5-NINT(xx-0.5)
  y=yy-0.5-NINT(yy-0.5)

  i0=mod(i0-1+nc,nc)+1
  j0=mod(j0-1+nc,nc)+1
  in1=mod(i0-1-1+nc,nc)+1
  i1=mod(i0+1-1+nc,nc)+1
  jn1=mod(j0-1-1+nc,nc)+1
  j1=mod(j0+1-1+nc,nc)+1
  dxn1=1./2*(3./2-(x+1))**2
  dx0=3./4-x**2
  dx1=1./2*(3./2-(1-x))**2
  dyn1=1./2*(3./2-(y+1))**2
  dy0=3./4-y**2
  dy1=1./2*(3./2-(1-y))**2

  rho=0
  rho=rho+den(i0,j0)*dx0*dy0
  rho=rho+den(in1,j0)*dxn1*dy0
  rho=rho+den(i1,j0)*dx1*dy0
  rho=rho+den(i0,j1)*dx0*dy1
  rho=rho+den(in1,j1)*dxn1*dy1
  rho=rho+den(i1,j1)*dx1*dy1
  rho=rho+den(i0,jn1)*dx0*dyn1
  rho=rho+den(in1,jn1)*dxn1*dyn1
  rho=rho+den(i1,jn1)*dx1*dyn1

  return
end subroutine tsc2d


subroutine tsc1d(den,nc,xx,rho)
  implicit none

  integer nc,i0,in1,i1
  real xx,rho,den(nc)
  real ncr,x,dxn1,dx0,dx1

  ncr=nc
  !i0=NINT(xv(1,ip))
  !x=xv(1,ip)-i0
  i0=modulo(NINT(xx-0.5),nc)+1
  x=xx-0.5-NINT(xx-0.5)

  i0=mod(i0-1+nc,nc)+1
  in1=mod(i0-1-1+nc,nc)+1
  i1=mod(i0+1-1+nc,nc)+1
  dxn1=1./2*(3./2-(x+1))**2
  dx0=3./4-x**2
  dx1=1./2*(3./2-(1-x))**2

  rho=0
  rho=rho+den(i0)*dx0
  rho=rho+den(in1)*dxn1
  rho=rho+den(i1)*dx1

  return
end subroutine tsc1d

subroutine tsc_cic_2d(den,nc,xx,yy,rho)
  implicit none

  integer nc,i0,in1,i1,j1,j2
  real xx,yy,rho,den(nc,nc)
  real ncr,x,y,dxn1,dx0,dx1,dy1,dy2

  ncr=nc
  i0=modulo(NINT(xx-0.5),nc)+1
  x=xx-0.5-NINT(xx-0.5)
  y=modulo(real(yy)-0.5+ncr,ncr)

  i0=mod(i0-1+nc,nc)+1
  in1=mod(i0-1-1+nc,nc)+1
  i1=mod(i0+1-1+nc,nc)+1
  dxn1=1./2*(3./2-(x+1))**2
  dx0=3./4-x**2
  dx1=1./2*(3./2-(1-x))**2
  j1=int(y)+1
  j2=mod(j1,nc)+1
  dy1=j1-y
  dy2=1-dy1

  rho=0
  rho=rho+den(i0,j1)*dx0*dy1
  rho=rho+den(in1,j1)*dxn1*dy1
  rho=rho+den(i1,j1)*dx1*dy1
  rho=rho+den(i0,j2)*dx0*dy2
  rho=rho+den(in1,j2)*dxn1*dy2
  rho=rho+den(i1,j2)*dx1*dy2

  return
end subroutine tsc_cic_2d

subroutine cic_tsc_2d(den,nc,xx,yy,rho)
  implicit none

  integer nc,j0,jn1,j1,i1,i2
  real xx,yy,rho,den(nc,nc)
  real ncr,x,y,dyn1,dy0,dy1,dx1,dx2
  
  ncr=nc
  x=modulo(real(xx)-0.5+ncr,ncr)
  j0=modulo(NINT(yy-0.5),nc)+1
  y=yy-0.5-NINT(yy-0.5)

  i1=int(x)+1
  i2=mod(i1,nc)+1
  dx1=i1-x
  dx2=1-dx1
  j0=mod(j0-1+nc,nc)+1
  jn1=mod(j0-1-1+nc,nc)+1
  j1=mod(j0+1-1+nc,nc)+1
  dyn1=1./2*(3./2-(y+1))**2
  dy0=3./4-y**2
  dy1=1./2*(3./2-(1-y))**2

  rho=0
  rho=rho+den(i1,j0)*dx1*dy0
  rho=rho+den(i1,jn1)*dx1*dyn1
  rho=rho+den(i1,j1)*dx1*dy1
  rho=rho+den(i2,j0)*dx2*dy0
  rho=rho+den(i2,jn1)*dx2*dyn1
  rho=rho+den(i2,j1)*dx2*dy1


  return
end subroutine cic_tsc_2d
