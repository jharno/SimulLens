!! integral function of z to get chi, the unit is Mpc

function chi(z,omegam,h)

  implicit none
  real(4) z,omegam,h,chi
  real(8) s,a,H_0,yita1,yitaa,D_L,cc
  cc=3*10**5
  H_0=h*100
  s=((1-omegam)/omegam)**(1.0/3)
  yita1=2*(s**3+1)**(0.5)*(1.-0.1540*s+0.4304*s**2+0.19097*s**3+&
        &0.066941*s**4)**(-1.0/8)
  a=1./(z+1.)
  yitaa=2*(s**3+1)**(0.5)*(1./a**4-0.1540*s/a**3+0.4304*s**2/a**2+&
       &0.19097*s**3/a+0.066941*s**4)**(-1.0/8.)
  D_L=cc/H_0*(1+z)*(yita1-yitaa)
  chi=a*D_L
  
end function chi
