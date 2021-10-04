! Calculate Redshift to a Comoving Distance in Mpc/h
! using fitting formula from Ue-Li Pen
! for matter + Dark Energy era
!
!Specify Chi in Mpc/h;

function InverseChi(chi,omegam)

  use Solve_Real_Poly

  implicit none
  real InverseChi,chi,omegam,s,c,H0,D_h, alpha1,alpha2,alpha3,alpha4,alpha5
  integer degree
  real, dimension(5)::p,zeror,zeroi
  logical fail

  degree = 5

  write(*,*) 'Calling InverseChi'

  s = ((1-omegam)/omegam)**(1./3.)

  c=3E5
  H0=100.
  D_h=c/H0

  alpha1=-0.1540*s;
  alpha2=0.4304*s*s;
  alpha3=0.19097*s*s*s;
  alpha4=0.066941*s*s*s*s;
  alpha5=( (eta(0,omegam)-1*chi/D_h)/(2*sqrt(s**3+1)))**(-8); 

  ! expand in power of Z, not Z+1
  !p = [1 (4+alpha1) (4+3*alpha1+alpha2) (4+3*alpha1+2*alpha2+alpha3) (1+alpha1+alpha2+alpha3+alpha4-alpha5)];

  ! expand in powers of Z+1
  p(1)= 1. 
  p(2)= alpha1
  p(3)= alpha2
  p(4)= alpha3
  p(5)= alpha4-alpha5

  !r = max(real(roots(p)));
  call rpoly(p,degree,zeror, zeroi,fail)  
  write(*,*)'root = ',zeror,zeroi,fail 

  InverseChi = maxval(zeror)-1;

end function InverseChi
