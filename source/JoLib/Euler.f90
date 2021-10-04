!! Active Rotation by the 3 Euler Angles.
! General       -> 1 Rotation about z_axis + 1 Rotation about y' axis + 1 Rotation about z'' 
! psi = 0       -> 1 Rotation about z-axis + 1 Rotation about y' axis
! psi,theta = 0 -> phi rotation about z-axis (clockwise)
! psi, phi = 0  -> Rotation about y axis (clockwise) by theta
! phi = pi/2, psi = -pi/2 -> Rotation about x-axis (counter clockwise) by theta

!function euler(vector,phi,theta,psi)
subroutine Rotate(vector,phi,theta,psi,Euler)

  implicit none
  real phi, theta, psi
  real, dimension(3) :: vector,euler
  real, dimension(3,3) :: Matrix

  Matrix(1,1) = cos(phi)*cos(theta)*cos(psi) - sin(phi)*sin(psi)
  Matrix(1,2) = sin(phi)*cos(theta)*cos(psi) + cos(phi)*sin(psi)
  Matrix(1,3) = -sin(theta)*cos(psi)

  Matrix(2,1) = -cos(phi)*cos(theta)*sin(psi) - sin(phi)*cos(psi)
  Matrix(2,2) = -sin(phi)*cos(theta)*sin(psi) + cos(phi)*cos(psi)
  Matrix(2,3) = sin(theta)*sin(psi)

  Matrix(3,1) = cos(phi)*sin(theta)
  Matrix(3,2) = sin(phi)*sin(theta)
  Matrix(3,3) = cos(theta)

  Euler(1) = Matrix(1,1)*vector(1) + Matrix(1,2)*vector(2) + Matrix(1,3)*vector(3)
  Euler(2) = Matrix(2,1)*vector(1) + Matrix(2,2)*vector(2) + Matrix(2,3)*vector(3)
  Euler(3) = Matrix(3,1)*vector(1) + Matrix(3,2)*vector(2) + Matrix(3,3)*vector(3)

!end function euler
end subroutine Rotate
