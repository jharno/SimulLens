
!! The code is to generate a Gaussian stochastic field with a given 3d power spectrum
!! written by T.T.Lu on Oct 13, 2006 
!! compile with 'efc GaussRandomField.f90 -L/opt/fftw-3.1.2_intel8/lib/ -lfftw3f -I/opt/fftw-3.1.2_intel8/include/ -o GaussRandomField.x' 
!! replace the fftw library by the corresponding one on the machine used

Subroutine GaussRandomField_3d_r2c(den, box, nc, kg, delta2, nk)
  implicit none

  !! nc is the number of grids on each dimension
  !! nt is the number of threads for openmp 
  !! nk is the number of points at the input power spectrum table
  !! box is the physical size of the box
  integer nc,nk
  integer, parameter :: nt=1

  real, dimension(nc+2,nc,nc) :: den
  real, dimension(nk) :: pk,kg,delta2
!  real, dimension(0:nk) :: factor
  real, dimension(0:nc) :: factor
  real, dimension(nc,nc,nc) :: x1,x2
  real box 

  integer hc,it,i,j,k,k1,k2,kx,ky,kz,kpt
  real pi,kr,w1,w2

  external sfft3_r2c

  write(*,*)'Called GaussRandomField_3d_r2c'


  pi=acos(-1.)
  hc=nc/2

  do i=1,nk
     pk(i)=delta2(i)/kg(i)**3*2*pi**2
  enddo


  !! generate whitenoise source
  call random_number(x1(:,:,:))
  call random_number(x2(:,:,:))

  den=0
  do k=1,nc
     do j=1,nc
        do i=1,nc
           den(i,j,k)=sqrt(-2*log(x1(i,j,k)))*cos(2*pi*x2(i,j,k))
        enddo
     enddo
  enddo
  !write(*,*) 'step3' 

  call sfft3_r2c(den,nc,1)
  !write(*,*) 'called  fft forward'

  !write(*,*) '1 :Max = ' , maxval(den), 'min = ', minval(den)


  factor=0
  do k=1,nk
     factor(k)=pk(k)/(box**3/nc**3)
  enddo
 
  kpt=nc/nt
  !!$omp parallel do default(shared) &
  !!$omp& private(it,i,j,k,kr,kri,kx,ky,kz,k1,k2,w1,w2,tmp)
  do it=1,nt
     do k=1+(it-1)*kpt,min(nc,it*kpt)
        if (k .le. hc+1) then
            kz=k-1
        else
            kz=k-1-nc
        endif
        do j=1,nc
            if (j .le. hc+1) then
                ky=j-1
            else
                ky=j-1-nc
            endif
            do i=1,nc+2,2
               kx=(i-1)/2
               kr=sqrt(kx**2.+ky**2+kz**2)
               if(kr.ne.0) then
                 k1=ceiling(kr)
                 k2=k1+1
                 w1=k1-kr
                 w2=1-w1
                 den(i,j,k)=den(i,j,k)*sqrt(factor(k1-1)*w1+factor(k2-1)*w2) !! interpolation causes problem
                 den(i+1,j,k)=den(i+1,j,k)*sqrt(factor(k1-1)*w1+factor(k2-1)*w2)
               else
                 den(i,j,k)=den(i,j,k)*sqrt(factor(0))
                 den(i+1,j,k)=den(i+1,j,k)*sqrt(factor(0))
               endif
           enddo
        enddo
     enddo
  enddo
  !!$omp end parallel do

!  write(*,*) '2 :Max = ' , maxval(den), 'min = ', minval(den)
  call sfft3_r2c(den,nc,-1)
!  write(*,*) '3 :Max = ' , maxval(den), 'min = ', minval(den)

!  write(*,*) 'called  fft backward'

  return
end subroutine GaussRandomField_3d_r2c

Subroutine GaussRandomField_3d_r2c_KSpace(den, box, nc, kg, delta2, nk)
  implicit none

  !! nc is the number of grids on each dimension
  !! nt is the number of threads for openmp 
  !! nk is the number of points at the input power spectrum table
  !! box is the physical size of the box
  integer nc,nk
  integer, parameter :: nt=1

  real, dimension(nc+2,nc,nc) :: den
  real, dimension(nk) :: pk,kg,delta2
!  real, dimension(0:nk) :: factor
  real, dimension(0:nc) :: factor
  real, dimension(nc,nc,nc) :: x1,x2
  real box 

  integer hc,it,i,j,k,k1,k2,kx,ky,kz,kpt
  real pi,kr,w1,w2

  external sfft3_r2c

  write(*,*)'Called GaussRandomField_3d_r2c_KSpace'


  pi=acos(-1.)
  hc=nc/2

  do i=1,nk
     pk(i)=delta2(i)/kg(i)**3*2*pi**2
  enddo


  !! generate whitenoise source
  call random_number(x1(:,:,:))
  call random_number(x2(:,:,:))

  den=0
  do k=1,nc
     do j=1,nc
        do i=1,nc
           den(i,j,k)=sqrt(-2*log(x1(i,j,k)))*cos(2*pi*x2(i,j,k))
        enddo
     enddo
  enddo
  !write(*,*) 'step3' 

  call sfft3_r2c(den,nc,1)
  !write(*,*) 'called  fft forward'

  !write(*,*) '1 :Max = ' , maxval(den), 'min = ', minval(den)


  factor=0
  do k=1,nk
     factor(k)=pk(k)/(box**3/nc**3)
  enddo
 
  kpt=nc/nt
  !!$omp parallel do default(shared) &
  !!$omp& private(it,i,j,k,kr,kri,kx,ky,kz,k1,k2,w1,w2,tmp)
  do it=1,nt
     do k=1+(it-1)*kpt,min(nc,it*kpt)
        if (k .le. hc+1) then
            kz=k-1
        else
            kz=k-1-nc
        endif
        do j=1,nc
            if (j .le. hc+1) then
                ky=j-1
            else
                ky=j-1-nc
            endif
            do i=1,nc+2,2
               kx=(i-1)/2
               kr=sqrt(kx**2.+ky**2+kz**2)
               if(kr.ne.0) then
                 k1=ceiling(kr)
                 k2=k1+1
                 w1=k1-kr
                 w2=1-w1
                 den(i,j,k)=den(i,j,k)*sqrt(factor(k1-1)*w1+factor(k2-1)*w2) !! interpolation causes problem
                 den(i+1,j,k)=den(i+1,j,k)*sqrt(factor(k1-1)*w1+factor(k2-1)*w2)
               else
                 den(i,j,k)=den(i,j,k)*sqrt(factor(0))
                 den(i+1,j,k)=den(i+1,j,k)*sqrt(factor(0))
               endif
           enddo
        enddo
     enddo
  enddo
  !!$omp end parallel do

!  write(*,*) '2 :Max = ' , maxval(den), 'min = ', minval(den)
!  call sfft3_r2c(den,nc,-1)
!  write(*,*) '3 :Max = ' , maxval(den), 'min = ', minval(den)

!  write(*,*) 'called  fft backward'

  return
end subroutine GaussRandomField_3d_r2c_KSpace

Subroutine GaussRandomField_3d_c2c(den, box, nc, kg, delta2, nk)
  implicit none
!#define logcic
!#define linearcic
#define ksquare

  !! nc is the number of grids on each dimension
  !! nt is the number of threads for openmp
  !! nk is the number of points at the input power spectrum table
  !! box is the physical size of the box
  integer nc,nk
  integer, parameter :: nt=1

!!  real, dimension(0:nk) :: factor
#ifdef ksquare
  real*8, dimension(0:3*nc**2/4) :: factor
  integer, dimension(2,3*nc**2/4) :: index_shell 
  integer, dimension(3*nc**2/4) :: arr_shell 
#else
  real*8, dimension(0:nc) :: factor
#endif
  real, dimension(nk) :: pk,kg,delta2
  real, dimension(nc,nc,nc) :: x1,x2
  real, dimension(nc,nc,nc) :: den
  complex, dimension(nc,nc,nc) :: cd

  integer hc,it,i,j,k,k1,k2,kx,ky,kz,kpt,indx
  real pi,box,kr,w1,w2
  real*8 tmp

  external sfft3_c2c

  pi=acos(-1.)
  hc=nc/2

  do i=1,nk
     pk(i)=delta2(i)/kg(i)**3*2*pi**2
  enddo

  !! generate whitenoise source
  call random_number(x1(:,:,:))
  call random_number(x2(:,:,:))
  !write(*,*) maxval(x1),minval(x1)
  !write(*,*) maxval(x2),minval(x2)

  den=0
  do k=1,nc
     do j=1,nc
        do i=1,nc
           den(i,j,k)=sqrt(-2*log(x1(i,j,k)))*cos(2*pi*x2(i,j,k))
        enddo
     enddo
  enddo
  !write(*,*) maxval(den),minval(den)

  cd=den+(0.,0.)
  call sfft3_c2c(cd,nc,1)
  !write(*,*) 'called  fft forward'

  factor=0
  do k=1,nk
     !factor(k)=pk(k)/(box**3/nc**3)
     factor(k)=pk(k)/(box/nc)**3
  enddo
  write(*,*) 'nk=',nk

#ifdef ksquare
  call index_k2_3d(index_shell,arr_shell,nc)
#endif

  kpt=nc/nt
  !!$omp parallel do default(shared) &
  !!$omp& private(it,i,j,k,kr,kri,kx,ky,kz,k1,k2,w1,w2,tmp,indx)
  do it=1,nt
     do k=1+(it-1)*kpt,min(nc,it*kpt)
        !if (k .le. hc+1) then
        if (k .lt. hc+1) then
            kz=k-1
        else
            kz=k-1-nc
        endif
        do j=1,nc
            !if (j .le. hc+1) then
            if (j .lt. hc+1) then
                ky=j-1
            else
                ky=j-1-nc
            endif
            do i=1,nc
               !if (i .le. hc+1) then
               if (i .lt. hc+1) then
                  kx=i-1
               else
                  kx=i-1-nc
               endif
               kr=sqrt(kx**2+ky**2+kz**2+0.)
               !tmp=nc**1.5/sqrt(cd(i,j,k)*conjg(cd(i,j,k)))
               tmp=sqrt(nc**3/(cd(i,j,k)*conjg(cd(i,j,k))))
               if(kr.ne.0) then
#ifdef ksquare
                 indx=index_shell(2,kx**2+ky**2+kz**2)
                 if(indx.le.nk)then 
                    cd(i,j,k)=cd(i,j,k)*sqrt(factor(indx))
                 else
                    cd(i,j,k)=0. 
                 endif
#else
                 k1=ceiling(kr)
                 k2=k1+1
                 !! cic in linear scale
#ifdef linearcic
                 w1=k1-kr
                 w2=1-w1
                 cd(i,j,k)=cd(i,j,k)*sqrt(factor(k1-1)*w1+factor(k2-1)*w2) !! interpolation causes problem
                 !! cic in log scale
#endif
#ifdef logcic
                 if(k1.gt.1.and.k2.le.nk)then
                    w1=(log(k2-1.)-log(kr+0.))/(log(k2-1.)-log(k1-1.))
                    w2=1-w1
                    cd(i,j,k)=cd(i,j,k)*sqrt(exp(log(factor(k1-1))*w1+log(factor(k2-1))*w2)) !! interpolation causes problem
                    !write(*,*) kx,ky,kz,kr,k1,k2,w1,w2,log(factor(k1-1))*w1+log(factor(k2-1))*w2
                    !pause
                 elseif(k1.eq.1)then
                    cd(i,j,k)=cd(i,j,k)*sqrt(factor(1))
                    !write(*,*) kx,ky,kz,kr,k1,k2,w1,w2,factor(1)
                    !pause
                 else
                    cd(i,j,k)=0
                 endif
#endif
#endif
                 !if(int(kr).le.nk)then 
                 !  cd(i,j,k)=cd(i,j,k)*sqrt(nc**3/box**3)!! TEST WITH WHITE NOISE
                 !else
                 !  cd(i,j,k)=0  
                 !endif
                 !if(factor(k1-1)*w1+factor(k2-1)*w2.lt.nc**3/box**3)then
                 !   write(*,*) kx,ky,kz,kr, k1,k2,factor(k1-1)*w1+factor(k2-1)*w2,nc**3/box**3
                 !   pause
                 !endif 
               else
                 cd(i,j,k)=cd(i,j,k)*sqrt(factor(0))
               endif
               !cd(i,j,k)=cd(i,j,k)*tmp
           enddo
        enddo
     enddo
  enddo
  !!$omp end parallel do
  
  !open(10,file='test.dat',form='binary')
  !write(10) real((cd/nc/nc/nc)*(conjg(cd)/nc/nc/nc))
  !close(10) 

  !write(*,*) 'sum(d(**2)= (k space 1)', sum(real(cd*conjg(cd)/nc/nc/nc,8)/nc/nc/nc)

  call sfft3_c2c(cd,nc,-1)
  !write(*,*) 'called  fft backward'
  den=real(cd)
  !write(*,*) maxval(den),minval(den)

  return
end subroutine GaussRandomField_3d_c2c

Subroutine GaussRandomField_2d_r2c(den, box, nc, kg, delta2, nk, x1,x2)
  implicit none

  !! nc is the number of grids on each dimension
  !! nt is the number of threads for openmp 
  !! nk is the number of points at the input power spectrum table
  !! box is the physical size of the box
  integer nc,nk
  integer, parameter :: nt=1

  real, dimension(nc+2,nc) :: den
  real, dimension(nk) :: pk,kg,delta2
  real, dimension(0:nk) :: factor
!  real, dimension(0:nc) :: factor
  real, dimension(nc,nc) :: x1,x2
  real box 

  integer hc,it,i,j,k,k1,k2,kx,ky,kpt
  real pi,kr,w1,w2

  external sfft2_r2c

  pi=acos(-1.)
  hc=nc/2

  !do i=1,nk
  !   pk(i)=delta2(i)/kg(i)**2*2*pi
  !enddo

  factor(0) = 0
  factor=delta2*real(nc)**2/box**2
  write(*,*) 'Already Cl:' !,factor

  write(*,*)'Got factor(k)'
 
 !! generate whitenoise source
  call random_number(x1(:,:))
  call random_number(x2(:,:))

  write(*,*)'Got random numbers'

  den=0
  do j=1,nc
     do i=1,nc
        den(i,j)=sqrt(-2*log(x1(i,j)))*cos(2*pi*x2(i,j))
     enddo
  enddo

  write(*,*)'Got Gaussian white noise', minval(den),maxval(den),sum(den)

  call sfft2_r2c(den,nc,1)
  write(*,*) 'called  fft forward'
  write(*,*)'Got Gaussian white density', minval(den),maxval(den),sum(den)

!  factor=0
!  do k=1,nk
!     factor(k)=pk(k)/(box/real(nc))**2
!     !factor(k)=pk(k)/(box**2/nc**2)
!  enddo
  kpt=nc/nt
  !!$omp parallel do default(shared) &
  !!$omp& private(it,i,j,k,kr,kx,ky,k1,k2,w1,w2)
  do it=1,nt
     do j=1+(it-1)*kpt,min(nc,it*kpt)
        if (j .le. hc+1) then
            ky=j-1
        else
            ky=j-1-nc
        endif
        do i=1,nc+2,2
           kx=(i-1)/2
           kr=sqrt(kx**2.+ky**2)
           if(kr.ne.0) then
              k1=ceiling(kr)
              k2=k1+1
              w1=k1-kr
              w2=1-w1
              den(i,j)=den(i,j)*sqrt(factor(k1-1)*w1+factor(k2-1)*w2) !! interpolation causes problem
              den(i+1,j)=den(i+1,j)*sqrt(factor(k1-1)*w1+factor(k2-1)*w2)
           else
              den(i,j)=den(i,j)*sqrt(factor(0))
              den(i+1,j)=den(i+1,j)*sqrt(factor(0))
           endif
        enddo
     enddo
  enddo
  !!$omp end parallel do

  write(*,*)'Got Gaussian Fourier density', minval(den),maxval(den),sum(den)

  call sfft2_r2c(den,nc,-1)
  write(*,*) 'called  fft backward'
  write(*,*)'Got Gaussian density', minval(den),maxval(den),sum(den)

  return
end subroutine GaussRandomField_2d_r2c

Subroutine GaussRandomField_2d_c2c(den, box, nc, kg, delta2, nk)
  implicit none

  !! nc is the number of grids on each dimension
  !! nt is the number of threads for openmp
  !! nk is the number of points at the input power spectrum table
  !! box is the physical size of the box
  integer nc,nk
  integer, parameter :: nt=1

!  real, dimension(0:nk) :: factor
  real, dimension(0:nc) :: factor
  real, dimension(nk) :: pk,kg,delta2
  real, dimension(nc,nc) :: x1,x2
  real, dimension(nc,nc) :: den
  complex, dimension(nc,nc) :: cd

  integer hc,it,i,j,k,k1,k2,kx,ky,kpt
  real box,pi,kr,w1,w2,tmp

  external sfft2_c2c

  pi=acos(-1.)
  hc=nc/2

  do i=1,nk
     pk(i)=delta2(i)/kg(i)**2*2*pi
  enddo

  !! generate whitenoise source
  call random_number(x1(:,:))
  call random_number(x2(:,:))

  den=0
  do j=1,nc
     do i=1,nc
        den(i,j)=sqrt(-2*log(x1(i,j)))*cos(2*pi*x2(i,j))
     enddo
  enddo

  cd=den+(0.,0.)
  call sfft2_c2c(cd,nc,1)
!  write(*,*) 'called  fft forward'

  factor=0
  do k=1,nk
     factor(k)=pk(k)/(box**2/nc**2)
  enddo

   kpt=nc/nt
  !!$omp parallel do default(shared) &
  !!$omp& private(it,i,j,k,kr,kx,ky,k1,k2,w1,w2)
  do it=1,nt
     do j=1+(it-1)*kpt,min(nc,it*kpt)
        if (j .le. hc+1) then
           ky=j-1
        else
           ky=j-1-nc
        endif
        do i=1,nc
           if (i .le. hc+1) then
              kx=i-1
           else
              kx=i-1-nc
           endif
           kr=sqrt(kx**2.+ky**2)
           tmp=nc**1.0/sqrt(cd(i,j)*conjg(cd(i,j)))
           if(kr.ne.0) then
              k1=ceiling(kr)
              k2=k1+1
              w1=k1-kr
              w2=1-w1
              cd(i,j)=cd(i,j)*sqrt(factor(k1-1)*w1+factor(k2-1)*w2) !! interpolation causes problem
           else
              cd(i,j)=cd(i,j)*sqrt(factor(0))
           endif
           !cd(i,j)=cd(i,j)*tmp
        enddo
     enddo
  enddo
  !!$omp end parallel do

  call sfft2_c2c(cd,nc,-1)
!  write(*,*) 'called  fft backward'
  den=real(cd)

  return
end subroutine GaussRandomField_2d_c2c
