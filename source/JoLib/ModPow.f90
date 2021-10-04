
subroutine ModPow_3d_r2c(den, box, nc, kg, delta2, nk) 
  implicit none
  
  integer nc,nk 
  integer, parameter :: nt=1

  real, dimension(nc+2,nc,nc) :: den
  real, dimension(nk) :: pk,kg,delta2
  real, dimension(0:nk) :: factor
  real, dimension(2,nc) :: ps
  real, dimension(2,nc,nt) :: pst
  real box

  integer hc,it,i,j,k,k1,k2,kx,ky,kz,kpt
  real pi,kr,w1,w2,pow

  external linear_interpolation, sfft3_r2c
  
  pi=acos(-1.)
  hc=nc/2

  do i=1,nk
     pk(i)=delta2(i)/kg(i)**3*2*pi**2
  enddo

  call sfft3_r2c(den,nc,1)
  
  kpt=nc/nt  
  !!$omp parallel do default(shared) &
  !!$omp& private(it,i,j,k,kr,kx,ky,kz,k1,k2,w1,w2,pow)
  do it=1,nt
     pst(:,:,it)=0
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
              kr=sqrt(kx**2*1.+ky**2+kz**2)
              if (kr .ne. 0) then
                  k1=ceiling(kr)
                  k2=k1+1
                  w1=k1-kr
                  w2=1-w1
                  pow=(den(i,j,k)/nc**3)**2+(den(i+1,j,k)/nc**3)**2
                  pst(1,k1,it)=pst(1,k1,it)+w1
                  pst(2,k1,it)=pst(2,k1,it)+w1*pow
                  pst(1,k2,it)=pst(1,k2,it)+w2
                  pst(2,k2,it)=pst(2,k2,it)+w2*pow
              endif
           enddo
        enddo
     enddo
  enddo
  !$omp end parallel do
   
  ps=0
  do it=1,nt
     ps=ps+pst(:,:,it)
  enddo
  
  factor=0
  do k=1,nk-1
     if (ps(1,k+1) .ne. 0) then
         !write(*,*) ps(1,k+1),ps(2,k+1) 
         ps(2,k+1)=4*pi*(k**3)*ps(2,k+1)/ps(1,k+1)
         ps(1,k+1)=2*pi*k/box
         call linear_interpolation(kg,delta2,nk,ps(1,k+1),factor(k))
         factor(k)=factor(k)/ps(2,k+1) !! note it's nonlinear to use the measured power ??? 
     endif
  enddo
  write(*,*) maxval(factor),minval(factor)
  
  !$omp parallel do default(shared) &
  !$omp& private(it,i,j,k,kr,kx,ky,kz,k1,k2,w1,w2)
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
                kr=sqrt(kx**2*1.+ky**2+kz**2)
                if (kr .ne. 0) then
                    k1=ceiling(kr)
                    k2=k1+1
                    w1=k1-kr
                    w2=1-w1
                    den(i,j,k)=den(i,j,k)*sqrt(factor(k1)*w1+factor(k2)*w2)
                    den(i+1,j,k)=den(i+1,j,k)*sqrt(factor(k1)*w1+factor(k2)*w2)
                endif
              enddo
          enddo
       enddo
  enddo
  !$omp end parallel do
       
  call sfft3_r2c(den,nc,-1)
  
  return
end subroutine ModPow_3d_r2c
