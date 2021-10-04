
!! Modified on Pengjie Zhang's code

subroutine ps2(a,b,ps,n)
  implicit none
  integer n,i,j,i1,j1,ir
  real(kind=4) a(n,n),b(n,n),ps(n),w,r,weight(n),tmp
  complex(kind=4) fa(n,n),fb(n,n)

  external sfft2_c2c

  fa=a+(0.,0.)
  fb=b+(0.,0.)
  call sfft2_c2c(fa,n,1)
  call sfft2_c2c(fb,n,1)

  fa=fa/n**2
  fb=fb/n**2
  weight=0
  ps=0
  do j=1,n
     j1=j-1
!     if(j>n/2) j1=j1-n
     if(j1>n/2) j1=j1-n
     do i=1,n
        i1=i-1
!        if(i>n/2) i1=i1-n
        if(i1>n/2) i1=i1-n
        r=(i1**2+j1**2+1./3)**0.5
        ir=r
        if(ir>=1) then
          tmp=fa(i,j)*conjg(fb(i,j))
          w=r-ir
          ps(ir)=ps(ir)+tmp*(1-w)
          ps(ir+1)=ps(ir+1)+tmp*w
          weight(ir)=weight(ir)+1-w
          weight(ir+1)=weight(ir+1)+w
        end if
     end do
  end do
  do i=1,n
     if(weight(i)>0) ps(i)=ps(i)/weight(i)*2*3.1415927*i**2
  end do

  return
end subroutine ps2

!! Modified on Hy Trac's code

subroutine ps3_r2c(d,ps,nc)
  implicit none

  integer nc,hc
  integer, parameter :: nt=4

  real, dimension(nc) :: ps,weight
  real, dimension(nc+2,nc,nc) :: d
  real, dimension(nc,nt) :: pst,weightt

  integer it,i,j,k,kpt
  real kr,kx,ky,kz,k1,k2,w1,w2,pow,pi

  pi=acos(-1.)
  hc=nc/2
  kpt=nc/nt

!  write(*,*) 'd overdensity, before fft', maxval(d),minval(d)
!  pause
  call sfft3_r2c(d,nc,1)
!  write(*,*) 'd overdensity, after fft', maxval(d),minval(d)
!  pause

  !$omp parallel do default(shared) &
  !$omp& private(it,i,j,k,kr,kx,ky,kz,k1,k2,w1,w2,pow)
  do it=1,nt
     pst(:,it)=0
     weightt(:,it)=0
     do k=1+(it-1)*kpt,min(nc,it*kpt)
        if (k .lt. hc+2) then
            kz=k-1
        else
            kz=k-1-nc
        endif
        do j=1,nc
           if (j .lt. hc+2) then
               ky=j-1
           else
               ky=j-1-nc
           endif
           do i=1,nc+2,2
              kx=(i-1)/2
              kr=sqrt(kx**2+ky**2+kz**2)
              if (kr .ne. 0) then
                  k1=ceiling(kr)
                  k2=k1+1
                  w1=k1-kr
                  w2=1-w1
                  pow=(d(i,j,k)/nc**3)**2+(d(i+1,j,k)/nc**3)**2
!                  write(*,*) 'kr,pow', kr,pow
!                  pause   
                  weightt(k1,it)=weightt(k1,it)+w1
                  pst(k1,it)=pst(k1,it)+w1*pow
                  weightt(k2,it)=weightt(k2,it)+w2
                  pst(k2,it)=pst(k2,it)+w2*pow
              endif
           enddo
        enddo
     enddo
  enddo
  !$omp end parallel do

  
  !! Merge power spectrum from threads
  ps=0
  weight=0
  do it=1,nt
     ps=ps+pst(:,it)
     weight=weight+weightt(:,it)
  enddo

  !! Divide by weights
  !! Convert P(k) to \Delta^2(k)
  !! Store k in ps(1,i)

  do k=1,nc-1
     if (weight(k+1) .ne. 0) then
         ps(k)=4*pi*(k**3)*ps(k+1)/weight(k+1)
     endif
  enddo

  return
end subroutine ps3_r2c

