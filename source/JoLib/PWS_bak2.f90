
!! Modified on Pengjie Zhang's code

subroutine ps2_x(a,b,ps,n)
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
  do j=1,1
     do i=1,n
        i1=i-1
        if(i1>n/2) i1=i1-n
        !ir=abs(i1)  
        ir=i1
        if(ir>=1) then
          tmp=fa(i,j)*conjg(fb(i,j))
          !if(ir.eq.1)write(*,*) 'after fft', fa(i,j),fa(i,j)*conjg(fb(i,j))  
          ps(ir)=ps(ir)+tmp
          weight(ir)=weight(ir)+1
        endif
     enddo
  enddo

  do i=1,n
     if(weight(i)>0) ps(i)=ps(i)/weight(i)*2*3.1415927*i**2
  enddo
  !pause

  return
end subroutine ps2_x

!! err bar, when different modes are independent
subroutine ps2_errbar(a,b,ps,ps_err,n)
  implicit none
  integer n,i,j,i1,j1,ir
  real(kind=4) a(n,n),b(n,n),ps(n),ps_S2(n),ps_err(n),w,r,weight(n),tmp,pi
  complex(kind=4) fa(n,n),fb(n,n)

  external sfft2_c2c

  pi=acos(-1.)
  fa=a+(0.,0.)
  fb=b+(0.,0.)
  call sfft2_c2c(fa,n,1)
  call sfft2_c2c(fb,n,1)
  
  fa=fa/n**2
  fb=fb/n**2
  weight=0
  ps=0
  ps_err=0
  do j=1,n
     j1=j-1
     if(j1>n/2) j1=j1-n
     do i=1,n
        i1=i-1
        if(i1>n/2) i1=i1-n
        r=(i1**2+j1**2)**0.5
        ir=r
        if(ir>=1) then
          tmp=fa(i,j)*conjg(fb(i,j))
          w=r-ir
          ps(ir)=ps(ir)+tmp*(1-w)
          ps(ir+1)=ps(ir+1)+tmp*w
          ps_S2(ir)=ps_S2(ir)+tmp**2*(1-w)
          ps_S2(ir+1)=ps_S2(ir+1)+tmp**2*w
          weight(ir)=weight(ir)+1-w
          weight(ir+1)=weight(ir+1)+w
        endif
     enddo
  enddo
  do i=1,n
     if(weight(i)>0)then
        ps(i)=ps(i)/weight(i) 
        ps_S2(i)=ps_S2(i)/weight(i)
        ps_err(i)=sqrt(2/weight(i))*ps(i)
        !write(*,*) sigma(i), sqrt(ps_S2(i)-ps(i)**2),ps_S2(i),ps(i) 
        ps(i)=ps(i)*2*pi*i**2
        ps_err(i)=ps_err(i)*2*pi*i**2
     endif 
  enddo

  return
end subroutine ps2_errbar

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
        !r=(i1**2+j1**2+1./3)**0.5
        r=(i1**2+j1**2)**0.5
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

  call sfft3_r2c(d,nc,1)

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

subroutine ps3_c2c(d,ps,nc)
  implicit none
!#define logcic
!#define linearcic
#define ksquare

  integer nc,hc
  integer, parameter :: nt=1

  real, dimension(nc,nc,nc) :: d
  complex, dimension(nc,nc,nc) :: cd

  external index_k2_3d

#ifdef ksquare
  real, dimension(3*nc**2/4+1) :: ps,weight
  real*8, dimension(3*nc**2/4+1,nt) :: pst,weightt
  integer, dimension(2,3*nc**2/4) :: index_shell
  integer, dimension(3*nc**2/4) :: arr_shell 
#else
  real*8, dimension(2*nc,nt) :: pst,weightt
  real, dimension(nc) :: ps,weight
#endif

  integer it,i,j,k,kpt,indx,k1,k2,kx,ky,kz
  real w1,w2,pi
  real*8 pow,kr

  pi=acos(-1.)
  hc=nc/2
  kpt=nc/nt

  cd=d+(0.,0.)
  !write(*,*) 'sum(d(**2)= (before FFT)', sum(real(cd*conjg(cd),8)/nc/nc/nc)
  
  call sfft3_c2c(cd,nc,1)
     
  !write(*,*) cd
  !write(*,*) cd(2,1,1)*conjg(cd(2,1,1))/nc**6 
  !pause
 
 
  !!write(*,*) 'sum(d(**2)= (after FFT)', sum(real(cd*conjg(cd)/nc/nc/nc,8)/nc/nc/nc)
  !write(*,*) 'sum(d(**2)= (after FFT)', sum(real(cd*conjg(cd)/nc/nc/nc,4)/nc/nc/nc)

#ifdef ksquare
  call index_k2_3d(index_shell,arr_shell,nc)
#endif

  !!$omp parallel do default(shared) &
  !!$omp& private(it,i,j,k,kr,kx,ky,kz,k1,k2,w1,w2,pow,indx)
  do it=1,nt
     pst(:,it)=0
     weightt(:,it)=0
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
               kr=sqrt(kx**2.+ky**2+kz**2)
               if(kr.ne.0) then
#ifdef ksquare
                 indx=index_shell(2,kx**2+ky**2+kz**2)
                 !arr_shell(indx)=kx**2+ky**2+kz**2 
                 !indx=kx**2+ky**2+kz**2
                 pow=(cd(i,j,k)/nc/nc/nc)*conjg(cd(i,j,k)/nc/nc/nc)
                 pst(indx+1,it)=pst(indx+1,it)+pow
                 weightt(indx+1,it)=weightt(indx+1,it)+1
                 !if(abs(kx)+abs(ky)+abs(kz).eq.1)then
                 !   write(*,*) kx,ky,kz, pow
                 !   pause
                 !endif
#else
                 k1=ceiling(kr)
                 k2=k1+1
                 !! cic in linear scale
#ifdef linearcic
                 w1=k1-kr
                 w2=1-w1
                 pow=(cd(i,j,k)/nc**3)*conjg(cd(i,j,k)/nc**3)
                 !! cic in log scale
#endif
#ifdef logcic
                 pow=(cd(i,j,k)/nc**3)*conjg(cd(i,j,k)/nc**3)
                 w1=0
                 w2=0
                 if(pow.eq.0)then
                    w1=0
                    w2=0 
                 elseif(k1.gt.1)then
                    w1=(log(k2-1.)-log(kr+0.))/(log(k2-1.)-log(k1-1.))
                    w2=1-w1
                    pow=log(pow)
                 elseif(kr.eq.1)then
                    w1=1
                    w2=0
                    pow=log(pow)
                 endif
#endif
                 !! ----------------------
                 weightt(k1,it)=weightt(k1,it)+w1
                 pst(k1,it)=pst(k1,it)+w1*pow
                 weightt(k2,it)=weightt(k2,it)+w2
                 pst(k2,it)=pst(k2,it)+w2*pow
#endif
              endif
           enddo
        enddo
     enddo
  enddo
  !!$omp end parallel do

  !! Merge power spectrum from threads
  ps=0
  weight=0
  do it=1,nt
     ps=ps+pst(:,it)
     weight=weight+weightt(:,it)
     !write(*,*) ps(10),pst(10,it),pst(10,it)/weightt(10,it),ps(10)/weight(10)
     !pause  
  enddo

  !! Divide by weights
  !! Convert P(k) to \Delta^2(k)
  !! Store k in ps(1,i)

#ifdef ksquare
  do k=1,3*nc**2/4
     if (weight(k+1) .ne. 0) then
         !k1=int(k/3)
         !k2=mod(k,3)
         !if(k2.eq.0)then
         !   kr=sqrt(3.*k1**2)
         !elseif(k2.eq.1)then
         !   kr=sqrt(2.*k1**2+(k1+1)**2) 
         !else
         !   kr=sqrt(k1**2+2.*(k1+1)**2) 
         !endif  
         kr=arr_shell(k)
         ps(k)=4*pi*(kr**1.5)*ps(k+1)/weight(k+1)
     endif   
  enddo
#else
  do k=1,nc-1
     if (weight(k+1) .ne. 0) then
         !write(*,*) k, ps(k+1), weight(k+1), exp(ps(k+1)/weight(k+1)) 
         !pause
#ifdef linearcic
         ps(k)=4*pi*(k**3)*ps(k+1)/weight(k+1)
#endif
#ifdef logcic
         ps(k)=4*pi*(k**3)*exp(ps(k+1)/weight(k+1))
#endif
     else
         ps(k)=0
     endif
     !write(*,*) k,k3**2,ps(k)
     !pause
  enddo
#endif

  return
end subroutine ps3_c2c
