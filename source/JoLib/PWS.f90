
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
        ps_err(i)=sqrt(2/weight(i))*ps(i)!Gaussian Error...
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

subroutine ps2_deconvolve(a,b,ps,n)
  implicit none
  integer n,i,j,i1,j1,ir
  real(kind=4) a(n,n),b(n,n),ps(n),w,r,weight(n),tmp,x,y,sinc_x,sinc_y,kernel,ell_ave(n)
  complex(kind=4) fa(n,n),fb(n,n)

  external sfft2_c2c

  fa=a+(0.,0.)
  fb=b+(0.,0.)
  call sfft2_c2c(fa,n,1)
  !call sfft2_c2c(fb,n,1)
  write(*,*) 'Assuming auto-correlation'
  fb = fa

  
  fa=fa/n**2
  fb=fb/n**2
  weight=0
  ps=0
  do j=1,n
     j1=j-1
     if(j1>n/2) j1=j1-n
     do i=1,n
        i1=i-1
        if(i1>n/2) i1=i1-n
        r=(i1**2+j1**2)**0.5
        !ir=ceiling(r)
        ir=r
        if(ir>=1) then

          x = 3.1415927*real(i1)/real(n)
          y = 3.1415927*real(j1)/real(n)

          if(x.eq.0) then
             sinc_x = 1
          else
             sinc_x = sin(x)/x
          endif
          if(y.eq.0) then
             sinc_y = 1
          else
             sinc_y = sin(y)/y
          endif

          kernel = sinc_x*sinc_y

          tmp=fa(i,j)*conjg(fb(i,j))/kernel**4

          w=r-ir ! CIC
          !w = 0.0!   ! NGP

          !if(ir .le. 3) write(*,*) i, i1, j, j1, r, ir

          ps(ir)=ps(ir)+tmp*(1.0-w)
          ps(ir+1)=ps(ir+1)+tmp*w
          weight(ir)=weight(ir)+1.0-w
          weight(ir+1)=weight(ir+1)+w
          ell_ave(ir)=ell_ave(ir)+(1.0-w)*r
          ell_ave(ir+1)=ell_ave(ir+1)+w*r
        end if
     end do
  end do
  do i=1,n
     if(weight(i)>0) then
        ps(i)=ps(i)/weight(i)*2*3.1415927*i**2
        ell_ave(i)= ell_ave(i)/weight(i)
        !write(*,*)i,weight(i),ell_ave(i)
     endif
  end do

  return
end subroutine ps2_deconvolve

!***********************************************************
!Get angular averaged 3D Delta^2(k) from a density matrix d(nc+2,nc,nc) 
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
              if(kx.eq.0 .and. ky <=0 .and. kz <=0) cycle;!write(*,*)'NO!!!';continue
              if(kx.eq.0 .and. ky >0 .and. kz <0) cycle;!write(*,*)'NO!!!';continue
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
     !write(*,*) weight(k+1)
     if (weight(k+1) .ne. 0) then
         ps(k)=4*pi*(k**3)*ps(k+1)/weight(k+1)
     endif
  enddo

!  write(*,*) 'Delta^2(k) : ', maxval(ps), ' > D > ',minval(ps)

  return
end subroutine ps3_r2c


subroutine ps3_r2c_err(d,ps,err,nc)
  implicit none

  integer nc,hc
  integer, parameter :: nt=4

  real, dimension(nc) :: ps,weight,err
  real, dimension(nc+2,nc,nc) :: d
  real, dimension(nc,nt) :: pst,weightt,errt

  integer it,i,j,k,kpt
  real kr,kx,ky,kz,k1,k2,w1,w2,pow,pi

  pi=acos(-1.)
  hc=nc/2
  kpt=nc/nt

  call sfft3_r2c(d,nc,1)

  !$omp parallel do default(shared) &
  !$omp& private(it,i,j,k,kr,kx,ky,kz,k1,k2,w1,w2,pow)
  do it=1,nt
     errt(:,it)=0
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
              if(kx.eq.0 .and. ky <=0 .and. kz <=0) cycle;!write(*,*)'NO!!!';continue
              if(kx.eq.0 .and. ky >0 .and. kz <0) cycle;!write(*,*)'NO!!!';continue
              if (kr .ne. 0) then
                  k1=ceiling(kr)
                  k2=k1+1
                  w1=k1-kr
                  w2=1-w1
                  pow=(d(i,j,k)/nc**3)**2+(d(i+1,j,k)/nc**3)**2
                  weightt(k1,it)=weightt(k1,it)+w1
                  pst(k1,it)=pst(k1,it)+w1*pow
                  errt(k1,it)=errt(k1,it)+w1*pow**2
                  weightt(k2,it)=weightt(k2,it)+w2
                  pst(k2,it)=pst(k2,it)+w2*pow
                  errt(k2,it)=errt(k2,it)+w2*pow**2
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
     err=err+errt(:,it)
     weight=weight+weightt(:,it)
  enddo

  !! Divide by weights
  !! Convert P(k) to \Delta^2(k)
  !! Store k in ps(1,i)

  do k=1,nc-1
     !write(*,*) weight(k+1)
     if (weight(k+1) .ne. 0) then
        ps(k) = ps(k+1)/weight(k+1)
        err(k)= err(k+1)/weight(k+1) - ps(k)**2       
        ps(k)=4*pi*(k**3)*ps(k)
        err(k)= 4*pi*(k**3)*err(k)
     endif
     !write(*,*) ps(k) ,err(k)
  enddo

!  write(*,*) 'Delta^2(k) : ', maxval(ps), ' > D > ',minval(ps)

  return
end subroutine ps3_r2c_err







!*********************************

subroutine ps3_c2c(d,ps,nc)
!subroutine ps3_c2c(d,ps,nc,psbin)
  implicit none
!#define logcic
!#define linearcic
!#define ksquare

  integer nc,hc
  integer, parameter :: nt=1

  real, dimension(nc,nc,nc) :: d
  complex, dimension(nc,nc,nc) :: cd

  external index_k2_3d

#ifdef ksquare
  real*8, dimension(3*nc**2/4+1) :: ps,weight
  real*8, dimension(3*nc**2/4+1,nt) :: pst,weightt
  integer, dimension(2,3*nc**2/4) :: index_shell
  integer, dimension(3*nc**2/4) :: arr_shell 
#else
  real, dimension(2*nc,nt) :: pst,weightt
  real, dimension(nc) :: ps,weight
#endif

  integer, parameter :: nbin=5
  real*8, dimension(nbin) :: psbin,logpsbin,countbin
  real*8  logkmax,logkmin,dlogk 
  integer ibin

  integer it,i,j,k,kpt,indx,k1,k2,kx,ky,kz
  real*8 w1,w2,pi
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

  logkmin=log(2*pi/300)!-3.865905       
  logkmax=0.2398363 !1.626555       
  dlogk=0.821148335933685 !1.09849202632904                

!  write(*,*) log(2*pi/300*nc)
!  pause

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
!                 !! TEST BIN
!                 ibin=floor((log(kr*2*pi/300)-logkmin)/dlogk)+1
!                 if(ibin.le.0)then
!                   write(*,*) kx**2+ky**2+kz**2, log(2*pi/300*1), logkmin, log(kr*2*pi/300)
!                   pause
!                 endif
!                 if(ibin.eq.5)then 
!                   write(*,*) 'ibin=5', index_shell(2,kx**2+ky**2+kz**2)-3*hc**2/4, log(kr*2*pi/300)   
!                   pause
!                 endif
!                 if (index_shell(2,kx**2+ky**2+kz**2).le.3*hc**2/4) then 
!                   logpsbin(ibin)=logpsbin(ibin)+log(pow)*index_shell(1,kx**2+ky**2+kz**2)
!                   countbin(ibin)=countbin(ibin)+index_shell(1,kx**2+ky**2+kz**2)
!                 endif
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

!  do ibin=1,nbin
!     write(*,*) ibin, countbin(ibin),logpsbin(ibin)
!     if(countbin(ibin).gt.0)then
!         logpsbin(ibin)=logpsbin(ibin)/countbin(ibin)
!         psbin(ibin)=exp(logpsbin(ibin)) 
!     endif
!  enddo
!  psbin=psbin*300.**3
  

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
  
#ifdef ksquare
  do k=1,3*nc**2/4
     if (weight(k+1) .ne. 0) then
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


!! Modified on Hy Trac's code

subroutine ps2_camb(a,b,ps, ps_camb,n)
  implicit none
  integer n,i,j,i1,j1,ir
  real(kind=4) a(n,n),b(n,n),ps(n),w,r,weight(n),tmp,x,y,sinc_x,sinc_y,kernel
  complex(kind=4) fa(n,n),fb(n,n)
  integer, parameter :: nk = 6000
  real ell_camb(nk), C_ell_camb(nk), dell_low, dell_high, d2ydx2(nk), tmp_camb, ps_camb(n), angle
   

  external sfft2_c2c
  
  angle = sqrt(60.)*3.14159/180.0


  fa=a+(0.,0.)
  fb=b+(0.,0.)
  call sfft2_c2c(fa,n,1)
  !call sfft2_c2c(fb,n,1)
  write(*,*) 'Assuming auto-correlation'
  fb = fa

  
  open(11,file='./elCl_par_WMAP9_SN_BAO_z0_z0.582_nokcut.dat', status = 'old')  
  do i=1,nk
     read(11,*) ell_camb(i), C_ell_camb(i)
  enddo
  close(11)
  
  dell_low  = (C_ell_camb(2) - C_ell_camb(1))/(ell_camb(2) - ell_camb(1))
  dell_high = (C_ell_camb(nk) - C_ell_camb(nk-1))/(ell_camb(nk) - ell_camb(nk-1))     

  call spline(ell_camb, C_ell_camb, nk, dell_low,dell_high,d2ydx2)

  fa=fa/n**2
  fb=fb/n**2
  weight=0
  ps=0
  do j=1,n
     j1=j-1
     if(j1>n/2) j1=j1-n
     do i=1,n
        i1=i-1
        if(i1>n/2) i1=i1-n
        r=(i1**2+j1**2)**0.5
        ir=r
        if(ir>=1) then

          x = 3.1415927*real(i1)/real(n)
          y = 3.1415927*real(j1)/real(n)

          if(x.eq.0) then
             sinc_x = 1
          else
             sinc_x = sin(x)/x
          endif
          if(y.eq.0) then
             sinc_y = 1
          else
             sinc_y = sin(y)/y
          endif

          kernel = sinc_x*sinc_y

          call splint(ell_camb, C_ell_camb, d2ydx2, nk, r*(2*3.14159/angle -1.0), tmp_camb)

          tmp=fa(i,j)*conjg(fb(i,j))/kernel**4
          w=r-ir
          ps(ir)=ps(ir)+tmp*(1-w)
          ps(ir+1)=ps(ir+1)+tmp*w
          ps_camb(ir)=ps_camb(ir)+tmp_camb*(1-w)
          ps_camb(ir+1)=ps_camb(ir+1)+tmp_camb*w
          weight(ir)=weight(ir)+1-w
          weight(ir+1)=weight(ir+1)+w
        end if
     end do
  end do
  do i=1,n
     if(weight(i)>0) ps(i)=ps(i)/weight(i)*2*3.1415927*i**2
     if(weight(i)>0) ps_camb(i)=ps_camb(i)/weight(i)*2*3.1415927*i**2
  end do

  return
end subroutine ps2_camb
