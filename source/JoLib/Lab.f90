!*********************************

!*****************************

! Angle-average a vector density
! in K-Space, using FFT ordering:
! [ky and kz] -> [0-hc,-hc+1,-1]
! [kx] -> odd  = Re[1-hc+1]
! [kx] -> even = Im[1-hc+1]


subroutine KAngleAverage(Vect,Ave,nc)

  implicit none

  integer nc,hc,space
!  integer, parameter :: nt=4

  real, dimension(nc) :: Ave,weight
  real, dimension(nc+2,nc,nc) :: Vect
!  real, dimension(nc,nt) :: Ave_t,weightt

  integer i,j,k, Nextreme!,it,kpt
  real kr,kx,ky,kz,k1,k2,w1,w2,Mag,pi, MAX

  !pi=acos(-1.)
  hc=nc/2
  !kpt=nc/nt
  
  MAX = 1000000000 !10^19

  Ave(:) = 0
 
  do k=1,nc
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
              
              !Mag=abs(Vect(i,j,k))
              Mag=sqrt(sum(Vect(i:i+1,j,k)**2))
              weight(k1)=weight(k1)+w1 ! counting how many times this magnitude was involved, for spatial averaging
              Ave(k1)=Ave(k1)+w1*Mag
              weight(k2)=weight(k2)+w2
              Ave(k2)=Ave(k2)+w2*Mag
              
              
              !Checks
              if(Ave(k1) > MAX)then
                 
                 write(*,*) 'Extremum at Coordinates = (' ,i,j,k, ')'
                 write(*,*) 'Vect = ', Vect(i,j,k),'Mag = ', Mag 
                 write(*,*) 'Ave(k1) = ', Ave(k1), 'w1 = ' , w1, 'w2 = ', w2
                 write(*,*) 'setting Ave(k1) to 0...'
                 
                 Nextreme = Nextreme + 1
                 Ave(k1) = 0
                 pause
              endif
              
              if(Ave(k2) > MAX)then
                 
                 write(*,*) 'Extremum at Coordinates = (' ,i,j,k, ')'
                 write(*,*) 'Vect = ', Vect(i,j,k),'Mag = ', Mag 
                 write(*,*) 'Ave(k2) = ', Ave(k2), 'w1 = ' , w1, 'w2 = ', w2
                 write(*,*) 'setting Ave(k2) to 0...'
                 
                 Nextreme = Nextreme + 1
                 write(*,*) 'Ave(k2) = ', Ave(k2), 'w1 = ' , w1, 'w2 = ', w2
                 write(*,*) 'setting Ave(k2) to 0...'
                 
                 Nextreme = Nextreme + 1
                 Ave(k2) = 0
                 pause
              endif
              
              !write(*,*) 'Coordinates = (' ,kx,ky,kz, ') yield kr = ', kr 
              !write(*,*) 'Vect = ', Vect(i,j,k),'Mag = ', Mag 
              !write(*,*) 'Ave(k1) = ', Ave(k1), 'w1 = ' , w1
              !write(*,*) 'Ave(k2) = ', Ave(k2), 'w2 = ', w2
              !pause
              
           endif
           
        enddo
     enddo
  enddo
  
  
  !write(*,*) 'Averaging Ave'
  
  do k = 1,nc-1
     if(weight(k+1) .ne. 0) then
        Ave(k) = Ave(k+1)/weight(k+1)
     endif
  enddo
  
  write(*,*) 'Averaged Covariance', maxval(Ave),' > CovAve > ',minval(Ave) 
  
  write(*,*) 'Nextreme = ', Nextreme
  
  return
end subroutine KAngleAverage


!*****************************

! Angle-average a vector density
! in X-Space, [x,y,z] -> [1-nc,1-nc,1-nc] 

subroutine XAngleAverage(Vect,Ave,nc)

  implicit none

  integer nc,hc,space
!  integer, parameter :: nt=4

  real, dimension(nc) :: Ave,weight
  real, dimension(nc+2,nc,nc) :: Vect
!  real, dimension(nc,nt) :: Ave_t,weightt

  integer i,j,k, Nextreme!,it,kpt
  real kr,kx,ky,kz,k1,k2,w1,w2,Mag,pi, MAX

  !pi=acos(-1.)
  hc=nc/2
  !kpt=nc/nt
  
  MAX = 1000000000 !10^19

  Ave(:) = 0
  
  write(*,*) 'Entering Average loop...'
  
  do k=1,hc! SHOULDN'T IT BE 1-nc????
     kz = k
     
     do j=1,hc
        ky = j
        
        do i=1,hc
                
           kx = i     
           
           kr=sqrt(kx**2+ky**2+kz**2)
           
           !write(*,*) 'Coordinates = (' ,kx,ky,kz, ') yield kr = ', kr 


           if (kr .ne. 0) then
              k1=ceiling(kr)
              k2=k1+1
              w1=k1-kr
              w2=1-w1
              
              !Mag=abs(Vect(i,j,k))
              Mag=Vect(i,j,k)
              weight(k1)=weight(k1)+w1 ! counting how many times this magnitude was involved, for spatial averaging
              Ave(k1)=Ave(k1)+w1*Mag
              weight(k2)=weight(k2)+w2
              Ave(k2)=Ave(k2)+w2*Mag
              
              
              !Checks
              if(Ave(k1) > MAX)then
                 
                 write(*,*) 'Extremum at Coordinates = (' ,i,j,k, ')'
                 write(*,*) 'Vect = ', Vect(i,j,k),'Mag = ', Mag 
                 write(*,*) 'Ave(k1) = ', Ave(k1), 'w1 = ' , w1, 'w2 = ', w2
                 write(*,*) 'setting Ave(k1) to 0...'
                 
                 Nextreme = Nextreme + 1
                 Ave(k1) = 0
                 pause
              endif
              
              if(Ave(k2) > MAX)then
                 
                 write(*,*) 'Extremum at Coordinates = (' ,i,j,k, ')'
                 write(*,*) 'Vect = ', Vect(i,j,k),'Mag = ', Mag 
                 write(*,*) 'Ave(k2) = ', Ave(k2), 'w1 = ' , w1, 'w2 = ', w2
                 write(*,*) 'setting Ave(k2) to 0...'
                 
                 Nextreme = Nextreme + 1
                 write(*,*) 'Ave(k2) = ', Ave(k2), 'w1 = ' , w1, 'w2 = ', w2
                 write(*,*) 'setting Ave(k2) to 0...'
                 
                 Nextreme = Nextreme + 1
                 Ave(k2) = 0
                 pause
              endif
              
              !write(*,*) 'Coordinates = (' ,kx,ky,kz, ') yield kr = ', kr 
              !write(*,*) 'Vect = ', Vect(i,j,k),'Mag = ', Mag 
              !write(*,*) 'Ave(k1) = ', Ave(k1), 'w1 = ' , w1
              !write(*,*) 'Ave(k2) = ', Ave(k2), 'w2 = ', w2
              !pause
              
           endif
           
        enddo
     enddo
  enddo
  

  
  !write(*,*) 'Averaging Ave'
  
  do k = 1,nc-1
     if(weight(k+1) .ne. 0)then
        Ave(k) = Ave(k+1)/weight(k+1)
     endif
  enddo
  
  write(*,*) 'Averaged Value: ', maxval(Ave),' > Ave(r) > ',minval(Ave) 
  
  write(*,*) 'Nextreme = ', Nextreme
  
  return
end subroutine XAngleAverage



!********************************

!! Get Covariance of power spectrum: Cov(dkx,dky,dkz) from power spectrum
!! d(i,j,k), 
!! where dkx, dky and dkz run from 1 to nc.

!! Don't forget that the real wave numbers run from -257 to 256...



subroutine Cov_c2r(d,Cov,nc)
  implicit none

  integer nc,hc
  integer, parameter :: nt=4

  !real, dimension(nc) :: ps,weight
  real, dimension(nc+2,nc,nc) :: d, Cov
  !real, dimension(nc,nt) :: pst,weightt
  real, dimension(nc+2,nc,nc, nt) :: Covt
  
  integer it,i,j,k,kpt,ii, Nextreme
  real kx,ky,kz,pi,pow,MAX !kr,k1,k2,w1,w2

  pi=acos(-1.)
  hc=nc/2
  kpt=nc/nt
  MAX = 100000000000


  write(*,*) 'd before fft : ', maxval(d),' > d > ' ,minval(d)

  !do ii = 1,nc
     !write(*,*) 'PPShell(', ii, ', hc, hc) = ', d(ii,hc,hc)  
  !enddo

  call sfft3_r2c(d,nc,-1) 

  write(*,*) 'd after fft : ', maxval(d),' > d > ' ,minval(d)

 
  ! Now d is in x-space.

  !$omp parallel do default(shared) &
  !$omp& private(it,i,j,k,kr,kx,ky,kz,k1,k2,w1,w2,pow)
  do it=1,nt
     Covt(:,:,:,it)=0
     !weightt(:,it)=0
     do k=1+(it-1)*kpt,min(nc,it*kpt)
        if (k .lt. hc+2) then
            kz=k-1
        else
            kz=k-1-nc
        endif
        ! I will shift values of ky and kz so that they run from 1 to nc
        kz = kz + hc
        
        do j=1,nc
           if (j .lt. hc+2) then
              ky=j-1
           else
              ky=j-1-nc
           endif
           ky = ky + hc

           do i=1,nc+2,2 !(i odd = real, i even = imaginary)
              kx=(i-1)/2 
              
              ! kx runs from 1 to hc. 
              ! I will shift it so that it runs from hc+1 to nc
 
              kx = kx + hc
                         
              !write(*,*) 'Test : kx = ', kx, 'ky = ', ky, 'kz = ',kz
              !pause

              if((kx==5*hc/4) .and. (ky==hc/4) .and. (kz==hc/4)) then

                 write(*,*) 'Re(PPshell) = ', d(i,j,k)
                 write(*,*) 'Im(PPshell) = ', d(i+1,j,k)
                 pause
              endif

              if((kx==5*hc/4) .and. (ky==3*hc/4) .and. (kz==3*hc/4)) then    
                 write(*,*) 'Re(PPshell) = ', d(i,j,k)
                 write(*,*) 'Im(PPshell) = ', d(i+1,j,k)
                 pause

              endif


              ! Consider only the real part:

              Covt(kx,ky,kz,it)=d(i,j,k)

              !kr=sqrt(kx**2+ky**2+kz**2)
              !if (kr .ne. 0) then
                  !k1=ceiling(kr)
                  !k2=k1+1
                  !w1=k1-kr
                  !w2=1-w1
                  !pow=(d(i,j,k)/nc**3)**2+(d(i+1,j,k)/nc**3)**2
                  !weightt(k1,it)=weightt(k1,it)+w1
                  !pst(k1,it)=pst(k1,it)+w1*pow
                  !weightt(k2,it)=weightt(k2,it)+w2
                  !pst(k2,it)=pst(k2,it)+w2*pow
              !endif
           enddo
        enddo
     enddo
  enddo
  !$omp end parallel do

  !! Recover the full nc^3 power spectrum by filling kx = 0 to hc
  !! and merge power spectrum from threads

  Cov=0
  !weight=0
  do it=1,nt
     do kz=1,nc-2  
        do ky=1,nc  
           do kx=1,nc
              if(kx < hc+1)then
                 Covt(kx,ky,kz,it)=Covt(nc-kx,nc-ky,nc-kz,it)
              endif
              Cov(kx,ky,kz)=Cov(kx,ky,kz)+Covt(kx,ky,kz,it)

              ! Place a maximum value to Power Spectrum...
              if( (Cov(kx,ky,kz) > MAX) .or. (Cov(kx,ky,kz) < -MAX) ) then
                 Cov(kx,ky,kz) = 0 ! = MAX ! ??
                 write(*,*) 'coordinates of extremum :', kx,ky,kz 
                 Nextreme = Nextreme + 1
              endif


           enddo
        enddo
     enddo
 
     !Cov=Cov+Covt(:,:,:,it)
     !weight=weight+weightt(:,it)
  enddo


  !do ii = 1,nc
     !write(*,*) 'Cov(', ii, ', 3hc/4, 3hc/4) = ', Cov(ii,3*hc/4,3*hc/4)  
  !enddo


  write(*,*) 'Test : Cov(-3*hc/4, 3*hc/4, 3*hc/4) = ', Cov(hc-3*hc/4,hc+3*hc/4,hc+3*hc/4)
  write(*,*) 'Test : Cov(3*hc/4, -3*hc/4, -3*hc/4) = ', Cov(hc+3*hc/4,hc-3*hc/4,hc-3*hc/4)
  write(*,*) 'Test : Cov(-hc/4, hc/4, hc/4) = ', Cov(hc-hc/4,hc+hc/4,hc+hc/4)
  write(*,*) 'Test : Cov(hc/4, -hc/4, -hc/4) = ', Cov(hc+hc/4,hc-hc/4,hc-hc/4)

  !This test was succesfull!!!

  write(*,*) 'N extreme = ', Nextreme
  write(*,*) 'Covariance : ', maxval(Cov), ' > Cov > ',minval(Cov)


  return
end subroutine Cov_c2r

!*********************************

subroutine CovAve_c2r(d,Cov,nc)
  implicit none

  integer nc,hc
  integer, parameter :: nt=4

  real, dimension(nc) :: weight, Cov
  real, dimension(nc+2,nc,nc) :: d !, Cov
  real, dimension(nc,nt) :: weightt, Covt
  !real, dimension(nc+2,nc,nc, nt) :: Covt
  
  integer it,i,j,k,kpt,ii, Nextreme
  real kx,ky,kz,pi,pow,MAX,kr,k1,k2,w1,w2

  pi=acos(-1.)
  hc=nc/2
  kpt=nc/nt
  MAX = 100000000000


  write(*,*) 'd before fft : ', maxval(d),' > d > ' ,minval(d)

  !do ii = 1,nc
     !write(*,*) 'PPShell(', ii, ', hc, hc) = ', d(ii,hc,hc)  
  !enddo

  call sfft3_r2c(d,nc,-1) 

  write(*,*) 'd after fft : ', maxval(d),' > d > ' ,minval(d)

 
  ! Now d is in x-space.

  !$omp parallel do default(shared) &
  !$omp& private(it,i,j,k,kr,kx,ky,kz,k1,k2,w1,w2,pow)
  do it=1,nt
     Covt(:,it)=0
     weightt(:,it)=0
     do k=1+(it-1)*kpt,min(nc,it*kpt)
        if (k .lt. hc+2) then
            kz=k-1
        else
            kz=k-1-nc
        endif
        ! I will shift values of ky and kz so that they run from 1 to nc
        !kz = kz + hc
        
        do j=1,nc
           if (j .lt. hc+2) then
              ky=j-1
           else
              ky=j-1-nc
           endif
           !ky = ky + hc

           do i=1,nc+2,2 !(i even = real, i odd = imaginary)
              kx=(i-1)/2 
              
              ! kx runs from 1 to hc. 
              ! I will shift it so that it runs from hc+1 to nc
 
              !kx = kx + hc
                         
              !1write(*,*) 'Test : kx = ', kx, 'ky = ', ky, 'kz = ',kz
              !pause

              if((kx==5*hc/4) .and. (ky==hc/4) .and. (kz==hc/4)) then

                 write(*,*) 'Re(PPshell) = ', d(i,j,k)
                 write(*,*) 'Im(PPshell) = ', d(i+1,j,k)
                 pause
              endif

              if((kx==5*hc/4) .and. (ky==3*hc/4) .and. (kz==3*hc/4)) then    
                 write(*,*) 'Re(PPshell) = ', d(i,j,k)
                 write(*,*) 'Im(PPshell) = ', d(i+1,j,k)
                 pause

              endif


              ! Consider only the real part:

              !Covt(kx,ky,kz,it)=d(i,j,k)

              kr=sqrt(kx**2+ky**2+kz**2)
              if (kr .ne. 0) then
                  k1=ceiling(kr)
                  k2=k1+1
                  w1=k1-kr
                  w2=1-w1
                  pow=(d(i,j,k)/nc**3)
                  weightt(k1,it)=weightt(k1,it)+w1
                  Covt(k1,it)=Covt(k1,it)+w1*pow
                  weightt(k2,it)=weightt(k2,it)+w2
                  Covt(k2,it)=Covt(k2,it)+w2*pow
              endif
           enddo
        enddo
     enddo
  enddo
  !$omp end parallel do

  !! Recover the full nc^3 power spectrum by filling kx = 0 to hc
  !! and merge power spectrum from threads

  Cov=0
  weight=0

  do it=1,nt
     Cov=Cov+Covt(:,it)
     weight=weight+weightt(:,it)
  enddo


  !do it=1,nt
     !do kz=1,nc-2  
        !do ky=1,nc  
           !do kx=1,nc
              !if(kx < hc+1)then
                 !Covt(kx,ky,kz,it)=Covt(nc-kx,nc-ky,nc-kz,it)
              !endif
              !Cov(kx,ky,kz)=Cov(kx,ky,kz)+Covt(kx,ky,kz,it)
   
              ! Place a maximum value to Power Spectrum...
              !if( (Cov(kx,ky,kz) > MAX) .or. (Cov(kx,ky,kz) < -MAX) ) then
              !   Cov(kx,ky,kz) = 0 ! = MAX ! ??
              !   write(*,*) 'coordinates of extremum :', kx,ky,kz 
              !   Nextreme = Nextreme + 1
              !endif


           !enddo
        !enddo
     !enddo
 
     !Cov=Cov+Covt(:,:,:,it)
     !weight=weight+weightt(:,it)
  !enddo


  !do ii = 1,nc
     !write(*,*) 'Cov(', ii, ', 3hc/4, 3hc/4) = ', Cov(ii,3*hc/4,3*hc/4)  
  !enddo


  !write(*,*) 'Test : Cov(-3*hc/4, 3*hc/4, 3*hc/4) = ', Cov(hc-3*hc/4,hc+3*hc/4,hc+3*hc/4)
  !write(*,*) 'Test : Cov(3*hc/4, -3*hc/4, -3*hc/4) = ', Cov(hc+3*hc/4,hc-3*hc/4,hc-3*hc/4)
  !write(*,*) 'Test : Cov(-hc/4, hc/4, hc/4) = ', Cov(hc-hc/4,hc+hc/4,hc+hc/4)
  !write(*,*) 'Test : Cov(hc/4, -hc/4, -hc/4) = ', Cov(hc+hc/4,hc-hc/4,hc-hc/4)

  !This test was succesfull!!!

  !write(*,*) 'N extreme = ', Nextreme
  write(*,*) 'Covariance : ', maxval(Cov), ' > Cov > ',minval(Cov)


  return
end subroutine CovAve_c2r

!*********************************




!********************************

! 1- Get Vectorial Power spectrum
! 2- Fill the symmetric missing part in negative kx due to FFT

! Getting the Vectorial Power Spectrum
! Requires : Re(d) is from kx = [hc+2,nc+2]
!            Im(d) is from kx = [1,hc+1]
! PS has dimension(nc+2,nc,nc), for later operations
! but is filled only as (nc,nc,nc)

subroutine GetPS(d,PS,nc)

  implicit none

  integer nc,hc
  integer, parameter :: nt=4


  real, dimension(nc+2,nc,nc) :: d,PS
  !real, dimension(nc,nt) :: pst,weightt
  !real, dimension(nc+2,nc,nc, nt) :: Covt
  
  integer it,i,j,k,kpt, Nextreme, Nnegative
  real pi,pow,MAX !kr,k1,k2,w1,w2

  pi=acos(-1.)
  hc=nc/2
  kpt=nc/nt
  MAX = 100000000000

  ! Get PS 
  PS=0
  do k=1,nc
     do j=1,nc
        do i=1,hc ! leaving away the 2 h+1 terms!

           ! PS = [Real(d)^2 + Im(d)^2]/V 
           PS(i+hc,j,k) = ( (d(i+hc+1,j,k))**2 + (d(i,j,k))**2 )/nc**3

           
           !PS(i,j,k) = ( (d(2*hc+1-i,j,k))**2 + (d(hc-i,nc-j,nc-k))**2 )/nc**3
  
           ! Explore negative values of Power Spectrum...
           if(PS(i+hc,j,k) < 0 ) then
              PS(i+hc,j,k) = 0 ! = MAX ! ??
              write(*,*) 'coordinates of negative values :', i,j,k 
              Nnegative = Nnegative + 1
           endif

           ! Place a maximum value to Power Spectrum...
           if( (PS(i+hc,j,k) > MAX) .or. (PS(i+hc,j,k) < -MAX) ) then
              PS(i+hc,j,k) = 0 ! = MAX ! ??
              write(*,*) 'coordinates of extremum :', i,j,k 
              Nextreme = Nextreme + 1
           endif

           ! PS is positive everywhere here

        enddo
     enddo
  enddo
 
  
  write(*,*) 'Before reassignment, N negative = ', Nnegative

 
  ! Reconstruct the missing half of PS
  do k=1,nc
     do j=1,nc
        do i=1,hc

           ! PS(-kx, ky, kz) = PS(kx, -ky, -kz) 
           PS(i,j,k) = PS(nc+1-i,nc+1-j,nc+1-k)
 
           ! Explore negative values of Power Spectrum...
           if(PS(i,j,k) < 0 ) then
              PS(i,j,k) = 0 ! = MAX ! ??
              !write(*,*) 'coordinates of negative values :', i,j,k               
              Nnegative = Nnegative + 1
           endif

               
        enddo
     enddo
  enddo

  !write(*,*) 'Test : PS(-50, 50, 50) = ', PS(hc-50,hc+50,hc+50)
  !write(*,*) 'Test : PS(50, -50, -50) = ', PS(hc+50,hc-50,hc-50)
  !write(*,*) 'Test : PS(-5, 5, 5) = ', PS(hc-5,hc+5,hc+5)
  !write(*,*) 'Test : PS(5, -5, -5) = ', PS(hc+5,hc-5,hc-5)

  !This test was succesfull!!!

  write(*,*) 'N extreme = ', Nextreme
  write(*,*) 'After reassignment, N negative = ', Nnegative
  write(*,*) 'Power spectrum : ', maxval(PS), ' > PS > ',minval(PS)

  return
end subroutine GetPS



!********************************

!! Get Vectorial power spectrum: psV(kx, ky, kz) from density d(i,j,k), 
!! where kx, ky and kz run from 1 to nc.

!! Don't forget that the real wave numbers run from -257 to 256...



subroutine VectPS_r2c(d,psV,nc)
  implicit none

  integer nc,hc
  integer, parameter :: nt=4

  !real, dimension(nc) :: ps,weight
  real, dimension(nc+2,nc,nc) :: d, psV
  !real, dimension(nc,nt) :: pst,weightt
  real, dimension(nc+2,nc,nc, nt) :: pstV
  
  integer it,i,j,k,kpt, Nextreme
  real kx,ky,kz,pi,pow,MAX !kr,k1,k2,w1,w2



  pi=acos(-1.)
  hc=nc/2
  kpt=nc/nt
  MAX = 100000000000

  write(*,*) 'd before fft : ', maxval(d),' > d > ' ,minval(d)

  call sfft3_r2c(d,nc,1) 

  write(*,*) 'd after fft : ', maxval(d),' > d > ' ,minval(d)



  ! Now d is in K space. Note that there should be only
  ! half as many k-modes as cells per dim, 
  ! So I should loop from 1 to nc/2 when writing power spectrum,

  !$omp parallel do default(shared) &
  !$omp& private(it,i,j,k,kr,kx,ky,kz,k1,k2,w1,w2,pow)
  do it=1,nt
     pstV(:,:,:,it)=0
     !weightt(:,it)=0
     do k=1+(it-1)*kpt,min(nc,it*kpt)
        if (k .lt. hc+2) then
            kz=k-1
        else
            kz=k-1-nc
        endif
        ! I will shift values of ky and kz so that they run from 1 to nc
        kz = kz + hc
        
        do j=1,nc
           if (j .lt. hc+2) then
              ky=j-1
           else
              ky=j-1-nc
           endif
           ky = ky + hc

           do i=1,nc+1,2 !(i odd = real, i even = imaginary)
              kx=(i-1)/2 
              
              ! kx runs from 0 to hc. 
              ! I will shift it so that it runs from hc to nc
 
              kx = kx + hc
                         
              !write(*,*) 'Test : kx = ', kx, 'ky = ', ky, 'kz = ',kz
              !pause

              ! power spectrum from threads: [Real(d)^2 + Im(d)^2]/V
              pstV(kx,ky,kz,it)=(d(i,j,k)**2+d(i+1,j,k)**2)/nc**3
              

              !kr=sqrt(kx**2+ky**2+kz**2)
              !if (kr .ne. 0) then
                  !k1=ceiling(kr)
                  !k2=k1+1
                  !w1=k1-kr
                  !w2=1-w1
                  !pow=(d(i,j,k)/nc**3)**2+(d(i+1,j,k)/nc**3)**2
                  !weightt(k1,it)=weightt(k1,it)+w1
                  !pst(k1,it)=pst(k1,it)+w1*pow
                  !weightt(k2,it)=weightt(k2,it)+w2
                  !pst(k2,it)=pst(k2,it)+w2*pow
              !endif
           enddo
        enddo
     enddo
  enddo
  !$omp end parallel do

  !! Recover the full nc^3 power spectrum by filling kx = 0 to hc
  !! and merge power spectrum from threads

  psV=0
  !weight=0
  do it=1,nt
     do kz=1,nc-2  
        do ky=1,nc  
           do kx=1,nc
              if(kx < hc+1)then
                 pstV(kx,ky,kz,it)=pstV(nc-kx,nc-ky,nc-kz,it)
              endif

              psV(kx,ky,kz)=psV(kx,ky,kz)+pstV(kx,ky,kz,it)

              ! Place a maximum value to Power Spectrum...
              if( (psV(kx,ky,kz) > MAX) .or. (psV(kx,ky,kz) < -MAX) ) then
                 psV(kx,ky,kz) = 0 ! = MAX ! ??
                 write(*,*) 'coordinates of extremum :', kx,ky,kz 
                 Nextreme = Nextreme + 1
              endif

           enddo
        enddo
     enddo
 
     !psV=psV+pstV(:,:,:,it)
     !weight=weight+weightt(:,it)
  enddo

  !write(*,*) 'Test : psV(-50, 50, 50) = ', psV(hc-50,hc+50,hc+50)
  !write(*,*) 'Test : psV(50, -50, -50) = ', psV(hc+50,hc-50,hc-50)
  !write(*,*) 'Test : psV(-5, 5, 5) = ', psV(hc-5,hc+5,hc+5)
  !write(*,*) 'Test : psV(5, -5, -5) = ', psV(hc+5,hc-5,hc-5)

  !This test was succesfull!!!

  write(*,*) 'N extreme = ', Nextreme
  write(*,*) 'Power spectrum : ', maxval(psV), ' > P > ',minval(psV)


  return
end subroutine VectPS_r2c
