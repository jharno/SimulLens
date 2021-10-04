
!======================================================
 Subroutine Get_KAverage(map,Ave,n,normp,normk)
!======================================================
! Calculates the Power Spectrum and Dumps it into a file


   Integer                         :: n
   Real*4                          :: normp,normk
   Real*4, Dimension(1:n+2,1:n,1:n)  :: map
   Character(Len=120)              :: file

   !Real*4, Dimension(1:n+2,1:n,1:n):: fft
   Real*4, Dimension(3,1:n)        :: pst
   Real*4, Dimension(2,1:n)        :: ps

   Real*4  :: w1,w2,kz,kx,ky,kr,pow
   Integer :: i,j,k,hn,k1,k2

   If (mod(n,2) .NE. 0) Then
      Write(*,*) ' n must be even'
      Stop
   End If


   !! == Dump power spectra
   !Open(unit=40,file=file,status='replace')
   hn  = n/2
   pst = 0.0
   Do k = 1,n
      If (k .Lt. hn+2) Then
         kz = k-1
      Else
         kz = k-1-n
      Endif
      Do j = 1,n
         If (j .Lt. hn+2) Then
            ky = j-1
         Else
            ky = j-1-n
         Endif
         Do i = 1,n+2,2
            kx = (i-1)/2
            kr = Sqrt(kx**2+ky**2+kz**2)

            If (kr .Ne. 0.) Then
               k1  = Ceiling(kr)
               k2  = k1+1
               w1  = k1-kr
               w2  = 1-w1
               pow = sqrt(Sum((map(i:i+1,j,k)/n**3)**2)) !Sum((fft(i:i+1,j,k)/n**3)**2)
               pst(1,k1)=pst(1,k1)+w1*pow
               pst(2,k1)=pst(2,k1)+w1*pow**2
               pst(3,k1)=pst(3,k1)+w1
               pst(1,k2)=pst(1,k2)+w2*pow
               pst(2,k2)=pst(2,k2)+w2*pow**2
               pst(3,k2)=pst(3,k2)+w2
            Endif
         Enddo
      Enddo
   End Do

   !! Divide by weights
   !! ps(1,k) stores p(k)
   !! ps(2,k) stores standard deviation
   Do k = 1,n
      If (pst(3,k) .Eq. 0) Then
         ps(:,k) = 0
      Else 
         ps(:,k) = pst(1:2,k)/pst(3,k)
         !ps(2,k) = Sqrt(Abs(ps(2,k)-ps(1,k)**2))
         !Write(40,'(3E18.5)') k*normk,ps(1,k)*normp,ps(2,k)*normp!4*3.141596254*(k-1)**3*ps(1:2,k)*normp
      Endif
   Enddo
   !Close(40)
   !Write(*,*) '3D power spectra written in ',Trim(file)
   Write(*,*) '3D Average Done '

   return

End Subroutine Get_KAverage

!==================================
!======================================================
 Subroutine Write_3D_powerspectra(map,n,file,normp,normk)
!======================================================
! Calculates the Power Spectrum and Dumps it into a file


   Integer                         :: n
   Real*4                          :: normp,normk
   Real*4, Dimension(1:n,1:n,1:n)  :: map
   Character(Len=120)              :: file

   Real*4, Dimension(1:n+2,1:n,1:n):: fft
   Real*4, Dimension(3,1:n)        :: pst
   Real*4, Dimension(2,1:n)        :: ps

   Real*4  :: w1,w2,kz,kx,ky,kr,pow
   Integer :: i,j,k,hn,k1,k2

     If (mod(n,2) .NE. 0) Then
      Write(*,*) 'Dump 2D powerspecta n must be even'
      Stop
   End If

   !! == Compute FFT
   fft(1:n,1:n,1:n) = map(1:n,1:n,1:n)
   Call sfftw_3d_even(fft,n,+1)

   !! == Dump power spectra
   Open(unit=40,file=file,status='replace')
   hn  = n/2
   pst = 0.0
   Do k = 1,n
      If (k .Lt. hn+2) Then
         kz = k-1
      Else
         kz = k-1-n
      Endif
      Do j = 1,n
         If (j .Lt. hn+2) Then
            ky = j-1
         Else
            ky = j-1-n
         Endif
         Do i = 1,n+2,2
            kx = (i-1)/2
            kr = Sqrt(kx**2+ky**2+kz**2)

            If (kr .Ne. 0.) Then
               k1  = Ceiling(kr)
               k2  = k1+1
               w1  = k1-kr
               w2  = 1-w1
               pow = Sum((fft(i:i+1,j,k)/n**3)**2)
               pst(1,k1)=pst(1,k1)+w1*pow
               pst(2,k1)=pst(2,k1)+w1*pow**2
               pst(3,k1)=pst(3,k1)+w1
               pst(1,k2)=pst(1,k2)+w2*pow
               pst(2,k2)=pst(2,k2)+w2*pow**2
               pst(3,k2)=pst(3,k2)+w2
            Endif
         Enddo
      Enddo
   End Do

   !! Divide by weights
   !! ps(1,k) stores p(k)
   !! ps(2,k) stores standard deviation
   Do k = 1,n
      If (pst(3,k) .Eq. 0) Then
         ps(:,k) = 0
      Else          
         ps(:,k) = pst(1:2,k)/pst(3,k)
         ps(2,k) = Sqrt(Abs(ps(2,k)-ps(1,k)**2))
         Write(40,'(3E18.5)') k*normk,ps(1,k)*normp,ps(2,k)*normp!4*3.141596254*(k-1)**3*ps(1:2,k)*normp
      Endif
   Enddo
   Close(40)
   Write(*,*) '3D power spectra written in ',Trim(file)
 
End Subroutine Write_3D_powerspectra

!==================================

!=======================================
 Subroutine  sfftw_3d_even(a,n,direction)
 !=======================================

   ! a  : real simple [1:n+2,1:n] array
   ! n  : dimension, even number
   ! Direction= +1 for forward, -1 for backward

   Include 'fftw3.f'
   Integer                       , Intent(In)    :: n,direction
   Real, Dimension(1:n+2,1:n,1:n), Intent(Inout) :: a

   Integer(8) :: plan,iplan

   If (mod(n,2) .NE. 0) Then
      Write(*,*) 'Dimension must be even in sfftw_3d_even...'
      Stop
   End If

   Select Case(direction)
   Case(+1) ! Forward
      Call sfftw_plan_dft_r2c_3d(plan,n,n,n,a,a,FFTW_ESTIMATE)
      Call sfftw_execute(plan)
      Call sfftw_destroy_plan(plan)
   Case(-1) ! Backward
      Call sfftw_plan_dft_c2r_3d(iplan,n,n,n,a,a,FFTW_ESTIMATE)
      Call sfftw_execute(iplan)
      a = a / Real(n)**3 ! Scale the FFT here
      Call sfftw_destroy_plan(iplan)
   Case Default
      Write(*,*) 'Wrong call to FFT_3D'
      Stop
   End Select
   Return

 End Subroutine sfftw_3d_even
 !===========================

