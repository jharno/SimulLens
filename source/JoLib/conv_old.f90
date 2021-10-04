 subroutine npconv_2d(a,b,ab,n)
! b is the window function, b(m)=w(|m-n+1|)
    implicit none

    integer n,n2
    real, dimension(0:2*n-1,0:2*n-1) :: aa,b
    real,dimension(0:n-1,0:n-1) :: a,ab

    n2=2*n
    aa=0
    aa(0:n-1,0:n-1)=a
    call pconv_2d(aa,b,aa,2*n)
    ab=aa(0:n-1,0:n-1)

    return
  end subroutine npconv_2d

subroutine pconv_2d(a,b,ab,n)
!c=a convolve with b, period condition
    implicit none
#ifndef CXML
#ifndef LOB6
        include '/opt/fftw-3.0.1/include/fftw3.f'
#endif
#ifdef LOB6 
        include '/opt/fftw-3.0.1_gcc/include/fftw3.f'
#endif
#endif
    integer n,status,zfft_2d
    complex, dimension(0:n/2,0:n-1) :: ca,cb,cab
    real, dimension(0:n-1,0:n-1) :: a,b,ab
    integer(8) :: plan,iplan
#ifdef CXML
    external sfft_2d
    call sfft_2d('R','C','F',a,ca,n,n,n+2,1,1)
    call sfft_2d('R','C','F',b,cb,n,n,n+2,1,1)
    cab=ca*cb
    call sfft_2d('C','R','B',cab,ab,n,n,n+2,1,1)
#else
    call sfftw_plan_dft_r2c_2d(plan,n,n,a,ca,FFTW_ESTIMATE)
    call sfftw_execute(plan)
    call sfftw_destroy_plan(plan)
    call sfftw_plan_dft_r2c_2d(plan,n,n,b,cb,FFTW_ESTIMATE)
    call sfftw_execute(plan)
    call sfftw_destroy_plan(plan)
    ca=ca*cb

    call sfftw_plan_dft_c2r_2d(iplan,n,n,ca,ab,FFTW_ESTIMATE)
    call sfftw_execute(iplan)
    ab = ab / real(n)**2
    call sfftw_destroy_plan(iplan)
#endif

    return
  end subroutine pconv_2d

subroutine fft_2d(a,n,direction)
! direction='f' or 'b'
#ifndef CXML
#ifndef LOB6
        include '/opt/fftw-3.0.1/include/fftw3.f'
#endif
#ifdef LOB6 
        include '/opt/fftw-3.0.1_gcc/include/fftw3.f'
#endif
#endif
    integer n,status,cfft_2d
    character*1 direction
    !complex, dimension(0:n/2,0:n-1) :: fa
    real, dimension(0:n+1,0:n-1) :: a,fa
    integer(8) :: plan,iplan
#ifdef CXML
    external sfft_2d
    if (direction .eq. 'f') then
       call sfft_2d('R','C','F',a,a,n,n,n+2,1,1)
    else
       call sfft_2d('C','R','B',a,a,n,n,n+2,1,1)
    endif
#else
    if(direction.eq.'f')then
      call sfftw_plan_dft_r2c_2d(plan,n,n,a,a,FFTW_ESTIMATE)
      call sfftw_execute(plan)
      call sfftw_destroy_plan(plan)
    else
      call sfftw_plan_dft_c2r_2d(iplan,n,n,a,a,FFTW_ESTIMATE)
      call sfftw_execute(iplan)
      a = a / real(n)**2
      call sfftw_destroy_plan(iplan)
    endif
#endif

    return
  end subroutine fft_2d

 subroutine npconv_1d(a,b,ab,n)
! b is the window function, b(m)=w(|m-n+1|)
    implicit none

    integer n,n2
    real, dimension(0:2*n-1) :: aa,bb,cc,b
    real,dimension(0:n-1) :: a,ab

    n2=2*n
    aa=0
    bb=0
    aa(0:n-1)=a
    bb=b
    call pconv_1d(aa,bb,cc,2*n)
    ab=cc(0:n-1)

    return
  end subroutine npconv_1d

 subroutine pconv_1d(a,b,ab,n)
!c=a convolve with b, period condition
    implicit none
#ifndef CXML
#ifndef LOB6
        include '/opt/fftw-3.0.1/include/fftw3.f'
#endif
#ifdef LOB6 
        include '/opt/fftw-3.0.1_gcc/include/fftw3.f'
#endif
#endif
    integer n,status
    complex, dimension(0:n/2) :: ca,cb,cab
    !real, dimension(0:n+1) :: fa,fb
    real, dimension(0:n-1) :: a,b,ab
    integer(8) :: plan,iplan
#ifdef CXML
    external sfft
    call sfft('R','C','F',a,ca,n,1)
    call sfft('R','C','F',b,cb,n,1)
    cab=ca*cb
    call sfft('C','R','B',cab,ab,n,1)
#else
    call sfftw_plan_dft_r2c_1d(plan,n,a,ca,FFTW_ESTIMATE)
    call sfftw_execute(plan)
    call sfftw_destroy_plan(plan)
    call sfftw_plan_dft_r2c_1d(plan,n,b,cb,FFTW_ESTIMATE)
    call sfftw_execute(plan)
    call sfftw_destroy_plan(plan)
    ca=ca*cb
    call sfftw_plan_dft_c2r_1d(iplan,n,ca,ab,FFTW_ESTIMATE)
    call sfftw_execute(iplan)
    ab = ab / real(n)**1
    call sfftw_destroy_plan(iplan)
#endif

    return
  end subroutine pconv_1d

 subroutine fft_1d(a,n,direction)
! direction='f' or 'b'
    implicit none
#ifndef CXML
#ifndef LOB6
        include '/opt/fftw-3.0.1/include/fftw3.f'
#endif
#ifdef  LOB6
        include '/opt/fftw-3.0.1_gcc/include/fftw3.f'
#endif
#endif
    integer n,status,cfft
    character*1 direction
    !complex, dimension(0:n/2) :: fa
    real, dimension(0:n+1) :: a,fa
    integer(8) :: plan,iplan
#ifdef CXML
    external sfft
    if (direction .eq. 'f') then
       call sfft('R','C','F',a,a,n,1)
    else
       call sfft('C','R','B',a,a,n,1)
    endif
#else
    if(direction.eq.'f')then
      call sfftw_plan_dft_r2c_1d(plan,n,a,a,FFTW_ESTIMATE)
      call sfftw_execute(plan)
      call sfftw_destroy_plan(plan)
    else
      call sfftw_plan_dft_c2r_1d(iplan,n,a,a,FFTW_ESTIMATE)
      call sfftw_execute(iplan)
      a = a / real(n)**1
      call sfftw_destroy_plan(iplan)
    endif
#endif

    return
  end subroutine fft_1d

subroutine convolution(t,w,c,n,k)
   implicit none
   integer n,k
   real t(n,n),w(-n:n,-n:n),c(n,n)
   integer i,j,ii,jj
        do i=1,n
         do j=1,n
          c(i,j)=-w(0,0)*t(i,j)
           do ii=-k,k
           do jj=-k,k
            if(((i-ii).gt.0).and.((j-jj).gt.0)&
             &.and.((i-ii).le.n).and.((j-jj).le.n))then
              c(i,j)=c(i,j)+w(ii,jj)*t(i-ii,j-jj)
            endif
           enddo
          enddo
         enddo
        enddo
    end subroutine convolution

   subroutine win(a1,a2,a3,s1,s2,s3,b,w,n)
   implicit none
   integer n
   real w(-n:n,-n:n)
   real a1,a2,a3,s1,s2,s3,b,r
   integer i,j
        do i=-n,n
          do j=-n,n
           r=sqrt((i*b)**2+(j*b)**2)
           w(i,j)=-a1*exp(-r**2*s1**2/2)*s1**4*(-2+r**2*s1**2)&
                           &-a2*exp(-r**2*s2**2/2)*s2**4*(-2+r**2*s2**2)&
                           &-a3*exp(-r**2*s3**2/2)*s3**4*(-2+r**2*s3**2)
          enddo
        enddo
    end subroutine win

   subroutine fft_3d(a,n,c)
   implicit none
#ifndef CXML
     include '/opt/fftw-3.0.1_gcc/include/fftw3.f'
#endif
    integer n
    character*1 c
    real a(n+2,n,n)
    integer(8) :: plan,iplan
#ifdef CXML
    external sfft_3d
    if (c .eq. 'f') then
       call sfft_3d('R','C','F',a,a,n,n,n,n+2,n,1,1,1)
    else
       call sfft_3d('C','R','B',a,a,n,n,n,n+2,n,1,1,1)
    endif
#else
    if (c .eq. 'f') then
       call sfftw_plan_dft_r2c_3d(plan,n,n,n,a,a,FFTW_ESTIMATE)
       call sfftw_execute(plan)
       call sfftw_destroy_plan(plan)
     endif
    if (c .eq. 'b') then
       call sfftw_plan_dft_c2r_3d(iplan,n,n,n,a,a,FFTW_ESTIMATE)
       call sfftw_execute(iplan)
       a = a / real(n)**3
   call sfftw_destroy_plan(iplan)
     endif
     if (c .ne. 'f' .and. c .ne. 'b') write(*,*) 'FFT direction error'
#endif
    return
  end subroutine fft_3d
 
