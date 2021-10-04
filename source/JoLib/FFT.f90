  
  subroutine sfft1_r2c(d,n,c)
    implicit none
    include 'fftw_f77.i'

    integer n
    integer c
    real d(n+2)
    integer(8) :: plan,iplan

    if (c .eq. 1) then
       call sfftw_plan_dft_r2c_1d(plan,n,d,d,FFTW_ESTIMATE)
       call sfftw_execute(plan)
       call sfftw_destroy_plan(plan)
    endif
    if (c .eq. -1) then
       call sfftw_plan_dft_c2r_1d(iplan,n,d,d,FFTW_ESTIMATE)
       call sfftw_execute(iplan)
       d = d / real(n)
       call sfftw_destroy_plan(iplan)
    endif
    if (c .ne. 1 .and. c .ne. -1) write(*,*) 'FFT direction error'

    return
  end subroutine sfft1_r2c

  subroutine sfft1_c2c(d,n,c)
    implicit none
    include 'fftw_f77.i'

    integer n
    integer c
    complex d(n)
    integer(8) :: plan,iplan

    if (c .eq. 1) then
       call sfftw_plan_dft_1d(plan,n,d,d,FFTW_FORWARD,FFTW_ESTIMATE)
       call sfftw_execute(plan)
       call sfftw_destroy_plan(plan)
    endif
    if (c .eq. -1) then
       call sfftw_plan_dft_1d(iplan,n,d,d,FFTW_BACKWARD,FFTW_ESTIMATE)
       call sfftw_execute(iplan)
       d = d / real(n)
       call sfftw_destroy_plan(iplan)
    endif
    if (c .ne. 1 .and. c .ne. -1) write(*,*) 'FFT direction error'

    return
  end subroutine sfft1_c2c

  subroutine sfft2_r2c(d,n,c)
    implicit none
    !include 'fftw_f77.i'
    include 'fftw3.f'

    integer n,nt
    integer c
    real d(n+2,n)
    integer(8) :: plan,iplan

    nt = 24

#ifdef OMP
    call sfftw_init_threads()
    call sfftw_plan_with_nthreads(nt)
    write(*,*)'Performing FFT with n_thread =',nt
#endif

    if (c .eq. 1) then
       call sfftw_plan_dft_r2c_2d(plan,n,n,d,d,FFTW_ESTIMATE)
       call sfftw_execute(plan)
       call sfftw_destroy_plan(plan)
    endif
    if (c .eq. -1) then
       call sfftw_plan_dft_c2r_2d(iplan,n,n,d,d,FFTW_ESTIMATE)
       call sfftw_execute(iplan)
       d = d / real(n)**2
       call sfftw_destroy_plan(iplan)
    endif
    if (c .ne. 1 .and. c .ne. -1) write(*,*) 'FFT direction error'
#ifdef OMP
   call sfftw_cleanup_threads()
#endif

    return
  end subroutine sfft2_r2c

  subroutine sfft2_c2c(d,n,c)
    implicit none
!    include 'fftw_f77.i'
    include 'fftw3.f'

    integer n,nt
    integer c
    complex d(n,n)
    integer(8) :: plan,iplan
    !nt = 24
    nt = 8

#ifdef OMP
    call sfftw_init_threads()
    call sfftw_plan_with_nthreads(nt)
    write(*,*)'Performing FFT with n_thread =',nt
#endif
    if (c .eq. 1) then
       call sfftw_plan_dft_2d(plan,n,n,d,d,FFTW_FORWARD,FFTW_ESTIMATE)
       call sfftw_execute(plan)
       call sfftw_destroy_plan(plan)
    endif
    if (c .eq. -1) then
       call sfftw_plan_dft_2d(iplan,n,n,d,d,FFTW_BACKWARD,FFTW_ESTIMATE)
       call sfftw_execute(iplan)
       d = d / real(n)**2
       call sfftw_destroy_plan(iplan)
    endif
    if (c .ne. 1 .and. c .ne. -1) write(*,*) 'FFT direction error'

#ifdef OMP
   call sfftw_cleanup_threads()
#endif
    return
  end subroutine sfft2_c2c

!------------------------

  subroutine sfft3_r2c(d,n,c)
    implicit none
    include 'fftw_f77.i'

    integer n,nt
    integer c
    real d(n+2,n,n)
    integer(8) :: plan,iplan
    nt=4

#ifdef OMP
    call sfftw_init_threads()
    call sfftw_plan_with_nthreads(nt)
    write(*,*)'using OMP in FFT'
#endif
    if (c .eq. 1) then
       call sfftw_plan_dft_r2c_3d(plan,n,n,n,d,d,FFTW_ESTIMATE)
       call sfftw_execute(plan)
       call sfftw_destroy_plan(plan)
    endif
    if (c .eq. -1) then
       call sfftw_plan_dft_c2r_3d(iplan,n,n,n,d,d,FFTW_ESTIMATE)
       call sfftw_execute(iplan)
       d = d / real(n)**3
       call sfftw_destroy_plan(iplan)
    endif
    if (c .ne. 1 .and. c .ne. -1) write(*,*) 'FFT direction error'
#ifdef OMP
   call sfftw_cleanup_threads()
#endif

    return
  end subroutine sfft3_r2c

  subroutine sfft3_c2c(d,n,c)
    implicit none
    include 'fftw_f77.i'

    integer n,nt
    integer c
    complex d(n,n,n)
    integer(8) :: plan,iplan
    nt=4

#ifdef OMP
    call sfftw_init_threads()
    call sfftw_plan_with_nthreads(nt)
#endif
    if (c .eq. 1) then
       call sfftw_plan_dft_3d(plan,n,n,n,d,d,FFTW_FORWARD,FFTW_ESTIMATE)
       call sfftw_execute(plan)
       call sfftw_destroy_plan(plan)
    endif
    if (c .eq. -1) then
       call sfftw_plan_dft_3d(iplan,n,n,n,d,d,FFTW_BACKWARD,FFTW_ESTIMATE)
       call sfftw_execute(iplan)
       d = d / real(n)**3
       call sfftw_destroy_plan(iplan)
    endif
    if (c .ne. 1 .and. c .ne. -1) write(*,*) 'FFT direction error'
#ifdef OMP
   call sfftw_cleanup_threads()
#endif

    return
  end subroutine sfft3_c2c
  
