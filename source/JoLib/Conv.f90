
  subroutine pconv2d_r2c(a,b,ab,n)
    implicit none
    include 'fftw3.f'

    integer n
    complex, dimension(0:n/2,0:n-1) :: ca,cb
    real, dimension(0:n-1,0:n-1) :: a,b,ab
    integer(8) :: plan,iplan

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

    return
  end subroutine pconv2d_r2c

  subroutine npconv2d_r2c(a,b,ab,n)
    implicit none

    integer n,n2
    real, dimension(0:2*n-1,0:2*n-1) :: aa,b
    real, dimension(0:n-1,0:n-1) :: a,ab

    n2=2*n
    aa=0
    aa(0:n-1,0:n-1)=a
    call pconv2d_r2c(aa,b,aa,2*n)
    !ab=aa(n-1:2*n-2,n-1:2*n-2)
    ab=aa(0:n-1,0:n-1)

    return
  end subroutine npconv2d_r2c

  subroutine xkspace_pconv2d_r2c(a,cb,ab,n)
    implicit none

    include 'fftw3.f'
    integer n,nt 
    complex, dimension(0:n/2,0:n-1) :: ca,cb
    real, dimension(0:n-1,0:n-1) :: a,ab
    integer(8) :: plan,iplan
    nt=4

#ifdef OMP
    call sfftw_init_threads()
    call sfftw_plan_with_nthreads(nt)
#endif
    call sfftw_plan_dft_r2c_2d(plan,n,n,a,ca,FFTW_ESTIMATE)
    call sfftw_execute(plan)
    call sfftw_destroy_plan(plan)
    ca=ca*cb
    call sfftw_plan_dft_c2r_2d(iplan,n,n,ca,ab,FFTW_ESTIMATE)
    call sfftw_execute(iplan)
    ab = ab / real(n)**2
    call sfftw_destroy_plan(iplan)
#ifdef OMP
    call sfftw_cleanup_threads()
#endif
    return
  end subroutine xkspace_pconv2d_r2c

  subroutine xkspace_pconv2d_c2c(a,cb,ab,n)
    implicit none

    include 'fftw3.f'
    integer n,nt 
    complex, dimension(0:n-1,0:n-1) :: ca,cb
    real, dimension(0:n-1,0:n-1) :: a,ab
    integer(8) :: plan,iplan
    integer i,j,kx,ky
!    real rkx,rky,sig,pi,box,map_ana(n,n) 
    
!    sig=1
!    pi=acos(-1.)
!    box=50
!   do j=1,n
!       ky=j-1
!       if(ky>n/2)ky=j-1-n
!       do i=1,n
!         kx=i-1
!          if(kx>n/2)kx=i-1-n
!          rkx=kx*2*pi/box
!          rky=ky*2*pi/box
!          !map_ana(i,j)=exp(-(rkx**2+rky**2)/sig**2)*(rkx**2+rky**2)
!       enddo
!    enddo
!    map_ana=map_ana*(n/box)**2


    nt=4
#ifdef OMP
    call sfftw_init_threads()
    call sfftw_plan_with_nthreads(nt)
#endif
    ca=a+(0.,0.) 
    call sfftw_plan_dft_2d(plan,n,n,ca,ca,FFTW_FORWARD,FFTW_ESTIMATE)
    call sfftw_execute(plan)
    call sfftw_destroy_plan(plan)
!    open(10,file='kmap_gamma1_num.dat',form='binary')
!    write(10)real(ca)
!    close(10)
!    open(20,file='kmap_gamma1_ana.dat',form='binary')
!    write(20)map_ana
!    close(20)
!    open(30,file='kmap_gamma1_diff.dat',form='binary')
!    write(30)map_ana-real(ca)
!    close(30)
!    stop 
    ca=ca*cb
    call sfftw_plan_dft_2d(iplan,n,n,ca,ca,FFTW_BACKWARD,FFTW_ESTIMATE)
    call sfftw_execute(iplan)
    ab = ca / real(n)**2
    call sfftw_destroy_plan(iplan)
#ifdef OMP
    call sfftw_cleanup_threads()
#endif
    return
  end subroutine xkspace_pconv2d_c2c


  subroutine xkspace_pconv3d_c2c(a,cb,ab,n)
    implicit none

    include 'fftw3.f'
    integer n,nt 
    complex, dimension(0:n-1,0:n-1,0:n-1) :: ca,cb
    real, dimension(0:n-1,0:n-1,0:n-1) :: a,ab
    integer(8) :: plan,iplan
    integer i,j,kx,ky
    
    nt=4
#ifdef OMP
    call sfftw_init_threads()
    call sfftw_plan_with_nthreads(nt)
#endif
    ca=a+(0.,0.) 
    call sfftw_plan_dft_3d(plan,n,n,n,ca,ca,FFTW_FORWARD,FFTW_ESTIMATE)
    call sfftw_execute(plan)
    call sfftw_destroy_plan(plan)
    !write(*,*) ca
    !write(*,*) (ca(1,0,0)/n**3)*(conjg(ca(1,0,0))/n**3) 
    !pause
    !write(*,*) cb  
    !pause 
    ca=ca*cb
    !write(*,*) sum(ca/n**3*conjg(ca)/n**3) 
      
    !ca=cb
    call sfftw_plan_dft_3d(iplan,n,n,n,ca,ca,FFTW_BACKWARD,FFTW_ESTIMATE)
    call sfftw_execute(iplan)
    ab = ca / real(n)**3
    call sfftw_destroy_plan(iplan)

    !write(*,*) sum((ca/n**3)*(conjg(ca)/n**3))/n**3  
    !write(*,*) sum(ab*ab)/n**3  
    !stop
    !open(100,file='windowx.dat')
    !do k=1,n
    ! do j=1,n
    !    do i=1,n
    !       k1=k-1
    !       if(k1>n/2)k1=k-1-n
    !       j1=j-1
    !       if(j1>n/2)j1=j-1-n
    !       i1=i-1
    !       if(i1>n/2)i1=i-1-n
    !       write(100,*) sqrt(i1**2+j1**2)*50/n,k1*50/n,cb 
    !    enddo
    ! enddo
    !enddo
 
#ifdef OMP
    call sfftw_cleanup_threads()
#endif
    return
  end subroutine xkspace_pconv3d_c2c


  subroutine xkspace_pconv3d_r2c(a,cb,ab,n)
    implicit none

    include 'fftw3.f'

    integer n,nt
    complex, dimension(0:n/2,0:n-1,0:n-1) :: ca,cb
    real, dimension(0:n-1,0:n-1,0:n-1) :: a,ab

    integer(8) :: plan,iplan
    nt=4

#ifdef OMP
    call sfftw_init_threads()
    call sfftw_plan_with_nthreads(nt)
#endif
    call sfftw_plan_dft_r2c_3d(plan,n,n,n,a,ca,FFTW_ESTIMATE)
    call sfftw_execute(plan)
    call sfftw_destroy_plan(plan)

    ca=ca*cb
    ab=0.
    call sfftw_plan_dft_c2r_3d(iplan,n,n,n,ca,ab,FFTW_ESTIMATE)
    call sfftw_execute(iplan)
    ab = ab / real(n)**3
    call sfftw_destroy_plan(iplan)

#ifdef OMP
   call sfftw_cleanup_threads()
#endif

    return
  end subroutine xkspace_pconv3d_r2c

