! for a given redshift, solve for z.
! chi in Mpc/h. Stops at chi_max for now...
! uses nrecip.f. Make sure NMAX = n_points in spline subroutine
!
! Accurate in the range 1e-5 < z < 5.87 (chi_max = 5800 Mpc/h)
! Lower redshifts are assigned zero

Function InverseChi(my_chi,omegam)
  implicit none


  integer, parameter:: n_points = 55 !79 !52
  real(4) omegam, chi, my_chi, z_int(n_points),chi_int(n_points), d2ydx2(n_points),z_min, z_max, chi_max, dchi_low, dchi_high,chi_out, dz(n_points)
  real(4) InverseChi, h
  integer(4) i
  logical, parameter :: ReadTable = .false.

 !---------------------
 if(my_chi .eq.0.0)then
    InverseChi = 0.0
    return
 elseif(my_chi .lt. 0.0) then
    write(*,*) '*** ERROR, chi negative!!! ***'
 else
 
    h = 0.7 ! will be factored out in the end anyway...

    ! Do the forwards calculation of chi(z) over n_points around the specified z, in two steps.
    ! Store the arrays [chi_int z] then linearly interpolate in reverse order with spline
    ! Do coarse, then fine.

    z_int(1) = 0.00001
    z_int(n_points) =  5.894
    chi_int(1) = chi(z_int(1), omegam, h)*h
    chi_max = 5800.0
    !chi_max = chi(z_int(n_points),omegam,h)*h
 
    !---------------------
    if(my_chi .gt. chi_max)then
       write(*,*) '*** ERROR, chi too large. Must modify the InverseChi code *** '
       !write(*,*) 'my_chi = ', my_chi, 'chi_max', chi_max
       InverseChi = -1
       return
    endif
    if(my_chi .le. 4.5e-2)then
       !write(*,*) '*** ERROR, chi too small. Must modify the InverseChi code *** '
       !write(*,*) 'my_chi = ', my_chi, 'chi_min = 4.5e-2'
       InverseChi = 0
       return
    endif
    !---------------------

    if(.not. ReadTable)then
       !Produce a table of the form:
       !write(*,*) 'Chi      z'

       !get linear bining in z:
       !open(11, file='ChiZ_TableWMAP9.dat', status = 'new')
       !write(11,*) chi_int(1), z_int(1)
       !do i = 2,n_points
       !   dz(i) = z_int(n_points)/real(n_points)
       !   z_int(i) = z_int(i-1) + dz(i)
       !enddo

       !get log bining in z, tuned for 52 bins:
       z_int =(/1e-5, 1.3e-5,1.6e-5, 2.2e-5, 2.8e-5, 3.7e-5, 4.8e-5, 6.2e-5, 8.0e-5, 1.0e-4, 1.3e-4, 1.7e-4, 2.2e-4, 3.0e-4, 3.8e-4, 5.0e-4,6.4e-4, 8.4e-4, 1.1e-3, 1.4e-3, 1.8e-3, &
            2.4e-3, 3.1e-3, 4.0e-3, 5.2e-3, 6.7e-3, 8.7e-3, 1.1e-2, 1.5e-2, 1.9e-2, 2.4e-2, 3.2e-2, 4.1e-2, 5.4e-2, 7.0e-2, 9.1e-2, 1.2e-1, 1.5e-1, 2.0e-1, 2.6e-1, 3.4e-1, 4.4e-1, &
             5.7e-1, 7.3e-1, 9.5e-1, 1.2, 1.6, 2.1, 2.7, 3.5, 4.0, 4.542, 5.0, 5.5, 5.894/)

       !get log bining + a few more points in z, tuned for 52 bins:
       !z_int =(/1.0e-5, 1.3e-5, 1.6e-5, 2.2e-5, 2.8e-5, 3.7e-5, 4.8e-5, 6.2e-5, 8.0e-5, 1.0e-4, 1.3e-4, 1.7e-4, 2.2e-4, 3.0e-4, 3.8e-4, 5.0e-4, 6.4e-4, 8.4e-4, 1.1e-3, 1.4e-3, 1.8e-3, &
       !         2.4e-3, 3.1e-3, 4.0e-3, 5.2e-3, 6.7e-3, 8.7e-3, 1.1e-2, 1.5e-2, 1.9e-2, 2.4e-2, 3.2e-2, 4.1e-2, 5.4e-2, 7.0e-2, 9.1e-2, 1.0e-1, 1.1e-1, 1.2e-1, 1.4e-1, 1.5e-1, 1.8e-1, &
       !         2.0e-1, 2.3e-1, 2.6e-1, 3.0e-1, 3.4e-1, 3.8e-1, 4.4e-1, 5.7e-1, 6.1e-1, 6.5e-1, 6.9e-1, 7.3e-1, 7.8e-1, 8.2e-1, 8.6e-1, 9.0e-1, 9.5e-1, 9.9e-1, &
       !         1.03,   1.05,   1.08,   1.14,   1.18,   1.21,   1.27,   1.34,   1.40,   1.50,   1.60,   1.80,   2.10,   2.40,   2.70,   3.10,   3.50,   4.54,   5.894/)

       !open(11, file='ChiZ_TableWMAP9_log_50.dat', status = 'new')
       do i = 2,n_points
          chi_int(i) = chi(z_int(i),omegam,h)*h 
       !   write(11,*) chi_int(i), z_int(i)
       enddo
       !close(11)
       !stop
       !----------
    else
       open(11, file='ChiZ_TableWMAP9.dat', status = 'old')
       do i = 1,n_points
          read(11,*) chi_int(i), z_int(i)
          !write(*,*) i,chi_int(i), z_int(i)
       enddo
       close(11)
       
    endif

    !----------------------------------------------
    !3- interpolate z onto  chi  
    !dchi_low  = (chi_int(2) - chi_int(1))/(z_int(2) - z_int(1))
    !dchi_high = (chi_int(n_points) - chi_int(n_points-1))/(z_int(n_points) - z_int(n_points-1))
     
    !call spline(z_int, chi_int, n_points, dchi_low,dchi_high,d2ydx2)
    !call splint(z_int, chi_int, d2ydx2, n_points, 0.042, chi_out)

    !----------------------------------------------
    !3- interpolate Chi onto  Z 
    dchi_low  = 1e-12!(z_int(2) - z_int(1))/(chi_int(2) - chi_int(1))
    dchi_high = 1e-12! (z_int(n_points) - z_int(n_points-1))/(chi_int(n_points) - chi_int(n_points-1))     
     
    call spline(chi_int, z_int, n_points, dchi_low,dchi_high,d2ydx2)
    call splint(chi_int, z_int, d2ydx2, n_points, my_chi, InverseChi)

    !-----------------
    ! Do again, but on a finer grid:

    !if(InverseChi .lt. 1.5e-5 ) 


    z_int(1) = InverseChi*0.8
    z_int(n_points) =  InverseChi*1.2
    
    !get linear bining in z:
    !write(11,*) chi_int(1), z_int(1)
    do i = 2,n_points
       dz(i) = (z_int(n_points)-z_int(1))/real(n_points)
       z_int(i) = z_int(i-1) + dz(i)
    enddo

    do i = 1,n_points
       chi_int(i) = chi(z_int(i),omegam,h)*h 
       !   write(11,*) chi_int(i), z_int(i)
    enddo

    call spline(chi_int, z_int, n_points, dchi_low,dchi_high,d2ydx2)
    call splint(chi_int, z_int, d2ydx2, n_points, my_chi, InverseChi)


    !InverseChi = 0.5

 endif
 return
end function InverseChi
