 
 subroutine PDF_3d(map,nc,PDF,rhobin,nbin)
   implicit none

   integer nc,nbin
   real, dimension(nc,nc,nc), intent(in) :: map
   real, dimension(nbin),intent(in) :: rhobin
   real, dimension(nbin),intent(out) :: PDF

   integer ibin,i,j,k
   
   PDF=0
if(0)then
   do k=1,nc
      do j=1,nc
         do i=1,nc
            ibin=1
            !do while(map(i,j,k).lt.rhobin(ibin).and.ibin.gt.1)   
            do while(map(i,j,k).ge.rhobin(ibin).and.ibin.lt.nbin)   
               ibin=ibin+1   
            enddo 
            PDF(ibin-1)=PDF(ibin-1)+1 
          enddo
       enddo
   enddo
endif
if(1)then
   do k=1,nc
       do j=1,nc
          do i=1,nc
             ibin=(map(i,j,k)-rhobin(1))/(rhobin(2)-rhobin(1))+1
             if(ibin.gt.nbin)ibin=nbin
             if(ibin.lt.1)ibin=1
             PDF(ibin)=PDF(ibin)+1
          enddo
       enddo
    enddo
endif

   return 
 end subroutine PDF_3d

! subroutine GaussPdf(map,nc,sigma)
!    implicit none
!
!    integer, parameter :: nbin=400
!    integer np=nc**3
!
!    integer, dimension(nbin+1) :: npbin
!    integer, dimension(np) :: rholist
!    integer, dimension(nbin) :: ncount1,ncount2
!    real, dimension(nc,nc,nc) :: map
!    real, dimension(np) :: rho_tmp,rhol
!    real, dimension(nbin+1) :: rhobin
!
!    integer i,j,k,ip
!    integer n1
!    real minrho,maxrho,sigma,rhomean,erf
!
!    equivalence(map,rhol)
!
!
!    !! Choose the parameters for the Gaussian distribution
!    !! 5 sigma deviation
!    maxrho=5*sigma
!    minrho=-5*sigma
!
!    rhol=0
!    rhobin=0
!    npbin=0
!    ncount1=0
!    ncount2=0
!    rholist=0
!    rho_tmp=0
!
!    ip=1
!    do k=1,nc
!       do j=1,nc
!          do i=1,nc
!             rhol(ip)=map(i,j,k)
!             ip=ip+1
!          enddo
!       enddo
!    enddo
!
!    do i=1,nbin
!       rhobin(i)=minrho+(maxrho-minrho)/nbin*(i-1)
!       npbin(i)=(erf((rhobin(i))/sigma/sqrt(2.))/2-erf((minrho)/sigma/sqrt(2.))/2)*np
!    enddo
!    rhobin(nbin+1)=maxrho
!    npbin(nbin+1)=np
!    do i=1,nbin
!       ncount1(i)=count(rhol.le.rhobin(i+1))-count(rhol.le.rhobin(i))
!       !!ncount2(i)=npbin(i+1)-npbin(i)
!    enddo

!    do ip=1,np
!       rholist(ip)=ip
!    enddo

!    !! SORTING
!    !! sort2 is NUMERICAl RECIPE
!    !! grade_up is from f90 (Compaq libary?)
!    !! indexsort is from Hugh Merz
!    !call sort2(np,rhol,rholist)
!    call indexedsort(np,rhol,rholist)
!    !rholist=grade_up(rhol,dim=1)
!
!    ip=0
!    do i=1,nbin-1
!       if(ip.lt.np) then
!          ip=npbin(i)+1
!          if(ip.lt.np)then
!             do while(rhol(ip+1).lt.rhobin(i+1).and.ip.lt.npbin(i+2))
!                ip=ip+1
!             enddo
!          endif
!          if(ip.lt.npbin(i+1)) then
!             rhol(ip:npbin(i+1))=(rhobin(i)+rhobin(i+1))/2
!          !elseif(i.lt.nbin) then
!          else
!             rhol(max(npbin(i+1),1):ip)=(rhobin(i+1)+rhobin(i+2))/2
!          endif
!          if(rhol(ip).lt.minrho)write(*,*) 'bug',rhol(ip),ip,minrho,i
!       endif
!    enddo
!    if(ip.lt.np) then
!       rhol(ip:np)=(rhobin(nbin)+rhobin(nbin+1))/2
!    endif
!    rho_tmp=rhol
!    do i=1,np
!       rhol(rholist(i))=rho_tmp(i)
!    enddo
!   
!     !! ncount1 means the PDF before the Gaussianization
!    !! ncount2 means the PDF after the Gaussianization
!    do i=1,nbin
!       ncount2(i)=count(rhol.le.rhobin(i+1))-count(rhol.le.rhobin(i))
!    enddo
!!    open(20,file='pdf_gauss.dat')
!!    do i=1,nbin
!!       write(20,*) (rhobin(i)+rhobin(i+1))/2,ncount2(i)
!!    enddo
!!    close(20)
!!    open(20,file='pdf_ngauss.dat')
!!    do i=1,nbin
!!       write(20,*) (rhobin(i)+rhobin(i+1))/2,ncount1(i)
!!    enddo
!!    close(20)
!
!    return
!  end subroutine GaussPdf
