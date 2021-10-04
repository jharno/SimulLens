
!   w1=0
!   do j=1,nc
!      do i=1,nc
!         r=sqrt((i-1.)**2+(j-1.)**2)*gbox/npc
!         w1(i,j)=exp(-r**2/2*sig2)*sig2  
!      enddo
!   enddo
!   w1=w1/2/pi*(box/nc)**2

subroutine PPadingWindow(w1,w2,n)
  implicit none
!!periodic field, padding

  integer n,n2,i,j
  real, dimension(n,n) :: w1,w2

  w2=0
  n2=n
  do j=1,n2
      do i=1,n2
         if((i<=n2/2).and.(j<=n2/2))then
            w2(i,j)=w1((i-1)+1,(j-1)+1)
         elseif((i>n2/2+1).and.(j<=n2/2))then
            w2(i,j)=w1(abs(i-n2-1)+1,(j-1)+1)
         elseif((i<=n2/2).and.(j>n2/2+1))then
            w2(i,j)=w1((i-1)+1,abs(j-n2-1)+1)
         elseif((i>n2/2+1).and.(j>n2/2+1))then
            w2(i,j)=w1(abs(i-n2-1)+1,abs(j-n2-1)+1)
         endif
      enddo
   enddo
  
  return
end subroutine PPadingWindow


subroutine NpPadingWindow(w1,w2,n)
  implicit none
!!periodic field, padding
  
  integer n,n2,i,j
  real, dimension(n,n) :: w1
  real, dimension(2*n,2*n) :: w2   

  n2=2*n
  w2=0
  do j=1,n2
      do i=1,n2
         if((i<=n2/2).and.(j<=n2/2))then
            w2(i,j)=w1((i-1)+1,(j-1)+1)
         elseif((i>n2/2+1).and.(j<=n2/2))then
            w2(i,j)=w1(abs(i-n2-1)+1,(j-1)+1)
         elseif((i<=n2/2).and.(j>n2/2+1))then
            w2(i,j)=w1((i-1)+1,abs(j-n2-1)+1)
         elseif((i>n2/2+1).and.(j>n2/2+1))then
            w2(i,j)=w1(abs(i-n2-1)+1,abs(j-n2-1)+1)
         endif
      enddo
   enddo
  
  return
end subroutine NpPadingWindow
