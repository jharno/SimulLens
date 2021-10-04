program testdriver
  implicit none

  integer, parameter :: nc=8
  integer, dimension(3*nc**2/4) :: index_shell,num_p
 
  external index_k2_3d

  integer i,j,k,hc,k2  

  call index_k2_3d(index_shell,nc)

  hc=nc/2
  num_p=0
  do k=0, hc
     do j=0, hc
        do i=0, hc
           k2=i**2+j**2+k**2
           if(k2.gt.0)num_p(k2)=num_p(k2)+1   
        enddo
     enddo
  enddo 

  do i=1, 3*nc**2/4
     write(*,*) num_p(i), index_shell(i)
  enddo

end program testdriver
