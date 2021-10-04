
subroutine index_k2_3d(index_shell,arr_shell,nc)
  implicit none
  
  integer, intent(in) :: nc
  integer, dimension(2,3*nc**2/4), intent(out) :: index_shell
  integer, dimension(3*nc**2/4), intent(out) :: arr_shell  

  integer, dimension(3*nc**2/4) :: num_shell  
  
  integer hc,k2,i_shell,i,j,k,kx,ky,kz
   
  hc=nc/2

  num_shell=0 

!  do k=0,hc
!      do j=0,hc
!         do i=0,hc
!            k2=i**2+j**2+k**2 
!            if(k2.gt.0)num_shell(k2)=num_shell(k2)+1
!         enddo
!      enddo
!  enddo     
  
  do k=1,nc
     kz=k-1
     if(kz>hc)kz=k-1-nc
     do j=1,nc
        ky=j-1
        if(ky>hc)ky=j-1-nc
        do i=1,nc
           kx=i-1
           if(kx>hc)kx=i-1-nc
           k2=kx**2+ky**2+kz**2 
           if(k2.gt.0)num_shell(k2)=num_shell(k2)+1
        enddo
     enddo
  enddo   


  index_shell(:,:)=0
  i_shell=0
  do i=1,3*nc**2/4
     if(num_shell(i).ne.0)then 
        i_shell=i_shell+1 
        index_shell(2,i)=i_shell
        index_shell(1,i)=num_shell(i)
     endif 
  enddo
   
  arr_shell(:)=0
  do k=1,nc
     kz=k-1
     if(kz>hc)kz=k-1-nc
     do j=1,nc
        ky=j-1
        if(ky>hc)ky=j-1-nc
        do i=1,nc
           kx=i-1
           if(kx>hc)kx=i-1-nc
           k2=kx**2+ky**2+kz**2 
           if(k2.gt.0)then
              arr_shell(index_shell(2,k2))=k2 
           endif
        enddo
     enddo
  enddo   

  return
end subroutine index_k2_3d
