module StringOpt
  
contains
  
  function digitn(i,n)
    implicit none
    integer, intent(in) :: i,n
    character(len=n) digitn

    integer, dimension(n) :: j
    integer k

    if(i.ge.0) then
      j(n)=i/10**(n-1)
      if(n>1) then
         do k=n-1,1,-1
            j(k)=mod(i,10**k)/10**(k-1)
         enddo
      endif
      digitn=char(ichar('0')+j(1))
      if(n>1) then
         do k=2,n
            digitn=char(ichar('0')+j(k))//adjustl(digitn)
         enddo
      endif
    else
      j(n-1)=abs(i)/10**(n-2)
      if(n-1>1) then
         do k=n-2,1,-1
            j(k)=mod(abs(i),10**k)/10**(k-1)
         enddo
      endif
      digitn=char(ichar('0')+j(1))
      if(n-1>1) then
         do k=2,n-1
            digitn=char(ichar('0')+j(k))//adjustl(digitn)
         enddo
      endif
      digitn='-'//adjustl(digitn)
    endif 
     

    return
  end function digitn

end module StringOpt

