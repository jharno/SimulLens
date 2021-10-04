program build_slabs
  implicit none
include 'project.fh'

!1- Read projections from all nodes

!integer, parameter, :: buffer
!real, dimension(nc+buffer,nc+buffer) :: small_proj
integer, parameter :: new_nodes_dim = 4
real, dimension(nc*new_nodes_dim, nc*new_nodes_dim) :: full_proj
integer ::  axis
integer, parameter :: MSL=150
!character(MSL)::filename, zstring, ifile
real(4)          :: current_z, ScaleFactor
!character(len=40):: node_str
character(*), parameter :: output_form=    'unformatted' !'binary'
character(*), parameter :: input_access=   'stream' !'sequential' or 'stream'
character(*), parameter :: output_access=  'stream' !'sequential' or 'stream'
character(*), parameter :: input_form=     'unformatted'
character(len=MSL) :: LOS, suffix
!integer, dimension(2) :: node_shift

! do xy slabs

current_z = 0.025

do axis = 1,3

   if(axis.eq.1)suffix='xy'
   if(axis.eq.2)suffix='xz'
   if(axis.eq.3)suffix='yz'

   call get_slab

enddo

contains

subroutine get_slab 
  implicit none
  include 'project.fh'


  character(MSL)::filename, zstring, ifile
  integer, dimension(2) :: node_shift
  integer, parameter :: buffer=1
  real, dimension(nc+buffer,nc+buffer) :: small_proj
  integer ::  i, min_x, min_y, max_x, max_y, j
  character(len=40):: node_str, slab_str

!loop over all slabs to construct
do j = 1,new_nodes_dim

  write(slab_str,'(i1)') j

full_proj = 0
!construct a slab
do i=1+ (j-1)*new_nodes_dim**2,new_nodes_dim**2 + (j-1)*new_nodes_dim**2

   if(i-1<10)then  
     write(node_str,'(i1)') i-1 !z3dps       ! Need (f6.3) for Z > 10
   else
     write(node_str,'(i2)') i-1 !z3dps       ! Need (f6.3) for Z > 10
   endif

   write(zstring,'(f5.3)') current_z !z3dps       ! Need (f6.3) for Z > 10
   filename='/scratch/jharno/Lensing_2048/'//LOS//'/'//trim(zstring)//'proj_'//trim(suffix)//trim(node_str)//'.dat'
   !filename=proj_path//'/'//trim(zstring)//'proj_xy'//trim(node_str)//'.dat'
   !write(*,*) 'opening ', filename
   open(10,file=filename,form=input_form)!, access = input_access)
   read(10) ScaleFactor
   read(10)  small_proj
   close(10)

   !write(*,*) 'Mean proj =', sum(small_proj)/max(1,size(small_proj))

   node_shift(1) = mod(i-1,4)*nc
   node_shift(2) = (i-1)/4*nc - (j-1)*new_nodes_dim*nc
   !write(*,*) i, shift(1),shift(2)
   min_x = 1+node_shift(1)
   min_y = 1+node_shift(2)

   max_x = min_x + nc +buffer -1
   max_y = min_y + nc +buffer -1

   if((max_x > nc*new_nodes_dim) .and. (max_y > nc*new_nodes_dim)) then
       max_x = max_x - buffer
       max_y = max_y - buffer
       full_proj(1:buffer,1:buffer) = full_proj(1:buffer,1:buffer) + small_proj(max_x+1:max_x+buffer,max_y+1:max_y+buffer)
       full_proj(1:buffer,min_y:max_y) = full_proj(1:buffer,min_y:max_y) + small_proj(max_x+1:max_x+buffer,min_y:max_y)
       full_proj(min_x:max_x,1:buffer) = full_proj(min_x:max_x,1:buffer) + small_proj(min_x:max_x,max_y+1:max_y+buffer)
   elseif(max_x > nc*new_nodes_dim) then 
       max_x = max_x - buffer
       full_proj(1:buffer,min_y:max_y) = full_proj(1:buffer,min_y:max_y) + small_proj(max_x+1:max_x+buffer,min_y:max_y) 
   elseif(max_y > nc*new_nodes_dim) then 
       max_y = max_y - buffer
       full_proj(min_x:max_x,1:buffer) = full_proj(min_x:max_x,1:buffer) + small_proj(min_x:max_x,max_y+1:max_y+buffer)
   endif

   !write(*,*) min_x, max_x, min_y, max_y

   full_proj(min_x:max_x,min_y:max_y) = full_proj(min_x:max_x,min_y:max_y) + small_proj

enddo   

if(j==1) zstring= '0.043' 
if(j==2) zstring= '0.031' 
if(j==3) zstring= '0.018' 
if(j==4) zstring= '0.006' 



!filename='/scratch/jharno/Lensing_2048/LOS10/'//trim(zstring)//'proj_xy_'//trim(slab_str)//'.dat'
filename='/scratch/jharno/Lensing_2048/'//LOS//'/'//trim(zstring)//'proj_'//trim(suffix)//'.dat'
open(11, file=filename, status='replace', form=output_form)!, access=output_access)
write(11) ScaleFactor
write(11) full_proj
close(11)

write(*,*) 'Wrote ', filename
write(*,*) 'Mean proj =', sum(full_proj)/max(1,size(full_proj))
write(*,*) 'proj(1:10,1)' ,full_proj(1:10,1)

enddo

return
end subroutine get_slab

end program build_slabs
