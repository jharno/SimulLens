! compile with ifort -shared-intel -fpp -DBINARY ourceapper.f90 -o Mapper
! Edit from pgm_proj.f90 - convert kappa maps to portable greymap format
program Mapper
  use StringOpt
implicit none
!include '../../parameters'
!#include 'CMBLens.fh'
#include 'Lens.fh'

!size of simulation mesh (nf_physical_dim in cubepm.par)

#ifdef test_read
integer, parameter :: n=nc !320
#else
integer, parameter :: n=npc !nc*2 !320
#endif

!file path to input density projections 

character(len=*),parameter :: ipath=Lens_output_path

!file path to output pgm files

character(len=*),parameter :: opath=Lens_output_path

!list of redshift to calculate pgms at (from cubepm.par)

!character(len=*),parameter :: projections='z_proj_WMAP5_140mix328Mpc'!cubepm_root//'input/projections_kappa'
character(len=*),parameter :: projections=fn_z!'Z_proj_WMAP5_328Mpc' 
!character(len=*),parameter :: projections=cubepm_root//'input/projections_test'

integer(4),parameter :: max_input=100
character(len=max_input) :: ofile,ifile
real(4),dimension(max_input) :: z_projection
integer(4) :: num_projections,cur_projection
integer(4) :: i,dim
integer(4) :: fstat
real(4) :: a
real(4), dimension(n,n) :: den 
!real(8), dimension(n,n) :: den 
integer(1), dimension(n,n) :: map

!! Read in projections to recompose

  open(11,file=projections,status='old',iostat=fstat)
  if (fstat /= 0) then
    write(*,*) 'error opening projections list file'
    write(*,*) 'file:',projections
    stop 
  endif
  do num_projections=1,max_input
    read(unit=11,err=51,end=41,fmt='(f20.10)') z_projection(num_projections)
  enddo
41  num_projections=num_projections-1
51  close(11)
  write(*,*) 'kappa maps to recompose:'
  do i=1,num_projections
    write(*,'(f5.1)') z_projection(i)
  enddo

do cur_projection=1,num_projections

!! Now start doing each projection
!do dim=1,3

  !if (dim == 1) then
#ifdef write_delta_maps
    if(z_projection(cur_projection)<10)write(ifile,'(f5.3,"delta.dat")') z_projection(cur_projection)
    if(z_projection(cur_projection)<10)write(ofile,'(f5.3,"delta.pgm")') z_projection(cur_projection)
#endif
#ifdef no_geometry
    if(z_projection(cur_projection)<10)write(ifile,'(f5.3,"kappa_no_weight.dat")') z_projection(cur_projection)
    if(z_projection(cur_projection)<10)write(ofile,'(f5.3,"kappa_no_weight.pgm")') z_projection(cur_projection)
#endif
#ifdef full_geometry
    if(z_projection(cur_projection)<10)write(ifile,'(f5.3,"kappa_weight.dat")') z_projection(cur_projection)
    if(z_projection(cur_projection)<10)write(ofile,'(f5.3,"kappa_weight.pgm")') z_projection(cur_projection)
#endif
#ifdef gamma1
    if(z_projection(cur_projection)<10)write(ifile,'(f7.3,"gamma1_weight.dat")') z_projection(cur_projection)
    if(z_projection(cur_projection)<10)write(ofile,'(f7.3,"gamma1_weight.pgm")') z_projection(cur_projection)
#endif
#ifdef gamma2
    if(z_projection(cur_projection)<10)write(ifile,'(f7.3,"gamma2_weight.dat")') z_projection(cur_projection)
    if(z_projection(cur_projection)<10)write(ofile,'(f7.3,"gamma2_weight.pgm")') z_projection(cur_projection)
#endif
#ifdef deflx
    if(z_projection(cur_projection)<10)write(ifile,'(f7.3,"deflx_weight.dat")') z_projection(cur_projection)
    if(z_projection(cur_projection)<10)write(ofile,'(f7.3,"deflx_weight.pgm")') z_projection(cur_projection)
#endif
#ifdef defly
    if(z_projection(cur_projection)<10)write(ifile,'(f7.3,"defly_weight.dat")') z_projection(cur_projection)
    if(z_projection(cur_projection)<10)write(ofile,'(f7.3,"defly_weight.pgm")') z_projection(cur_projection)
#endif


#ifdef write_delta_maps
    if(z_projection(cur_projection)>10 .and. z_projection(cur_projection)<100 )write(ifile,'(f6.3,"delta.dat")') z_projection(cur_projection)
    if(z_projection(cur_projection)>10 .and. z_projection(cur_projection)<100 )write(ofile,'(f6.3,"delta.pgm")') z_projection(cur_projection)
#endif
#ifdef no_geometry
    if(z_projection(cur_projection)>10 .and. z_projection(cur_projection)<100 )write(ifile,'(f6.3,"kappa_no_weight.dat")') z_projection(cur_projection)
    if(z_projection(cur_projection)>10 .and. z_projection(cur_projection)<100 )write(ofile,'(f6.3,"kappa_no_weight.pgm")') z_projection(cur_projection)
#endif
#ifdef full_geometry
    if(z_projection(cur_projection)>10 .and. z_projection(cur_projection)<100 )write(ifile,'(f6.3,"kappa_weight.dat")') z_projection(cur_projection)
    if(z_projection(cur_projection)>10 .and. z_projection(cur_projection)<100 )write(ofile,'(f6.3,"kappa_weight.pgm")') z_projection(cur_projection)
#endif
#ifdef gamma1
    if(z_projection(cur_projection)>10 .and. z_projection(cur_projection)<100)write(ifile,'(f7.3,"gamma1_weight.dat")') z_projection(cur_projection)
    if(z_projection(cur_projection)>10 .and. z_projection(cur_projection)<100)write(ofile,'(f7.3,"gamma1_weight.pgm")') z_projection(cur_projection)
#endif
#ifdef gamma2
    if(z_projection(cur_projection)>10 .and. z_projection(cur_projection)<100)write(ifile,'(f7.3,"gamma2_weight.dat")') z_projection(cur_projection)
    if(z_projection(cur_projection)>10 .and. z_projection(cur_projection)<100)write(ofile,'(f7.3,"gamma2_weight.pgm")') z_projection(cur_projection)
#endif
#ifdef deflx
    if(z_projection(cur_projection)>10 .and. z_projection(cur_projection)<100)write(ifile,'(f7.3,"deflx_weight.dat")') z_projection(cur_projection)
    if(z_projection(cur_projection)>10 .and. z_projection(cur_projection)<100)write(ofile,'(f7.3,"deflx_weight.pgm")') z_projection(cur_projection)
#endif
#ifdef defly
    if(z_projection(cur_projection)>10 .and. z_projection(cur_projection)<100)write(ifile,'(f7.3,"defly_weight.dat")') z_projection(cur_projection)
    if(z_projection(cur_projection)>10 .and. z_projection(cur_projection)<100)write(ofile,'(f7.3,"defly_weight.pgm")') z_projection(cur_projection)
#endif

#ifdef write_delta_maps
    if(z_projection(cur_projection)>100)write(ifile,'(f7.3,"delta.dat")') z_projection(cur_projection)
    if(z_projection(cur_projection)>100)write(ofile,'(f7.3,"delta.pgm")') z_projection(cur_projection)
#endif
#ifdef no_geometry
    if(z_projection(cur_projection)>100)write(ifile,'(f7.3,"kappa_no_weight.dat")') z_projection(cur_projection)
    if(z_projection(cur_projection)>100)write(ofile,'(f7.3,"kappa_no_weight.pgm")') z_projection(cur_projection)
#endif
#ifdef full_geometry
    if(z_projection(cur_projection)>100)write(ifile,'(f7.3,"kappa_weight.dat")') z_projection(cur_projection)
    if(z_projection(cur_projection)>100)write(ofile,'(f7.3,"kappa_weight.pgm")') z_projection(cur_projection)
#endif
#ifdef gamma1
    if(z_projection(cur_projection)>100)write(ifile,'(f7.3,"gamma1_weight.dat")') z_projection(cur_projection)
    if(z_projection(cur_projection)>100)write(ofile,'(f7.3,"gamma1_weight.pgm")') z_projection(cur_projection)
#endif
#ifdef gamma2
    if(z_projection(cur_projection)>100)write(ifile,'(f7.3,"gamma2_weight.dat")') z_projection(cur_projection)
    if(z_projection(cur_projection)>100)write(ofile,'(f7.3,"gamma2_weight.pgm")') z_projection(cur_projection)
#endif
#ifdef deflx
    if(z_projection(cur_projection)>100)write(ifile,'(f7.3,"deflx_weight.dat")') z_projection(cur_projection)
    if(z_projection(cur_projection)>100)write(ofile,'(f7.3,"deflx_weight.pgm")') z_projection(cur_projection)
#endif
#ifdef defly
    if(z_projection(cur_projection)>100)write(ifile,'(f7.3,"defly_weight.dat")') z_projection(cur_projection)
    if(z_projection(cur_projection)>100)write(ofile,'(f7.3,"defly_weight.pgm")') z_projection(cur_projection)
#endif


#ifdef test_read
    write(ifile,'(f5.3,"test.dat")') z_projection(cur_projection)
    write(ofile,'(f5.3,"test.pgm")') z_projection(cur_projection)
#endif

    print *,'opening ',ifile
  !elseif (dim == 2) then
  !  write(ifile,'(f7.3,"proj_xz.dat")') z_projection(cur_projection)
  !  write(ofile,'(f7.3,"proj_xz.pgm")') z_projection(cur_projection)
  !  print *,'calculating xz'
  !else
  !  write(ifile,'(f7.3,"proj_yz.dat")') z_projection(cur_projection)
  !  write(ofile,'(f7.3,"proj_yz.pgm")') z_projection(cur_projection)
  !  print *,'calculating yz'
  !endif
  ifile=adjustl(ifile)
  ofile=adjustl(ofile)

!#ifdef BINARY
!  open(10,file=ipath//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"Mpc_"//digitn(int(nc),4)//"/"//ifile,form='binary',iostat=fstat)
open(10,file=ipath//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(nc),4)//"_bicubic/"//ifile,form='binary',iostat=fstat)
!open(10,file=ipath//ifile,form='binary',iostat=fstat)

!#else
!  open(10,file=ipath//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"Mpc_"//digitn(int(nc),4)//"/"//ifile,form='unformatted',iostat=fstat)
!#endif
if (fstat /= 0) then
  print *,'error opening density file:',ipath//ifile
  stop
endif
!read(10) a
read(10) den 
close(10)
!print *,'scalefactor=',a,'redshift=',1.0/a-1.0
print *,'total density =',sum(den)

if(sum(den).ne.0.0)then
den=den-minval(den)
den=den/maxval(den)
endif

! one may wish to change the mapping function 

map=127*sqrt(sqrt(den))

open(10,file=ipath//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(nc),4)//"_bicubic/"//ofile,status='replace',iostat=fstat)
!open(10,file=opath//ofile,status='replace',iostat=fstat)
if (fstat /= 0) then
  print *,'error opening output pgm:',opath//ofile
  stop
endif
write(10,'(2hP5)')
write(10,*) n,n
write(10,*) 127
close(10)
!#ifdef BINARY
open(10,file=ipath//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(nc),4)//"_bicubic/"//ofile,access='append',form='binary',iostat=fstat)
!open(10,file=opath//ofile,access='append',form='binary',iostat=fstat)
!#else
!!! NOTE:  THIS IS INCORRECT AND WILL NOT WORK PROPERLY DUE TO FORTRAN RECORDS
!open(10,file=opath//ofile,access='append',form='unformatted',iostat=fstat)
!#endif
if (fstat /= 0) then
  print *,'error opening pgm file to append:',ofile
  stop
endif
write(10) map
close(10)
!enddo
enddo
end program Mapper
