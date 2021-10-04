! Code to get the global rank from halos in the light cone.
! Compile with: 
! make GetSortIDLightCone
! Run with :
! ./GetSortIDLightCone 74 /data_borussia/jharno/KiDS_Halos/ /data_borussia/jharno/KiDS_LightCone_Halos_10deg/
program getrank
  !use StringOpt
  implicit none
  include '../Lens.fh'


  character(*), parameter :: output_form=   'binary'
  character(*), parameter :: output_access= 'sequential'
  character(*), parameter :: input_form=    'binary'

  ! for 28 columns: x_peak(1:3), x_cm(4:6), v_cm(7:9), l(10:12), v_disp(13:15), radius_scale(16),
  ! halo_mass(17), halo_mass_pp(18), halo_mass1(19), var_x(20:22),I_ij(23:28)
  !----------------

  integer(4), parameter :: ncol = 28 ! choose 17 or 28
  integer(4), parameter :: nhalo_max = 6000000 ! max number of halos in sim box
  real(4), dimension(ncol,nhalo_max) :: halo_list,halo_list_out

  real, dimension(nslice) :: z_write,z_write_s
  integer i,j,k,fu1,fu2,i1,j1,j2,i3,icount, file_status, ii
  real(4), dimension(ncol) :: halo_input_buffer
  real(4), dimension(ncol+1) :: halo_output_buffer
  character (len=180) :: halofn, new_halofn, path1, path2
  character(len=7) z_string, LOS_str

  integer(4), dimension(nhalo_max)::isorthalo
  real(4), dimension(nhalo_max)::halomass!, sorted_list

  integer(4) nhalo_LC, nhalo_box 
  real(4) halomass_LC

  fu1 = 31
  fu2 = 32

  !----------------------
  ! Read redshift tables:
  !-----------
  !read z-lens
  open(11,file=fn_z)
     do i=1,nslice
         read(11,*) z_write(i)
     enddo
  close(11)
  !write(*,*)'z_lens = ', z_write

  ! read LOS
  call GetArg(1,LOS_str)
  call GetArg(2,path1)
  call GetArg(3,path2)

  path1 = adjustl(path1)
  path2 = adjustl(path2)

  do i = 1,nslice

     nhalo_box = 0
     nhalo_LC = 0 

     write(z_string,'(f5.3)') z_write(i)
     z_string = adjustl(z_string)
     write(*,*)'******************'
     write(*,*)'z_lens = ', z_string
     write(*,*)'******************'

     halofn=path1(1:len_trim(path1))//z_string(1:len_trim(z_string))//'halo.dat_LOS'//trim(LOS_str)  
     !halofn='./halo_box_file'
     write(*,*)'SimBox halo path:', halofn

     open(unit=fu2,file=halofn,form='formatted',status='old',iostat=file_status) ! for reading merged catalogs
     if (file_status.gt.0) then
        close(fu2)
        write(*,*) '***** Failed to open file ', halofn, '*******'
        write(*,*) 'file_status',file_status
     else
        write(*,*) 'Opened',  halofn
     endif

     do
        read(fu2,'(28f20.4)',end=112,err=113) halo_input_buffer
        nhalo_box = nhalo_box + 1
        if(nhalo_box > nhalo_max) then
           write(*,*) 'To many halos, increase nhalo_max'    
           nhalo_box = nhalo_box-1 
           exit
        endif

#ifdef HALO_PATCH_CORRUPTED_FILES
        !write(*,*) '************************'
        !write(*,*) '**** PATCH FOR CORRUPTED halo.dat FILES'
        !write(*,*) halo_input_buffer

        halo_input_buffer = (/halo_input_buffer(2:4), halo_input_buffer(27:28), halo_input_buffer(1), halo_input_buffer(5:26) /)
#endif

        halo_list(:,nhalo_box)=halo_input_buffer

        halomass(nhalo_box) = halo_input_buffer(17)
        !write(*,*) 'Halo mass = ', halomass(nhalo_box)
 
        ! Remove zeromass entries for 1/halomass sort
        if(halomass(nhalo_box)==0)halomass(nhalo_box) = 0.0001

113     continue
     enddo
112  close(fu2)
     write(*,*) 'Closed', halofn  

     write(*,*) '**************'
     write(*,*) 'Sorting halos '
     write(*,*) 'There are ', nhalo_box, 'halos to sort'

     !-----------------
     ! using indexsort
     isorthalo(:nhalo_box)=(/ (k,k=1,nhalo_box) /)
     print *,'finished isort_init'
     !!!!call indexedsort(nhalo_box,halomass,isorthalo)
     call indexedsort(nhalo_box,1.0/halomass(1:nhalo_box),isorthalo)
     print *,'finished isort', minval(isorthalo(1:nhalo_box))
     halomass(1:nhalo_box)=halomass(isorthalo(:nhalo_box))
 
     halo_list_out(:,:nhalo_box)=halo_list(:,isorthalo(:nhalo_box))
     print *,'finished shuffle'
     !pause

     !do ii = 1,min(nhalo_box,200)
     !   write(*,*) halo_list_out(17,ii)
     !enddo
     !pause

     !write(*,*),'***************'
     !do ii = 1,min(nhalo_box,200)
     !   write(*,*) halomass(ii), isorthalo(ii)
     !enddo
     
     !do ii = -10,min(nhalo_box,200)
     !   write(*,*) halomass(nhalo_box-ii) ii
     !enddo
     !pause

     !----------------------
     ! Get light cone halos
     !new_halofn='./halo_sort_file'
     new_halofn=path2(1:len_trim(path2))//z_string(1:len_trim(z_string))//'halo_sort.dat_LOS'//trim(LOS_str)  
     !write(*,*)'Halo sort path:', new_halofn
     !open(unit=fu1,file=new_halofn, access = output_access, form = output_form, status = 'old',iostat=file_status)
     !open(unit=fu1,file=new_halofn, form='unformatted', status='old',access='stream')
     !open (fu1,file=new_halofn, access='sequential', FORM='binary', status = 'old')
     open(unit=fu1,file=new_halofn,form='formatted',status='replace',iostat=file_status) ! for reading merged catalogs
     if (file_status.gt.0) then
        close(fu1)
        write(*,*) '***** Failed to open file ', new_halofn, '*******'
        write(*,*) 'file_status',file_status
     else
        write(*,*) 'Opened',  new_halofn
     endif
     do j=1,nhalo_box
         write(fu1,'(28f20.10)') halo_list_out(:,j)
     enddo
     close(fu1)
     print *,'finished write'

    !Read LC halo file
     !do while (nhalo_LC .le. nhalo_max)
        
     !   read(fu1,end=7) halo_input_buffer
     !   nhalo_LC = nhalo_LC + 1
     !   halomass_LC = halo_input_buffer(17)
        !write(*,*) 'halomass_LC=',halomass_LC
        !pause

        ! Find the rank of the halo.
        !do ii=1,nhalo_box
        !   !write(*,*) ii,halomass_LC, halomass(ii)
        !   if(halomass_LC .ge. halomass(ii))exit
        !enddo
        !write(*,*) nhalo_LC, halomass_LC, ii
        !pause
     !enddo

     !-----------------
     !WRITE (*,*)  'TOO MANY HALOS !!!'
     !STOP 
!7  continue
     !-----------------
!     close(fu1)

    write(*,*) 'Done z = ', z_write(i)
  enddo

  !write(*,*) 'nhalo_LC = ', nhalo_LC

  !***

end program getrank
