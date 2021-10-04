! Lensing code written by TingTing Lu, modified by Joachim Harnois-Deraps
! to accomodate more than one box size, to compute the geometry with source
! redshifts not equal to lens redshifts -- i.e. sources are behind the box, 
! while lenses are at the mid plane. No galaxy PDF yet.
!
! Organized to read 2D slabs written by cubep3m, either proj_ij.dat or
! proj_half_ij.dat, with 'ij' = [xy,xz,yz].
!
! This code also zooms on the halo catalogue and produces a 'zoomed catalogue',
! for exact matching with lenses in the light cone.
! For halos, read Z.ZZZhalo.dat catalogues with either 17 or 28 entries per halos.
!
!
! To be run on the TCS or GPC cluster at SciNet

program SimulLens
  use StringOpt
  use Lensing
  implicit none
  include 'Lens.fh'

  integer, dimension(num_nbody_runs) :: sublist_run,count_run

!--------
#ifndef halo_only
  ! Arrays on original grid
  real, dimension(nc,nc) :: input_map 
  real, dimension(2*nc,2*nc):: zoom_map
  complex, dimension(nc,nc):: map_cplx

#ifdef calshear
  real, dimension(nc,nc) :: phi,tmp_map 
  type(Spin2), dimension(nc,nc) :: shear  
  type(vec2D), dimension(nc,nc) :: defl 
#endif

  ! Arrays on pixels
  real, dimension(npc,npc,nslice) :: map_3D
  real, dimension(npc,npc) :: cumul_kappa

#ifdef calshear
  real, dimension(npc,npc) :: cumul_gamma1,cumul_gamma2,cumul_deflx,cumul_defly,tmp_pix_map 
  type(vec2D), dimension(npc,npc,nslice) :: newdefl
  type(vec2D), dimension(npc,npc) :: CorrBornDefl
  type(Spin2), dimension(npc,npc,nslice) :: newshear
#endif
#ifdef power_spectrum
  !real, dimension(npc*2, npc*2) :: zeropad_map
#endif

#endif
!-------

  ! Other stuff
  real, dimension(npc) :: power
  real, dimension(nslice) :: z_write,z_write_s
  type(vec2D) shift
  integer i,j,fu,i1,j1,j2,i3,icount,kx,ky,ir,index_nbody_run,file_status , newLOS,patch!,redshift, z
  real lense_weight,frac,angle,dang,chi,pi,random_index_run,box,kernel ,ScaleFactor
  real(kind=8) rhomean
  character(len=180) :: fn,fn1,fp,my_status, test_str
  character(len=7) z_string, index_str, newLOS_str,LOS_str

  !-------- Multi-plateform code. For projection maps :
  ! in/output_form =  TCS: 'unformatted' ; GPC: 'binary'
  ! output_access  =  TCS: 'stream'      ; GPC: 'sequential'

  ! GPC:
  character(*), parameter :: output_form=   'binary'
  character(*), parameter :: output_access= 'sequential' 
  character(*), parameter :: input_form=    'binary' 

  ! TCS:
  !character(*), parameter :: output_form=   'unformatted'
  !character(*), parameter :: output_access= 'stream' 
  !character(*), parameter :: input_form=    'unformatted' 

  ! For halos: unmerged = input_format
  !              merged = formatted 


  ! For random shiftings
  !character(len=max_path) :: seedfile
  integer(4) :: seedsize
  integer(4), allocatable, dimension(:) :: iseed,old
  integer k_seed, clock
  integer, dimension(8) :: values
  real  rand
  integer, dimension(nslice) :: rand_LOS

  real, parameter :: z_switch = 2.51 ! 2.51!1.0 !2.4 = CMBLens ! 1.0 = CFHT Lens, 2.51 = RCSLenS/KiDS

!--------------------
#ifdef HALO_CATALOGUE
  
  !not up to date with two different boxe sizes

  !----------------
  ! for 17 columns: x_peak(1:3), x_cm(4:6), v_cm(7:9), l(10:12), v_disp(13), radius_scale(14), 
  ! halo_mass(15), halo_mass_pp(16), halo_mass1(17)
  !
  ! for 28 columns: x_peak(1:3), x_cm(4:6), v_cm(7:9), l(10:12), v_disp(13:15), radius_scale(16), 
  ! halo_mass(17), halo_mass_pp(18), halo_mass1(19), var_x(20:22),I_ij(23:28)
  !----------------

  integer(4), parameter :: ncol = 28 ! choose 17 or 28

  real(4), dimension(ncol) :: halo_input_buffer, halo_output_buffer

  !--------------
  integer nh_total, nploc(nn), ip, fu2,fu3,ii,jj,kk,ll, nh_final, nh_final_local, n_out, n_behind, n_replicate, GottaNaN,BadHalo
  integer(kind=8),dimension(10) :: halo_pid
  character (len=4) :: rank_s
  character (len=180) :: halofn, new_halofn
  real, dimension(3) ::  I_ij_vector_in, I_ij_vector_out, temp_vector
  real, dimension(3,3) ::  I_matrix, M_matrix, Rotated_I_matrix, Rotated_M_matrix
  real :: mag, chi_halo, dummy1, dummy2, theta_rot, phi_rot, psi_rot, radius_scale  
  real, dimension(3,nslice):: shake_offset

#endif
!--------------------
! Only write shift to file for kappa maps. Else, read it.

#ifdef write_kappa
  integer, parameter ::  write_shift = 0 ! Mead, fergus
  integer, parameter ::  read_shift = 1  ! Mead, fergus
  !integer, parameter ::  write_shift = 1 ! KiDS
  !integer, parameter ::  read_shift = 0  ! KiDS
#else
  integer, parameter ::  write_shift = 0 ! KiDS
  integer, parameter ::  read_shift = 1  ! KiDS
  !integer, parameter ::  write_shift = 1 ! Use these two lines for multiple random shifts 
  !integer, parameter ::  read_shift = 0
#endif

  integer, parameter ::  no_shift = 0

  ! Not working yet...
  !------------
  !method =  1: Write shift to a file
  !method = -1: Read shift from file
  !method =  0: No shift
  integer, parameter ::  method = -1
  !-----------

#ifndef halo_only
!-----------------------------------
  !equivalence assume npc >= nc
  !equivalence (cumul_kappa, input_map)
  !write(*,*) '********** Assuming npc > nc ********'

  !equivalence assume nc >= npc
  equivalence (input_map, cumul_kappa)
  write(*,*) '********** Assuming nc > npc ********'

  !equivalence (zoom_map, zeropad_map)

#ifdef calshear
!  equivalence (phi, cumul_gamma1)
#endif  
#endif

  !common angle, dang !map_3D,newshear,newdefl,CorrBornDefl



  box=lbox
  pi=acos(-1.)

  call init_var


#ifndef halo_only
  !--- Init variables here 
  map_3D(:,:,:)  =0.
#ifdef calshear
  newshear(:,:,:)%p  =0.
  newshear(:,:,:)%c  =0.
  newdefl(:,:,:)%x  =0.
  newdefl(:,:,:)%y  =0.

  CorrBornDefl(:,:)%x=0.
  CorrBornDefl(:,:)%y=0.
#endif
#endif

 !Field of view and resolution:
  !angle=5./180.*pi  ! For CMBLens
  !angle = sqrt(12.84)/180*pi  ! For UBC Lens
  angle = sqrt(Area)/180.*pi  ! For UBC Lens
  !angle=box/chi(zs,omegam,h)/h
  dang=angle/npc

  write(*,*) 'Field of view      :',(angle*180./pi)**2,'sq.degrees'
  write(*,*) 'Resolution         :',(dang*180./pi)*60,'arc min.'
  write(*,*) 'Number of z-slices :',nslice

  !-----------
  !read z-lens
  open(11,file=fn_z)  
  do i=1,nslice
     read(11,*) z_write(i)
  enddo
  close(11)
  !write(*,*)'z_lens = ', z_write

  !read z_sources
  open(11,file=fn_z_s)  
  do i=1,nslice
     read(11,*) z_write_s(i)
  enddo
  close(11)
  !write(*,*)'z_source = ', z_write_s
  !------------


  do ir=1,nr
     if(nr .lt. 1) write(*,*) 'Computing Lens', ir , 'out of', nr
     call random_seed()
     call random_seed(size=seedsize)

     allocate(iseed(seedsize))
     allocate(old(seedsize))

     call date_and_time(values=values)
     !if(rank==0) write(*,*) values(7:8), iseed      
     call random_seed(put=values(7:8)) 


     !call random_number(shift%x)
     !call random_number(shift%y)

     !icount=0 ! for testing halos straight with z = 0.042 or 0.317
     !icount=2 ! for testing halos straight with z = 0.130
     icount=1 ! for full KiDS runs or z = 0.221
     
     !Loop over the slices
     write(*,*) '****************'
     write(*,*) 'sbox = ', sbox
     write(*,*) '****************'



     !----------------------------
     ! Start redshift do-loop here
     !----------------------------

     do j=1,nslice      ! full ray tracing. Set icount to 1 above 
 
 
     !write(*,*) '*************** Doing only a few slices for testing *****************'
     !do j = 18,nslice   ! do only the first slab. set icount = 0 above 
     !do j = 16,nslice    ! do only the first three slabs. set icount = 1 above 
     !do j = 17,nslice    ! do only the first two slabs. set icount = 2 above 
     !do j = 15,nslice   ! Do only the first 4 slabs. Set icount = 0 above   
     !do j = 10,nslice   ! Do only the first 4 slabs. Set icount = 0 above   
     !do j = 1,1   ! Do only the last slabs. Set icount = X above   

        if(j==nslice_far+1)then  
           write(*,*) '****************'
           write(*,*) 'switching to lbox = ', lbox
           write(*,*) '****************'
        endif

        i=mod(icount,3)

        ! Get the Redshift
        write(z_string,'(f7.3)') z_write(j)
        z_string=adjustl(z_string)
 
        write(*,*) '******************'
        write(*,*) 'Processing z = ',  z_string

        ! for mixed size of simulations
        if(z_write(j) > z_switch) box=sbox
        if(z_write(j) <= z_switch) box=lbox


        ! Get the projection
!################################################
#ifdef mix_nbody_runs
        write(*,*) 'Mixing different N-Body runs'

        sublist_run=0
        !call random_number(random_index_run)
         call random_number(rand) ! throw away the first one
         call random_number(rand)
        !rdm_LOS = floor(rand*numLOS)        
        index_nbody_run = floor(rand*num_nbody_runs) +1
        index_nbody_run = list_run(index_nbody_run)

        if(index_nbody_run<10)then 
            index_str=digitn(index_nbody_run,1)
        elseif(index_nbody_run<100)then
            index_str=digitn(index_nbody_run,2)
        else
            index_str=digitn(index_nbody_run,3)
        endif

        write(*,*) 'Currently set to read CFTH sims.'
        if(z_write(j) <= z_switch)then 
           fn = '/scratch/jharno/Lensing_2048/LOS'//trim(index_str)//'/'//z_string(1:len_trim(z_string)) 
        endif
        if(z_write(j) >= z_switch)then 
           fn = '/scratch/jharno/NewLensing_2048/LOS'//trim(index_str)//'/'//z_string(1:len_trim(z_string)) 
        endif

        !write(*,*) 'Random Number = ',  rand, rand*num_nbody_runs, floor(rand*num_nbody_runs) +1, list_run(floor(rand*num_nbody_runs) +1)
        write(*,*) 'Selected LOS  = ',  index_nbody_run 
        rand_LOS(j) = index_nbody_run

#else        
        ! #### No Mixing Runs (default way) #####
        call GetArg(1,LOS_str)  


        !if(z_write(j) <= z_switch)fn=proj_path//trim(LOS_str)//"/"//z_string(1:len_trim(z_string)) !Scinet
        !if(z_write(j) > z_switch)fn=proj_path2//trim(LOS_str)//"/"//z_string(1:len_trim(z_string)) !
        if(z_write(j) <= z_switch)fn=proj_path//z_string(1:len_trim(z_string))                      !Jade
        if(z_write(j) > z_switch)fn=proj_path2//z_string(1:len_trim(z_string))                      ! 


#endif
!########################################

        fn=adjustl(fn)

        !if(z_write(j) < 0.050)then
           !use always the same orientation, say xy for now...
           !fp=fn(1:len_trim(fn))//'proj_xy.dat'         
        !else

        !if(i.eq.0)fp=fn(1:len_trim(fn))//'proj_half_xy.dat'     ! Scinet 
        !if(i.eq.1)fp=fn(1:len_trim(fn))//'proj_half_yz.dat'
        !if(i.eq.2)fp=fn(1:len_trim(fn))//'proj_half_xz.dat'

        !if(i.eq.0)fp=fn(1:len_trim(fn))//'proj_xy.dat'     !Scinet
        !if(i.eq.1)fp=fn(1:len_trim(fn))//'proj_yz.dat'
        !if(i.eq.2)fp=fn(1:len_trim(fn))//'proj_xz.dat'

        !if(i.eq.0)fp=fn(1:len_trim(fn))//'proj_half_xy.dat_LOS'//trim(LOS_str) ! jade
        !if(i.eq.1)fp=fn(1:len_trim(fn))//'proj_half_yz.dat_LOS'//trim(LOS_str)
        !if(i.eq.2)fp=fn(1:len_trim(fn))//'proj_half_xz.dat_LOS'//trim(LOS_str)


        newLOS_str = LOS_str      
        !LOS_str='74'
        !write(*,*) 'Enforcing to read always LOS', LOS_str

        !----------
        ! For Ian:
        !write(z_string,'(i3)') j-1
        !z_string=adjustl(z_string)
        !fn=proj_path//z_string(1:len_trim(z_string))
        !fp=proj_path//'DMONLY_L400N1024_WMAP7_cone_0_dm_box_'//z_string(1:len_trim(z_string))//'.dat '
        !write(*,*) 'Using projections from Ian McCarthy:'
        !write(*,*) 'Opening', fp

        !-------------
        ! For Mead and fergus:
        fp=fn(1:len_trim(fn))//'proj_half_finer_xy.dat'
        ! or 
        !fp=fn(1:len_trim(fn))//'proj_half_finer_xy.dat_LOS'//trim(LOS_str)
        write(*,*) 'Enforcing xy projections only'

        !-----------
        ! For KiDS : 
        !if(i.eq.0)fp=fn(1:len_trim(fn))//'proj_half_finer_xy.dat_LOS'//trim(LOS_str)
        !if(i.eq.1)fp=fn(1:len_trim(fn))//'proj_half_finer_yz.dat_LOS'//trim(LOS_str)
        !if(i.eq.2)fp=fn(1:len_trim(fn))//'proj_half_finer_xz.dat_LOS'//trim(LOS_str)

        !if(i.eq.0)write(*,*)'xy, no rotation'
        !if(i.eq.1)write(*,*)'yz, rotation'
        !if(i.eq.2)write(*,*)'xz, rotation'
        
        LOS_str = newLOS_str


        !endif

        fu=10+i

        patch=0 !0  ! -1 means don't keep track of shift. 
                     ! 0 means write shift to file


#ifndef halo_only

        !###############################################

        open(unit=fu,file=fp,form=input_form,status='old',iostat=file_status)

        !################## PATCHING ?!?!?!? ############
        if (file_status.gt.0) then 
            close(fu)
            write(*,*) 'Failed to open file ', fp
            write(*,*) 'file_status',file_status
            write(*,*) 'Trying to open another file, modifying random shifting...'

            ! Used for RCS/KiDS
            newLOS_str = '1'
            write(*,*) 'Switching ', LOS, ' for ',newLOS_str
            if(z_write(j) <= z_switch)fn=proj_path//z_string(1:len_trim(z_string))                      !Jade
            if(z_write(j) > z_switch)fn=proj_path2//z_string(1:len_trim(z_string))                      ! 

            if(i.eq.0)fp=fn(1:len_trim(fn))//'proj_half_finer_xy.dat_LOS'//trim(newLOS_str) !jade
            if(i.eq.1)fp=fn(1:len_trim(fn))//'proj_half_finer_yz.dat_LOS'//trim(newLOS_str)
            if(i.eq.2)fp=fn(1:len_trim(fn))//'proj_half_finer_xz.dat_LOS'//trim(newLOS_str)

            if(i.eq.0)write(*,*)'xy, no rotation'
            if(i.eq.1)write(*,*)'yz, rotation'
            if(i.eq.2)write(*,*)'xz, rotation'

            ! Used for CFHT organization:
            !write(*,*) 'Currently set for CFHT sims...'
            !call GetArg(2,newLOS_str)
            !if(z_write(j) <= z_switch)then
            !   fn = '/scratch/p/pen/jharno/Lensing_2048/LOS'//trim(newLOS_str)//'/'//z_string(1:len_trim(z_string))
            !endif
            !if(z_write(j) >= z_switch)then
            !   fn = '/scratch/p/pen/jharno/NewLensing_2048/LOS'//trim(newLOS_str)//'/'//z_string(1:len_trim(z_string))
            !endif

            !if(i.eq.0)fp=fn(1:len_trim(fn))//'proj_xy.dat'
            !if(i.eq.1)fp=fn(1:len_trim(fn))//'proj_yz.dat'
            !if(i.eq.2)fp=fn(1:len_trim(fn))//'proj_xz.dat'

            open(unit=fu,file=fp,form=input_form,status='old',iostat=file_status)
            write(*,*) 'Opened', fp
            patch=1 

        !################## OR NO PATCHING : ############
        else
           write(*,*) 'Opened', fp
        endif        

        !--------------
        !Read the files
        read(fu) ScaleFactor
        write(*,*) 'Got Scalefactor  = ', ScaleFactor

        !write(*,*) '******** Skipped the reading of the scale factor. Only OK for maps from McCarthy!!! ************'



        !if(z_write(j) <= z_switch) then
            read(fu)  input_map(:,:) ! rho_pxy_close
            !read(fu)  map_close(:,:,j-nslice_far) ! rho_pxy_close
            write(*,*) 'Read', fp
            close(fu)
        !else
        !    read(fu)  map_far(:,:,j) !rho_pxy_far
        !    !write(*,*) 'Read', fp
        !    close(fu)
        !endif 

        write(*,*) 'Got projection!'

        !###########################

!---------------
#ifdef test_read
        write (fn,'(f5.3,"test.dat")') z_write(j)
        open(10,file=Lens_output_path//trim(LOS_str)//"/"//fn, access = output_access, form = output_form, status = 'replace')
        !open(10,file=Lens_output_path//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(npc),4)//"/"//fn, access = output_access, form = output_form, status = 'replace')
        write(10) input_map(:,:)
        !write(10) map_close(:,:,j)
        close(10)
        write(*,*)'Wrote test Map for z=',z_write(j)
        stop
#endif
!----------------

#endif
!end of '*not* halo_only' block

!stop

!------------------------
! Randomly shift the map/halos
!------------------------
        if(write_shift.eq.1)then
           call random_number(shift%x)
           call random_number(shift%y)           

           if(patch==0)then
             !---------------------------
             ! Replace the old shift file
             if(j.eq.1)then
                open(100,file=Lens_output_path//"random_shift_LOS"//trim(LOS_str))
                !open(100,file=Lens_output_path//trim(LOS_str)//"/random_shift")
                write(100,*) shift%x, shift%y
                close(100)
             else
                open(100,file=Lens_output_path//"random_shift_LOS"//trim(LOS_str), position = 'append')
                !open(100,file=Lens_output_path//trim(LOS_str)//"/random_shift", position = 'append')
                write(100,*) shift%x, shift%y
                close(100)
             endif
             !----------------------------

           elseif(patch .eq. 1)then
              !-----------------------------------
              ! must read the other file's random shift and change it first!
              open(101,file=Lens_output_path//"/random_shift_LOS"//trim(newLOS_str), status='old')              
              do i3=1,j
                 read(101,*) shift%x, shift%y
              enddo
              close(101)

              ! Must ensure the random shift does not overlap with the orignal
              ! patch
              if (shift%x .le. 0.5) then
                  shift%x = shift%x + 0.5!0.5
              else 
                   shift%x = shift%x - 0.5!0.5
              endif
              if (shift%y .le. 0.5) then 
                  shift%y = shift%y + 0.5!0.5
              else
                  shift%y = shift%y - 0.5!0.5
              endif

              ! Now write to file
              open(100,file=Lens_output_path//"/random_shift_LOS"//trim(LOS_str), position = 'append')
              write(100,*) shift%x, shift%y
              close(100)
              !--------------------------------------

           else !patch = -1
             write(*,*) 'Not keeping track of random shift'
           endif


        elseif(read_shift.eq.1)then
           open(100,file=Lens_output_path//"random_shift_LOS"//trim(LOS_str), status='old',position='asis')
           !open(100,file=Lens_output_path//trim(LOS_str)//"/random_shift", status='old',position='asis')
           ! Sloppy way to get the right entry number
           do j2= 1,j          
              read(100,*) shift%x, shift%y
           enddo           
           close(100)

           write(*,*) 'Read random shift from file'

           !------------------------------------------------------------
           !write(*,*)'****** MODIFYING SHIFT!!!! USED ONLY FOR Mead and fergus ****'
           ! Here, the shift that is read is exactly in simulation cells units. 
           
           !First convert these in finer cells, then in fractions of nc
           shift%x = (4.0*shift%x)/real(nc)
           shift%y = (4.0*shift%y)/real(nc)
        
           !write(*,*)'Random shift before BC', shift%x, shift%y

           ! enforce BC
           if(shift%x < 0.0 ) shift%x = 1.0+shift%x
           if(shift%y < 0.0 ) shift%y = 1.0+shift%y

           !------------------------------------------------------------

        elseif(no_shift.eq.1)then
           shift%x = 0;
           shift%y = 0;
           write(*,*) 'No random shift at all'
        else
           write(*,*) 'Must specify a shift type. (No shift,  read from or write to file)'
           stop
        endif     
           
        write(*,*)'Random shift ', shift%x, shift%y
        !-----------------
        ! Got shift vector
        !-----------------        


!********************
#ifdef HALO_CATALOGUE
!********************

        ip = 0
        nh_total = 0
        nh_final = 0
        n_out = 0
        n_behind = 0
        n_replicate = 0
        !do ii=1,nn ! for un-merged halo catalogues
           ii = 1    ! for merged halo catalogues

           nh_final_local=0

           !---------------------------------------
           ! The following is for parallel catalogs
           !write(*,*) 'Reading Node ',ii           
           !write(rank_s,'(i4)') ii-1
           !rank_s=adjustl(rank_s)
           !write(zstring,'(f5.3)') z3dps       ! Need (f6.3) for Z > 10
           !halofn=fn(1:len_trim(fn))//'halo'//rank_s(1:len_trim(rank_s))//'.dat'        
           !new_halofn = Lens_output_path//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(npc),4)//"/"//z_string(1:len_trim(z_string))//'zoomed_halo'//rank_s(1:len_trim(rank_s))//".dat"
           !----------------------------------------           
           !
           ! This is for merged catalogs
           !
           halofn=halo_path//z_string(1:len_trim(z_string))//'halo.dat_LOS'//trim(LOS_str)        
           !new_halofn = '/data/jharno/KiDS_Halos_tmp/'//z_string(1:len_trim(z_string))//'zoomed_halo.dat_LOS'//trim(LOS_str)
           new_halofn = halo_path2//z_string(1:len_trim(z_string))//'zoomed_halo.dat_LOS'//trim(LOS_str)
           !new_halofn = Lens_output_path//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(npc),4)//"/"//z_string(1:len_trim(z_string))//'zoomed_halo.dat'
           !-----------------------------------------

           fu2=20+ii
           fu3=30+ii

           !open(unit=fu2,file=halofn,form=input_form,status='old',iostat=file_status) ! for read by node
           open(unit=fu2,file=halofn,form='formatted',status='old',iostat=file_status) ! for reading merged catalogs
           if (file_status.gt.0) then 
              close(fu2)
              write(*,*) 'Failed to open file ', halofn
              write(*,*) 'file_status',file_status
           else
              write(*,*) 'Opened', halofn
           endif
           open(unit=fu3,file=new_halofn, access = output_access, form = output_form, status = 'replace',iostat=file_status)
           write(*,*) 'Opened', new_halofn

           !------------          
           ! Merged halos do not have a header, 
           !------------
           ! Unmerged halos have a header, so uncomment the following:
           ! read(fu2) nploc(ii)
           ! write(fu3) nploc(ii)
           ! write(*,*) 'nploc(ii) =',nploc(ii)
           !------------


           !read halo_offset
            open(11,file='/data_borussia/jharno/KiDS_offsets/halo_offset.log_LOS'//trim(LOS_str), status='old',iostat=file_status)  
            !open(11,file='/data_borussia/jharno/KiDS_offsets/halo_offset_new.log_LOS'//trim(LOS_str), status='old',iostat=file_status)  
            if (file_status .gt.0) then
               write(*,*) 'Failed to open shake file'
               shake_offset(:,i3) = 0
               close(11)
            else
               do i3=1,nslice
               read(11,*) shake_offset(:,i3)
               write(*,*)'shake_offset = ', shake_offset(:,i3)
               enddo
               close(11)
            endif
           !------------
           !stop
      
           BadHalo = 0

           do

!#ifdef debug_halo_intense
              !write(*,*) 'Reading halo'
!#endif
#ifdef HALOPID
              !read(fu2,'(17f20.4)',end=112,err=113) halo_input_buffer, halo_pid           ! for reading merged catalogues
              read(fu2,'(28f20.4)',end=112,err=113) halo_input_buffer, halo_pid           ! for reading merged catalogues
              !read(fu2,end=112,err=113) halo_input_buffer, halo_pid                      ! for read by node
#else
              if(ncol.eq.17) then
                 read(fu2,'(17f20.4)', end=112,err=113) halo_input_buffer                    ! for reading merged catalogues, 17 columns
              elseif(ncol.eq.28)then
                 read(fu2,'(28f20.4)', end=112,err=113) halo_input_buffer                    ! for reading merged catalogues, 28 columns
              else
                 write(*,*) 'Non-valid number of columns, should be ncol = 17 or 28'
                 stop
              endif
              !read(fu2, end=112,err=113) halo_input_buffer                               ! for read by nodes
#endif        

#ifdef debug_halo_intense
              write(*,*) 'Got halo'
#endif

              !----- Look for NaN
              GottaNaN = 0
              do jj = 1,ncol 
                 if(isnan(halo_input_buffer(jj)) .or. (halo_input_buffer(jj)) .gt. HUGE(0.0) .or. halo_input_buffer(jj) .lt. -HUGE(0.0)) then
                    GottaNaN = GottaNaN + 1; 
                 endif
              enddo

              if(GottaNaN .gt.0)then
                 write(*,*) 'Caught a NaN in the input halo.dat file'!, halo mass=',halo_input_buffer(17)
                 BadHalo = BadHalo+1
                 !exit              
                 cycle
              elseif(maxval(halo_input_buffer)==0) then
!#ifdef debug_halo
                 write(*,*) 'No halos in  node', ii 
                 write(*,*) 'nh_total =' , nh_total
!#endif
                 !cycle
                 exit
              else
                 nh_total=nh_total+1
              endif

              !-------
              !write(*,*)' Undo the shake offset from the halo finder:', shake_offset(1:3,j)
              ! Only for the first node: first quadrant: 0 <= x,y,z < 3072/4 = 768
              !-------
              if(halo_input_buffer(1).lt. 768.0 .and. halo_input_buffer(2) .lt.768.0 .and. halo_input_buffer(3) .lt. 768.0) then
                 halo_input_buffer(1:3) = halo_input_buffer(1:3) + shake_offset(1:3,j)
                 halo_input_buffer(4:6) = halo_input_buffer(4:6) + shake_offset(1:3,j)
              endif

              !-----------
              ! Correct for the bugged  I_zz
              ! Empirical coorection:
              ! mp = 2.8792e9

              halo_input_buffer(28) = (sum(halo_input_buffer(20:21))*(halo_input_buffer(18)/8.0 -1) - 3.3e-19*(halo_input_buffer(17)*2.8792e9)**(1.61))/2.0 + 3.3e-19*(halo_input_buffer(17)*2.8792e9)**(1.61)



              !  Tricks for debug and test: 
              !-----
              !write(*,*) 'For quick testing, keep only 2000 halos'
              !if(nh_total>2000)exit
              !----
              !
              !if(nh_final_local .gt. 1718600) then
              !    write(*,*) 'Entering halo', nh_final_local+1
              !    write(*,*) halo_input_buffer(:)
              !endif
              !if(nh_final_local .lt. 1718600) then
              !    nh_final_local = nh_final_local + 1
              !    cycle
#ifdef debug_halo
                  write(*,*) '**** Entering halo', nh_final_local+1, '****'
              !    !write(*,*) halo_input_buffer(:)
#endif
              !else
              !    write(*,*) 'Entering halo', nh_final_local
              !    pause
              !endif

              
              !********************************
              ! 1- Rotate the halo catalogues *
              !********************************
#ifdef debug_halo_intense
!                 write(*,*)'Before rotation:'
                 write(*,*)'halo_input_buffer(4:6)', halo_input_buffer(4:6)
#endif                 
              
              !----------------------------------------------------------------------------------------
              ! I_ij = m*sum[ r**2 * kronecker_delta_ij -x_i*x_j ] = Tr(M) * krnoecker_delta_ij - M_ij 
              ! M_ij = m*sum[ x_i*x_j ] = 0.5*Tr(I)*kronecker_delta_ij - I_ij 

              I_matrix(1,1:3) = halo_input_buffer(23:25)
              I_matrix(2,2:3) = halo_input_buffer(26:27)
              I_matrix(3,3)   = halo_input_buffer(28)
              I_matrix(2,1)   = I_matrix(1,2)
              I_matrix(3,1)   = I_matrix(1,3)
              I_matrix(3,2)   = I_matrix(2,3)

              M_matrix = - I_matrix
              M_matrix(1,1) = M_matrix(1,1) + 0.5*(I_matrix(1,1) + I_matrix(2,2) + I_matrix(3,3))
              M_matrix(2,2) = M_matrix(2,2) + 0.5*(I_matrix(1,1) + I_matrix(2,2) + I_matrix(3,3))
              M_matrix(3,3) = M_matrix(3,3) + 0.5*(I_matrix(1,1) + I_matrix(2,2) + I_matrix(3,3))

              !write(*,*) '******* I_matrix ********'
              !write(*,*) I_matrix(1:3,1)
              !write(*,*) I_matrix(1:3,2)
              !write(*,*) I_matrix(1:3,3)

              !write(*,*) '******* M_matrix ********'
              !write(*,*) M_matrix(1:3,1)
              !write(*,*) M_matrix(1:3,2)
              !write(*,*) M_matrix(1:3,3)
              !pause

              !***************
              if(i .eq. 0)then  ! proj_xy
                 !no rotation
                 !write(*,*) 'xy-plane'
                 halo_output_buffer = halo_input_buffer
              endif

              !***************
              if(i .eq. 1)then  ! proj_yz      
                 !rotate these vectors about y axis by pi/2, clockwise:
                 !x_peak, x_cm, v_cm, l, and, for 28 columns: v_disp 3D, var_x, I_ij
                 theta_rot = pi/2.0
                 phi_rot = 0.0
                 psi_rot = 0.0

                 !write(*,*) 'yz-plane'

                 ! go to yz plane, but missing a second rotation to match with correct xy projections
                 call Rotate(halo_input_buffer(1:3),phi_rot,theta_rot,psi_rot,halo_output_buffer(1:3) )     !x_peak
                 call Rotate(halo_input_buffer(4:6),phi_rot,theta_rot,psi_rot,halo_output_buffer(4:6) )     !x_cm
                 call Rotate(halo_input_buffer(7:9),phi_rot,theta_rot,psi_rot,halo_output_buffer(7:9) )     !v_cm
                 call Rotate(halo_input_buffer(10:12),phi_rot,theta_rot,psi_rot,halo_output_buffer(10:12) ) !l

                 !&&&&&&&&&&   The following if for 28 columns: &&&&&&
                 if(ncol.eq.28) then
                    call Rotate(halo_input_buffer(13:15),phi_rot,theta_rot,psi_rot,halo_output_buffer(13:15) ) !v_disp 3D
                    call Rotate(halo_input_buffer(20:22),phi_rot,theta_rot,psi_rot,halo_output_buffer(20:22) ) !var_x
                    halo_output_buffer(20:22) = abs(halo_output_buffer(20:22))                   !take absolute value
                    halo_output_buffer(13:15) = abs(halo_output_buffer(13:15))                   !Added 11/2013


                    !------------------------------------------------------------------------------
                    ! For the I_ij matrix, I need to rotate the M_matrix first, then solve for I_ij
                    ! Rotated_M_ij = sum_alpha sum_beta R_i_alpha R_j_beta M_alpha_beta 

                    ! do x:
                    call Rotate(M_matrix(1:3,1), phi_rot, theta_rot, psi_rot, I_ij_vector_out)
                    temp_vector(1) = I_ij_vector_out(1)
                    call Rotate(M_matrix(1:3,2), phi_rot, theta_rot, psi_rot, I_ij_vector_out)
                    temp_vector(2) = I_ij_vector_out(1)
                    call Rotate(M_matrix(1:3,3), phi_rot, theta_rot, psi_rot, I_ij_vector_out)
                    temp_vector(3) = I_ij_vector_out(1)

                    call Rotate(temp_vector, phi_rot, theta_rot, psi_rot, Rotated_M_matrix(1:3,1))
                    
                    ! do y:
                    call Rotate(M_matrix(1:3,1), phi_rot, theta_rot, psi_rot, I_ij_vector_out)
                    temp_vector(1) = I_ij_vector_out(2)
                    call Rotate(M_matrix(1:3,2), phi_rot, theta_rot, psi_rot, I_ij_vector_out)
                    temp_vector(2) = I_ij_vector_out(2)
                    call Rotate(M_matrix(1:3,3), phi_rot, theta_rot, psi_rot, I_ij_vector_out)
                    temp_vector(3) = I_ij_vector_out(2)

                    call Rotate(temp_vector, phi_rot, theta_rot, psi_rot, Rotated_M_matrix(1:3,2))

                    ! do z:
                    call Rotate(M_matrix(1:3,1), phi_rot, theta_rot, psi_rot, I_ij_vector_out)
                    temp_vector(1) = I_ij_vector_out(3)
                    call Rotate(M_matrix(1:3,2), phi_rot, theta_rot, psi_rot, I_ij_vector_out)
                    temp_vector(2) = I_ij_vector_out(3)
                    call Rotate(M_matrix(1:3,3), phi_rot, theta_rot, psi_rot, I_ij_vector_out)
                    temp_vector(3) = I_ij_vector_out(3)

                    call Rotate(temp_vector, phi_rot, theta_rot, psi_rot, Rotated_M_matrix(1:3,3))


                    !write(*,*) '******* Rotated M_matrix ********'
                    !write(*,*) Rotated_M_matrix(1:3,1)
                    !write(*,*) Rotated_M_matrix(1:3,2)
                    !write(*,*) Rotated_M_matrix(1:3,3)
                    !pause
                    
                 endif
                 !&&&&&&&&&&& end of 28 columns block &&&&&&&&&&

                 !rotate by pi/2 clockwise about the new z-axis to match the projections

                 theta_rot = 0.0
                 phi_rot = pi/2.0
                 psi_rot = 0.0
                 call Rotate(halo_output_buffer(1:3),phi_rot,theta_rot,psi_rot,halo_output_buffer(1:3) )     !x_peak
                 call Rotate(halo_output_buffer(4:6),phi_rot,theta_rot,psi_rot,halo_output_buffer(4:6) )     !x_cm
                 call Rotate(halo_output_buffer(7:9),phi_rot,theta_rot,psi_rot,halo_output_buffer(7:9) )     !v_cm
                 call Rotate(halo_output_buffer(10:12),phi_rot,theta_rot,psi_rot,halo_output_buffer(10:12) ) !l


                 !&&&&&&&&&&   The following if for 28 columns: &&&&&&
                 if(ncol.eq.28) then
                    call Rotate(halo_output_buffer(13:15),phi_rot,theta_rot,psi_rot,halo_output_buffer(13:15) ) !v_disp 3D
                    call Rotate(halo_output_buffer(20:22),phi_rot,theta_rot,psi_rot,halo_output_buffer(20:22) ) !var_x
                    halo_output_buffer(20:22) = abs(halo_output_buffer(20:22))                   !take absolute value
                    halo_output_buffer(13:15) = abs(halo_output_buffer(13:15))                   !Added 11/2013


                    !------------------------------------------------------------------------------
                    ! For the I_ij matrix, I need to rotate the M_matrix first, then solve for I_ij
                    ! Rotated_M_ij = sum_alpha sum_beta R_i_alpha R_j_beta M_alpha_beta 

                    M_matrix = Rotated_M_matrix

                    ! do x:
                    call Rotate(M_matrix(1:3,1), phi_rot, theta_rot, psi_rot, I_ij_vector_out)
                    temp_vector(1) = I_ij_vector_out(1)
                    call Rotate(M_matrix(1:3,2), phi_rot, theta_rot, psi_rot, I_ij_vector_out)
                    temp_vector(2) = I_ij_vector_out(1)
                    call Rotate(M_matrix(1:3,3), phi_rot, theta_rot, psi_rot, I_ij_vector_out)
                    temp_vector(3) = I_ij_vector_out(1)

                    call Rotate(temp_vector, phi_rot, theta_rot, psi_rot, Rotated_M_matrix(1:3,1))
                    
                    ! do y:
                    call Rotate(M_matrix(1:3,1), phi_rot, theta_rot, psi_rot, I_ij_vector_out)
                    temp_vector(1) = I_ij_vector_out(2)
                    call Rotate(M_matrix(1:3,2), phi_rot, theta_rot, psi_rot, I_ij_vector_out)
                    temp_vector(2) = I_ij_vector_out(2)
                    call Rotate(M_matrix(1:3,3), phi_rot, theta_rot, psi_rot, I_ij_vector_out)
                    temp_vector(3) = I_ij_vector_out(2)

                    call Rotate(temp_vector, phi_rot, theta_rot, psi_rot, Rotated_M_matrix(1:3,2))

                    ! do z:
                    call Rotate(M_matrix(1:3,1), phi_rot, theta_rot, psi_rot, I_ij_vector_out)
                    temp_vector(1) = I_ij_vector_out(3)
                    call Rotate(M_matrix(1:3,2), phi_rot, theta_rot, psi_rot, I_ij_vector_out)
                    temp_vector(2) = I_ij_vector_out(3)
                    call Rotate(M_matrix(1:3,3), phi_rot, theta_rot, psi_rot, I_ij_vector_out)
                    temp_vector(3) = I_ij_vector_out(3)

                    call Rotate(temp_vector, phi_rot, theta_rot, psi_rot, Rotated_M_matrix(1:3,3))


                    !write(*,*) '******* Rotated M_matrix ********'
                    !write(*,*) Rotated_M_matrix(1:3,1)
                    !write(*,*) Rotated_M_matrix(1:3,2)
                    !write(*,*) Rotated_M_matrix(1:3,3)
                    

                    ! I_ij = -M_ij + Tr(M)*kronecker_delta_ij
                    Rotated_I_matrix = - Rotated_M_matrix
                    Rotated_I_matrix(1,1) = Rotated_I_matrix(1,1) + Rotated_M_matrix(1,1) + Rotated_M_matrix(2,2) + Rotated_M_matrix(3,3)
                    Rotated_I_matrix(2,2) = Rotated_I_matrix(2,2) + Rotated_M_matrix(1,1) + Rotated_M_matrix(2,2) + Rotated_M_matrix(3,3)
                    Rotated_I_matrix(3,3) = Rotated_I_matrix(3,3) + Rotated_M_matrix(1,1) + Rotated_M_matrix(2,2) + Rotated_M_matrix(3,3)

 
                    !write(*,*) '******* Rotated I_matrix ********'
                    !write(*,*) Rotated_I_matrix(1:3,1)
                    !write(*,*) Rotated_I_matrix(1:3,2)
                    !write(*,*) Rotated_I_matrix(1:3,3)
                    !pause
                    
                    halo_output_buffer(23) = Rotated_I_matrix(1,1) !I_xx
                    halo_output_buffer(24) = Rotated_I_matrix(1,2) !I_xy
                    halo_output_buffer(25) = Rotated_I_matrix(1,3) !I_xz
                    halo_output_buffer(26) = Rotated_I_matrix(2,2) !I_yy
                    halo_output_buffer(27) = Rotated_I_matrix(2,3) !I_yz
                    halo_output_buffer(28) = Rotated_I_matrix(3,3) !I_zz

                 endif
                 !&&&&&&&&&&& end of 28 columns block &&&&&&&&&&


                 !---- scalars: ------
                 !Halo radius and masses (28 colums):

                 if(ncol.eq.28) then
                    halo_output_buffer(16:19) =  halo_input_buffer(16:19)
                 else 
                    halo_output_buffer(13:17) =  halo_input_buffer(13:17)
                 endif
                 !--------------------

              endif

              !***************
              if(i .eq. 2)then      ! proj_xz
                 !rotate about x axis by pi/2, clockwise
                 theta_rot = -pi/2.0
                 phi_rot = pi/2.0
                 psi_rot = -pi/2.0

                 !write(*,*) 'xz-plane'

                 call Rotate(halo_input_buffer(1:3),phi_rot,theta_rot,psi_rot,halo_output_buffer(1:3) )     !x_peak
                 call Rotate(halo_input_buffer(4:6),phi_rot,theta_rot,psi_rot,halo_output_buffer(4:6) )     !x_cm
                 call Rotate(halo_input_buffer(7:9),phi_rot,theta_rot,psi_rot,halo_output_buffer(7:9) )     !v_cm
                 call Rotate(halo_input_buffer(10:12),phi_rot,theta_rot,psi_rot,halo_output_buffer(10:12) ) !l_cm

                 !&&&&&&&&&&   The following if for 28 columns: &&&&&&
                 if(ncol.eq.28) then
                    call Rotate(halo_input_buffer(13:15),phi_rot,theta_rot,psi_rot,halo_output_buffer(13:15) ) !v_disp 3D
                    call Rotate(halo_input_buffer(20:22),phi_rot,theta_rot,psi_rot,halo_output_buffer(20:22) ) !var_x
                    halo_output_buffer(20:22) = abs(halo_output_buffer(20:22))                   !take absolute value
                    halo_output_buffer(13:15) = abs(halo_output_buffer(13:15))                   !Added 11/2013


                    !------------------------------------------------------------------------------
                    ! For the I_ij matrix, I need to rotate the M_matrix first, then solve for I_ij
                    ! Rotated_M_ij = sum_alpha sum_beta R_i_alpha R_j_beta M_alpha_beta 

                    ! do x:
                    call Rotate(M_matrix(1:3,1), phi_rot, theta_rot, psi_rot, I_ij_vector_out)
                    temp_vector(1) = I_ij_vector_out(1)
                    call Rotate(M_matrix(1:3,2), phi_rot, theta_rot, psi_rot, I_ij_vector_out)
                    temp_vector(2) = I_ij_vector_out(1)
                    call Rotate(M_matrix(1:3,3), phi_rot, theta_rot, psi_rot, I_ij_vector_out)
                    temp_vector(3) = I_ij_vector_out(1)

                    call Rotate(temp_vector, phi_rot, theta_rot, psi_rot, Rotated_M_matrix(1:3,1))
                    
                    ! do y:
                    call Rotate(M_matrix(1:3,1), phi_rot, theta_rot, psi_rot, I_ij_vector_out)
                    temp_vector(1) = I_ij_vector_out(2)
                    call Rotate(M_matrix(1:3,2), phi_rot, theta_rot, psi_rot, I_ij_vector_out)
                    temp_vector(2) = I_ij_vector_out(2)
                    call Rotate(M_matrix(1:3,3), phi_rot, theta_rot, psi_rot, I_ij_vector_out)
                    temp_vector(3) = I_ij_vector_out(2)

                    call Rotate(temp_vector, phi_rot, theta_rot, psi_rot, Rotated_M_matrix(1:3,2))

                    ! do z:
                    call Rotate(M_matrix(1:3,1), phi_rot, theta_rot, psi_rot, I_ij_vector_out)
                    temp_vector(1) = I_ij_vector_out(3)
                    call Rotate(M_matrix(1:3,2), phi_rot, theta_rot, psi_rot, I_ij_vector_out)
                    temp_vector(2) = I_ij_vector_out(3)
                    call Rotate(M_matrix(1:3,3), phi_rot, theta_rot, psi_rot, I_ij_vector_out)
                    temp_vector(3) = I_ij_vector_out(3)

                    call Rotate(temp_vector, phi_rot, theta_rot, psi_rot, Rotated_M_matrix(1:3,3))


                    !write(*,*) '******* Rotated M_matrix ********'
                    !write(*,*) Rotated_M_matrix(1:3,1)
                    !write(*,*) Rotated_M_matrix(1:3,2)
                    !write(*,*) Rotated_M_matrix(1:3,3)
                 
                 endif
                 !&&&&&&&&&&& end of 28 columns block &&&&&&&&&&

                 !write(*,*) 'Flipping the new z-axis'
                 halo_output_buffer(3) = real(nc)/4. + halo_output_buffer(3)
                 halo_output_buffer(6) = real(nc)/4. + halo_output_buffer(6)
                 !halo_output_buffer(9) = real(nc)/4. + halo_output_buffer(9)    !V_z and L_z (about CM) are translation invariant 
                 !halo_output_buffer(12) = real(nc)/4. + halo_output_buffer(12)  !Bug fixed on Nov. 25 2013.

                 ! For 28 columns:
                 ! var and M_(i==j)  are even under flipping z
                 ! but M_xz and M_yz must flip sign:
                 if(ncol.eq.28) then
                     Rotated_M_matrix(1,3) = -Rotated_M_matrix(1,3) 
                     Rotated_M_matrix(2,3) = -Rotated_M_matrix(2,3) 
                     Rotated_M_matrix(3,1) = -Rotated_M_matrix(3,1) 
                     Rotated_M_matrix(3,2) = -Rotated_M_matrix(3,2) 

                     !write(*,*) '******* Rotated M_matrix ********'
                     !write(*,*) Rotated_M_matrix(1:3,1)
                     !write(*,*) Rotated_M_matrix(1:3,2)
                     !write(*,*) Rotated_M_matrix(1:3,3)
                    

                     ! I_ij = -M_ij + Tr(M)*kronecker_delta_ij
                     Rotated_I_matrix = - Rotated_M_matrix
                     Rotated_I_matrix(1,1) = Rotated_I_matrix(1,1) + Rotated_M_matrix(1,1) + Rotated_M_matrix(2,2) + Rotated_M_matrix(3,3)
                     Rotated_I_matrix(2,2) = Rotated_I_matrix(2,2) + Rotated_M_matrix(1,1) + Rotated_M_matrix(2,2) + Rotated_M_matrix(3,3)
                     Rotated_I_matrix(3,3) = Rotated_I_matrix(3,3) + Rotated_M_matrix(1,1) + Rotated_M_matrix(2,2) + Rotated_M_matrix(3,3)

 
                     !write(*,*) '******* Rotated I_matrix ********'
                     !write(*,*) Rotated_I_matrix(1:3,1)
                     !write(*,*) Rotated_I_matrix(1:3,2)
                     !write(*,*) Rotated_I_matrix(1:3,3)
                     !pause
                    
                     halo_output_buffer(23) = Rotated_I_matrix(1,1) !I_xx
                     halo_output_buffer(24) = Rotated_I_matrix(1,2) !I_xy
                     halo_output_buffer(25) = Rotated_I_matrix(1,3) !I_xz
                     halo_output_buffer(26) = Rotated_I_matrix(2,2) !I_yy
                     halo_output_buffer(27) = Rotated_I_matrix(2,3) !I_yz
                     halo_output_buffer(28) = Rotated_I_matrix(3,3) !I_zz
                 endif


                 !---- scalars: ------
                 !Halo radius and masses:
                 if(ncol.eq.28) then
                    halo_output_buffer(16:19) =  halo_input_buffer(16:19)
                    radius_scale=halo_output_buffer(16)
                 else
                    halo_output_buffer(13:17) =  halo_input_buffer(13:17)
                    radius_scale=halo_output_buffer(14)
                 endif
                 !--------------------

              endif
              ! end of rotation block 
#ifdef debug_halo_intense
              write(*,*)'after rotation'
              write(*,*)'halo_output_buffer(4:6)', halo_output_buffer(4:6)
#endif




              !------------------------------------------------------ 
              ! Apply random shift to position, then use periodic BC:             
              do jj = 1,6

                 if(jj==1 .or. jj==4)halo_output_buffer(jj) = halo_output_buffer(jj) - shift%x*real(nc)/4.0  ! Divide nc by 4 for finer proj.
                 if(jj==2 .or. jj==5)halo_output_buffer(jj) = halo_output_buffer(jj) - shift%y*real(nc)/4.0 

                 if((jj==1 .or. jj==2 .or. jj==4 .or. jj==5) .and. halo_output_buffer(jj) < -15.0) then
                    halo_output_buffer(jj) = halo_output_buffer(jj) + nc/4 ! 1/4 for finer projections
                 elseif((jj==3 .or. jj==6) .and. halo_output_buffer(jj) < 0.0) then
                    halo_output_buffer(jj) = halo_output_buffer(jj) + nc/4 ! 1/4 for finer projections
                 elseif(halo_output_buffer(jj) >nc/4) then                 ! idem  ! This case should not arise since we subtract a number...
                    halo_output_buffer(jj) = halo_output_buffer(jj) - nc/4 ! idem
                 endif
              enddo

              ! Now all haloes have their (x,y,z) coordinates contained in the range [1,nc/4]
              !------------------------------------------------------ 

#ifdef debug_halo_intense
              write(*,*)'after shift and mod:'
              write(*,*)'halo_output_buffer(4:6)', halo_output_buffer(4:6)
#endif

              !*********************************
              ! 2- Find magnification and zoom *
              !*********************************

              ! Get chi halo = chi lens + local z coordinate
              !chi_halo = chi(z_write(j),omegam,h)*h + (halo_output_buffer(6) - nc/2)*(box/nc)         ! full box projected (1/2)
              chi_halo = chi(z_write(j),omegam,h)*h + (halo_output_buffer(6) - (nc/4)/4.0)*(box/(nc/4))  ! half box projected (1/4)
                                                                                                       ! extra 1/4 for finer projections

#ifdef debug_halo_intense
              write(*,*) 'chi(z_s) =',chi(z_write_s(j),omegam,h)*h, 'Mpc/h'
              write(*,*) 'chi(z_i) =',chi(z_write(j),omegam,h)*h, 'Mpc/h'
              write(*,*) 'chi_halo =',chi_halo, 'Mpc/h'
#endif

              ! Do we keep or reject halos behind the source?
              !if(chi_halo>chi(zs,omegam,h)*h) then
              !   !write(*,*)'Halo behind source'
              !   cycle
              !endif


              ! ------------
              ! Reject halos on the upper half of the box in 1/2 box projections (xy and yz projection ) or lower half (xz proj)
              if(i.lt.2) then
                 !if(halo_output_buffer(6)>(nc/4)/2) then ! extra 1/4 for finer projections

                 !***allow for large halos just on the other side of both half box edge to contribute, using radius_calc. 14th or 16th column ***
                 if((halo_output_buffer(6) - radius_scale .gt. (nc/4)/2) .and. &
                        (halo_output_buffer(6) + radius_scale .lt. nc/4)) then ! extra 1/4 for finer projections
#ifdef debug_halo_intense
                    write(*,*)'Halo in upper half:', halo_output_buffer(6), 'Boundary is ', (nc/4)/2, 'cells'
                    write(*,*)'Halo in upper half:', chi_halo, chi(z_write_s(j),omegam,h)*h, 'Mpc/h' 
#endif                    
                    n_behind=n_behind+1
                    cycle
                 endif
              elseif(i.eq.2)then ! keep only upper half for this rotation
                 !if(halo_output_buffer(6)<(nc/4)/2) then ! extra 1/4 for finer projections
                 if((halo_output_buffer(6) + radius_scale .lt. (nc/4)/2) .and. & 
                        (halo_output_buffer(6) - radius_scale .gt. 0 )) then ! extra 1/4 for finer projections
#ifdef debug_halo_intense
                    write(*,*)'Halo in lower half:', halo_output_buffer(6), (nc/4)/2
                    write(*,*)'Halo in lower half:', chi_halo, chi(z_write_s(j),omegam,h)*h, 'Mpc/h' 
#endif
                    n_behind=n_behind+1
                    cycle
                 endif

              endif
              !-------------

  
              !*** Conical geometry ***!
              ! Interpolate npc pixels in a continuous conical volume, 
              ! Find distance to center of lens, then zoom on that volume
              !mag = ((angle/2)*chi_halo*(nc/box))

              !*** Blocks geometry ***
              ! Interpolate npc pixels in a volume that matches the projected
              ! densities. Use the lens distance for all halos.

              ! Not sure if that is good...
              !mag =  ((angle/2)*chi(z_write(j),omegam,h)*h*((nc/4)/box)) !in grid cells units. Extra 1/4 for finer projections

              frac=angle/(box / chi(z_write(j),omegam,h) / h) 
              ! halos with x>frac*nc are excluded,
              ! others are rescaled onto npc pixels.
              mag = frac * (real(nc)/4.0) ! 1/4 because of finer projections

              ! Apply empirical correction factor
              !write(*,*) 'mag=',mag
              mag = mag/(1.0+132.6332e-12*(chi(z_write(j), omegam,h)*h)**2 - 1.2417e-6*chi(z_write(j), omegam, h)*h + 4.0923e-3) 
              !write(*,*) 'mag=',mag, z_write(j), chi(z_write(j),omegam,h) * h
 
              !stop

#ifdef debug_halo_intense
              write(*,*) 'frac = ',frac
              write(*,*) 'mag = ',mag, 'grid cells'              
              write(*,*)'x_before mag =',halo_output_buffer(4) 
              write(*,*)'y_before mag =',halo_output_buffer(5) 
              !write(*,*)'var_x_before mag=',halo_output_buffer(20)
              !write(*,*)'var_y_before mag=',halo_output_buffer(21)
#endif



              !*********************************************************
              ! 3- Shift [x,y]  w.r.t centre of lens-plane at [x,y]=[hc,hc]
              !  - Zoom and shift back [x,y]
              !  - Zoom on [x,y] velocity, v_disp, radius_calc,l_x,l_y *
              !  - Zoom**2 on var_x,var_y, l_z               
              !*********************************************************

              !*** Conical geometry ***!
              !halo_output_buffer(1)  = hc/(hc - mag)*(halo_output_buffer(1) - hc) + hc 
              !halo_output_buffer(2)  = hc/(hc - mag)*(halo_output_buffer(2) - hc) + hc  
              !halo_output_buffer(4)  = hc/(hc - mag)*(halo_output_buffer(4) - hc) + hc 
              !halo_output_buffer(5)  = hc/(hc - mag)*(halo_output_buffer(5) - hc) + hc  

              !halo_output_buffer(7)  = hc/(hc - mag)*halo_output_buffer(7)
              !halo_output_buffer(8)  = hc/(hc - mag)*halo_output_buffer(8)
              !halo_output_buffer(10) = hc/(hc - mag)*halo_output_buffer(10)
              !halo_output_buffer(11) = hc/(hc - mag)*halo_output_buffer(11)
              !halo_output_buffer(13) = hc/(hc - mag)*halo_output_buffer(13) ! not sure though... v_disp is more complicated...
              !halo_output_buffer(14) = hc/(hc - mag)*halo_output_buffer(14)

              !halo_output_buffer(12) = (hc/(hc - mag))**2*halo_output_buffer(12)
              !halo_output_buffer(18) = (hc/(hc - mag))**2*halo_output_buffer(18)
              !halo_output_buffer(19) = (hc/(hc - mag))**2*halo_output_buffer(19)
             

               
              !*** Blocks geometry ***
              !halo_output_buffer(1)  = hc/(mag)*(halo_output_buffer(1) - hc) + hc
              !halo_output_buffer(2)  = hc/(mag)*(halo_output_buffer(2) - hc) + hc
              !halo_output_buffer(4)  = hc/(mag)*(halo_output_buffer(4) - hc) + hc
              !halo_output_buffer(5)  = hc/(mag)*(halo_output_buffer(5) - hc) + hc              

              !halo_output_buffer(7)  = hc/(mag)*halo_output_buffer(7)
              !halo_output_buffer(8)  = hc/(mag)*halo_output_buffer(8)
              !halo_output_buffer(10) = hc/(mag)*halo_output_buffer(10)
              !halo_output_buffer(11) = hc/(mag)*halo_output_buffer(11)
              !halo_output_buffer(13) = hc/(mag)*halo_output_buffer(13) ! not sure though... v_disp is more complicated...
              !halo_output_buffer(14) = hc/(mag)*halo_output_buffer(14)

              !halo_output_buffer(16) = (hc/(mag))**2*halo_output_buffer(16)

              !halo_output_buffer(20) = (hc/(mag))**2*halo_output_buffer(20)
              !halo_output_buffer(21) = (hc/(mag))**2*halo_output_buffer(21)

               ! interpolate onto zoom pixels the positions, but leave the other kinematical variables in simulation units.
               halo_output_buffer(1) = halo_output_buffer(1)*real(npc)/mag
               halo_output_buffer(2) = halo_output_buffer(2)*real(npc)/mag
               halo_output_buffer(4) = halo_output_buffer(4)*real(npc)/mag
               halo_output_buffer(5) = halo_output_buffer(5)*real(npc)/mag

#ifdef debug_halo_intense
              write(*,*)'x_after mag=',halo_output_buffer(4) 
              write(*,*)'y_after mag=',halo_output_buffer(5) 
              !write(*,*)'var_x_after=',halo_output_buffer(20)
              !write(*,*)'var_y_after=',halo_output_buffer(21)              
#endif

! If halo peak is inside the light cone, allowing for 15 pixel buffers, write. Else, don't write!              
              if((halo_output_buffer(1) < npc +15.) .and. (halo_output_buffer(1).ge.-15.) .and. &
                   (halo_output_buffer(2) < npc + 15.) .and. (halo_output_buffer(2).ge.-15.))then
              !if((halo_output_buffer(4) < npc) .and. (halo_output_buffer(4).ge.0) .and. &
              !     (halo_output_buffer(5) < npc) .and. (halo_output_buffer(5).ge.0))then

                 nh_final_local = nh_final_local+1

#ifdef HALOPID
                 write(fu3) halo_output_buffer, halo_pid!(:,nh_total)           
#else
                 write(fu3) halo_output_buffer
#endif        

                 !----------------------
                 !if frac > X.X, replicate in 8 regions with periodic BC: 
                 ! x+, x+y+, y+, y+x-, x-, x-y-, y-, x+y-
                 if (frac>0.0) then


                    ! 1- 
                    ! Replicate along x+ axis:

                    halo_output_buffer(1) = halo_output_buffer(1) + (nc/4)*npc/mag 
                    halo_output_buffer(4) = halo_output_buffer(4) + (nc/4)*npc/mag 

                    ! write to file
                    ! if(halo_output_buffer(1) .le. npc + 15. .and. halo_output_buffer(2) .le. npc + 15.) then
                   if((halo_output_buffer(1) < npc +15.) .and. (halo_output_buffer(1).ge.-15.) .and. &
                   (halo_output_buffer(2) < npc + 15.) .and. (halo_output_buffer(2).ge.-15.))then
#ifdef HALOPID
                       write(fu3) halo_output_buffer, halo_pid
#else
                       write(fu3) halo_output_buffer
#endif        

#ifdef debug_halo_intense
                       write(*,*) 'Replicated halo',frac, halo_output_buffer(1:2)
#endif
                       n_replicate = n_replicate + 1
                       nh_final_local = nh_final_local + 1
                    endif

                    ! 2- 
                    ! Replicate along x+y+ axes (need only to apply y+ now, x+
                    ! replica was added above already):

                    halo_output_buffer(2) = halo_output_buffer(2) + (nc/4)*npc/mag 
                    halo_output_buffer(5) = halo_output_buffer(5) + (nc/4)*npc/mag 

                    !if(halo_output_buffer(1) .le. npc +15. .and. halo_output_buffer(2) .le. npc + 15.) then
                   if((halo_output_buffer(1) < npc +15.) .and. (halo_output_buffer(1).ge.-15.) .and. &
                   (halo_output_buffer(2) < npc + 15.) .and. (halo_output_buffer(2).ge.-15.))then
#ifdef HALOPID
                       write(fu3) halo_output_buffer, halo_pid
#else
                       write(fu3) halo_output_buffer
#endif        

#ifdef debug_halo_intense
                       write(*,*) 'Replicated halo',frac, halo_output_buffer(1:2)
#endif
                       n_replicate = n_replicate + 1
                       nh_final_local = nh_final_local + 1
                    endif

                    
                    ! 3- 
                    ! Replicate along y+ axis (need only to remove x+ do go from x+y+

                    halo_output_buffer(1) = halo_output_buffer(1) - (nc/4)*npc/mag 
                    halo_output_buffer(4) = halo_output_buffer(4) - (nc/4)*npc/mag

                    !if(halo_output_buffer(1) .le. npc + 15. .and. halo_output_buffer(2) .le. npc + 15.) then
                   if((halo_output_buffer(1) < npc +15.) .and. (halo_output_buffer(1).ge.-15.) .and. &
                   (halo_output_buffer(2) < npc + 15.) .and. (halo_output_buffer(2).ge.-15.))then
#ifdef HALOPID
                       write(fu3) halo_output_buffer, halo_pid
#else
                       write(fu3) halo_output_buffer
#endif        

#ifdef debug_halo_intense
                       write(*,*) 'Replicated halo',frac, halo_output_buffer(1:2)
#endif
                       n_replicate = n_replicate + 1
                       nh_final_local = nh_final_local + 1
                    endif


                    ! 4- 
                    ! Replicate along y+x- axis (need to remove x+ now, since  y+
                    ! replicas were added above already):

                    halo_output_buffer(1) = halo_output_buffer(1) - (nc/4)*npc/mag 
                    halo_output_buffer(4) = halo_output_buffer(4) - (nc/4)*npc/mag
                    !halo_output_buffer(2) = halo_output_buffer(2) - (nc/4)*npc/mag 
                    !halo_output_buffer(5) = halo_output_buffer(5) - (nc/4)*npc/mag 

                    !if(halo_output_buffer(1) .le. npc + 15. .and. halo_output_buffer(2) .le. npc + 15.) then
                   if((halo_output_buffer(1) < npc +15.) .and. (halo_output_buffer(1).ge.-15.) .and. &
                   (halo_output_buffer(2) < npc + 15.) .and. (halo_output_buffer(2).ge.-15.))then
#ifdef HALOPID
                       write(fu3) halo_output_buffer, halo_pid
#else
                       write(fu3) halo_output_buffer
#endif        

#ifdef debug_halo_intense
                       write(*,*) 'Replicated halo',frac, halo_output_buffer(1:2)
#endif
                       n_replicate = n_replicate + 1
                       nh_final_local = nh_final_local + 1
                    endif

                    ! 5- 
                    ! Replicate along x- axis (need to remove y+ now, since
                    ! y+x- replicas were added above already):

                    halo_output_buffer(2) = halo_output_buffer(2) - (nc/4)*npc/mag 
                    halo_output_buffer(5) = halo_output_buffer(5) - (nc/4)*npc/mag 

                    !if(halo_output_buffer(1) .le. npc + 15. .and. halo_output_buffer(2) .le. npc + 15.) then
                   if((halo_output_buffer(1) < npc +15.) .and. (halo_output_buffer(1).ge.-15.) .and. &
                   (halo_output_buffer(2) < npc + 15.) .and. (halo_output_buffer(2).ge.-15.))then
#ifdef HALOPID
                       write(fu3) halo_output_buffer, halo_pid
#else
                       write(fu3) halo_output_buffer
#endif        

#ifdef debug_halo_intense
                       write(*,*) 'Replicated halo',frac, halo_output_buffer(1:2)
#endif
                       n_replicate = n_replicate + 1
                       nh_final_local = nh_final_local + 1
                    endif

                    ! 6- 
                    ! Replicate along x-y- axis (need only to remove y- now, & x-
                    ! replicas were added above already):

                    halo_output_buffer(2) = halo_output_buffer(2) - (nc/4)*npc/mag 
                    halo_output_buffer(5) = halo_output_buffer(5) - (nc/4)*npc/mag

                    !if(halo_output_buffer(1) .le. npc + 15. .and. halo_output_buffer(2) .le. npc + 15.) then
                   if((halo_output_buffer(1) < npc +15.) .and. (halo_output_buffer(1).ge.-15.) .and. &
                   (halo_output_buffer(2) < npc + 15.) .and. (halo_output_buffer(2).ge.-15.))then
#ifdef HALOPID
                       write(fu3) halo_output_buffer, halo_pid
#else
                       write(fu3) halo_output_buffer
#endif        

#ifdef debug_halo_intense
                       write(*,*) 'Replicated halo',frac, halo_output_buffer(1:2)
#endif
                       n_replicate = n_replicate + 1
                       nh_final_local = nh_final_local + 1
                    endif

                    ! 7-
                    ! Replicate along y- axis (need only to remove x- now (or add x+), x-y-
                    ! replicas were added above already):

                    halo_output_buffer(1) = halo_output_buffer(1) + (nc/4)*npc/mag 
                    halo_output_buffer(4) = halo_output_buffer(4) + (nc/4)*npc/mag

                    !if(halo_output_buffer(1) .le. npc + 15. .and. halo_output_buffer(2) .le. npc + 15.) then
                   if((halo_output_buffer(1) < npc +15.) .and. (halo_output_buffer(1).ge.-15.) .and. &
                   (halo_output_buffer(2) < npc + 15.) .and. (halo_output_buffer(2).ge.-15.))then
#ifdef HALOPID
                       write(fu3) halo_output_buffer, halo_pid
#else
                       write(fu3) halo_output_buffer
#endif        

#ifdef debug_halo_intense
                       write(*,*) 'Replicated halo',frac, halo_output_buffer(1:2)
#endif
                       n_replicate = n_replicate + 1
                       nh_final_local = nh_final_local + 1
                    endif

                    ! 8- 
                    ! Replicate along y-x+ axis (need only to add x+, y-
                    ! replicas were added above already):

                    halo_output_buffer(1) = halo_output_buffer(1) + (nc/4)*npc/mag 
                    halo_output_buffer(4) = halo_output_buffer(4) + (nc/4)*npc/mag

                    !if(halo_output_buffer(1) .le. npc + 15. .and. halo_output_buffer(2) .le. npc + 15.) then
                   if((halo_output_buffer(1) < npc +15.) .and. (halo_output_buffer(1).ge.-15.) .and. &
                   (halo_output_buffer(2) < npc + 15.) .and. (halo_output_buffer(2).ge.-15.))then
#ifdef HALOPID
                       write(fu3) halo_output_buffer, halo_pid
#else
                       write(fu3) halo_output_buffer
#endif        

#ifdef debug_halo_intense
                       write(*,*) 'Replicated halo',frac, halo_output_buffer(1:2)
#endif
                       n_replicate = n_replicate + 1
                       nh_final_local = nh_final_local + 1
                    endif

                 endif
              else                 
                 n_out = n_out+1     
              endif

#ifdef debug_halo_intense
              write(*,*)'peak=', halo_output_buffer(1:3)
              write(*,*)'x_cm=',halo_output_buffer(4:6)
              write(*,*)'v_cm=',halo_output_buffer(7:9)
              write(*,*)'l=',halo_output_buffer(10:12)

              if(ncol.eq.17) then
                 write(*,*)'v_disp=',halo_output_buffer(13)
                 write(*,*)'radius_scale=',halo_output_buffer(14)
                 write(*,*)'halo_mass=',halo_output_buffer(15)
                 write(*,*)'halo_mass_pp=',halo_output_buffer(16)
                 write(*,*)'halo_mass1=',halo_output_buffer(17)
              else
              !28 colomns
                 write(*,*)'v_disp=',halo_output_buffer(13:15)
                 write(*,*)'radius_scale=',halo_output_buffer(16)
                 write(*,*)'halo_mass=',halo_output_buffer(17)
                 write(*,*)'halo_mass_pp=',halo_output_buffer(18)
                 write(*,*)'halo_mass1=',halo_output_buffer(19)
                 write(*,*)'var=',halo_output_buffer(20:22)              
                 write(*,*)'I_ij=', halo_output_buffer(23:28)             
              endif
#ifdef HALOPID
              write(*,*)'PID=',halo_pid!(:,nh_total)   
#endif              

#endif
              !stop
              !pause

113           continue
           enddo
112        close(fu2)
           close(fu3)
           write(*,*)'closed original and new halo files'

           nh_final = nh_final + nh_final_local
           write(*,*) 'nh final cumulative = ', nh_final ,'nh final local = ', nh_final_local
           write(*,*) 'n_rejected  = ',n_behind +  n_out
           write(*,*) 'n_replicate = ', n_replicate
           write(*,*) 'n_bad = ', BadHalo
        !enddo ! end loop over nodes

        !call cic(xv,d,nc,ip)

#endif
!********************
! done HALO_CATALOGUE      
!********************

#ifndef halo_only   

#ifdef cubepm
        !! substract the mean to get the overdensity
        write(*,*) 'substracting mean'

        !write(*,*) '******** Impose mean: Only OK for map sfrom McCarthy!!! ************'
        ! Ian
        !!!!rhomean=input_map(:,:)/654.5 
        !rhomean=sum(real(input_map(:,:),kind=8))/nc/nc 
        !write(*,*) 'Mean before subtraction= ',rhomean
        !input_map(:,:)=input_map(:,:)-rhomean  
        !rhomean=sum(real(input_map(:,:),kind=8))/nc/nc 
        !!!divide by mean rho 3D -- useless in cubep3m since rho 3D = 1.0
        !input_map(:,:) = input_map(:,:)/0.6392
        !write(*,*) 'Mean after subtraction= ',rhomean

        ! TCS
        !if(z_write(j) <= z_switch) then
           !rhomean=sum(real(map_close(:,:,j-nslice_far),kind=8))/nc/nc 
           !map_close(:,:,j-nslice_far)=map_close(:,:,j-nslice_far)-rhomean  

        ! SLICS   
        rhomean=sum(real(input_map(:,:),kind=8))/nc/nc 
        write(*,*) 'Mean before subtraction= ',rhomean
        input_map(:,:)=input_map(:,:)-rhomean  
        rhomean=sum(real(input_map(:,:),kind=8))/nc/nc 
        write(*,*) 'Mean after subtraction= ',rhomean

        ! TCS
        !else
        !   rhomean=sum(real(map_far(:,:,j),kind=8))/nc/nc
        !   map_far(:,:,j)=map_far(:,:,j)-rhomean !nc 
        !endif

#ifdef test_read
        !write (fn,'(f5.3,"test.dat")') z_write(j)
        !open(10,file=Lens_output_path//"zs_"//digitn(int(zs),1)//"_"//digitn(int(box),3)//"Mpc_"//digitn(int(npc),4)//"/"//fn,form='binary')
        !write(10) map1(:,:,j)
        !close(10)
        !write(*,*)'Wrote test Map for z=',z_write(j)
#endif
#endif

!************
! Get shear :
!************       
#ifdef calshear
!        if(z_write(j) <= z_switch) then
           !write(*,*)'before kappa2defl delta mean,min,max:' , sum(map_close(:,:,j-nslice_far)/nc/nc), minval(map_close(:,:,j-nslice_far)),maxval(map_close(:,:,j-nslice_far)) 
           !call kappa_to_defl(map_close(:,:,j-nslice_far),defl,shear,phi,nc) !Using discrete differentiation
           !write(*,*)'after kappa2defl delta mean,min,max:' , sum(map_close(:,:,j-nslice_far)/nc/nc), minval(map_close(:,:,j-nslice_far)),maxval(map_close(:,:,j-nslice_far)) 
           !call kappa2gamma(map_close(:,:,j-nslice_far),defl,shear,phi,nc)  !Using continuous differentiation
           write(*,*)'before kappa2defl delta mean,min,max:' , sum(input_map(:,:)/nc/nc), minval(input_map(:,:)),maxval(input_map(:,:)) 
           call kappa_to_defl(input_map, defl, shear, phi, map_cplx,nc) !Using discrete differentiation
           !call kappa_to_shear(input_map, defl, shear, phi, map_cplx, nc)
           write(*,*)'after kappa2defl delta mean,min,max:' , sum(input_map(:,:)/nc/nc), minval(input_map(:,:)),maxval(input_map(:,:)) 
           !call kappa2gamma(map_close(:,:,j-nslice_far),defl,shear,phi,nc)  !Using continuous differentiation
!        else
!           write(*,*)'before kappa2defl delta mean,min,max:' , sum(map_far(:,:,j)/nc/nc), minval(map_far(:,:,j)),maxval(map_far(:,:,j)) 
!           call kappa_to_defl(map_far(:,:,j),defl,shear,phi,nc) !Using discrete differentiation
           !call kappa2gamma(map_far(:,:,j),defl,shear,phi,nc)  !Using continuous differentiation
!           write(*,*)'after kappa2defl delta mean,min,max:' , sum(map_far(:,:,j)/nc/nc), minval(map_far(:,:,j)),maxval(map_far(:,:,j)) 
!        endif
        write(*,*) 'Calculated shear and deflection'
        !write(*,*)' Unzoomed shear p mean,min,max:' , sum(shear(:,:)%p/nc/nc), minval(shear(:,:)%p),maxval(shear(:,:)%p) 
        !write(*,*)' Unzoomed shear c mean,min,max:' , sum(shear(:,:)%c/nc/nc), minval(shear(:,:)%c),maxval(shear(:,:)%c) 

        write(*,*)' Unzoomed defl x mean,min,max:' , sum(defl(:,:)%x/nc/nc), minval(defl(:,:)%x),maxval(defl(:,:)%x) 
#endif

        !******************************
        ! Randomly shift the origin and  
        ! zoom the slices by 'frac = 
        !(theta_source/theta_lens)' to fit the field of view
        !******************************
        
        if(chi(z_write(j),omegam,h) .eq.0) then 
           frac=angle/(box / 0.00001) 
        else
           frac=angle/(box / chi(z_write(j),omegam,h) / h) 
        endif
        write(*,*) 'theta_source / theta_lens =', frac
        write(*,*) 'Calling Zoomshiftmap'


        !if(z_write(j) <= z_switch) then
#ifdef calshear
        call zoomshiftmap(input_map,map_3D(:,:,j),zoom_map,nc,npc,shift,newdefl,frac)
#else
        call zoomshiftmap_nodefl(input_map,map_3D(:,:,j),zoom_map,nc,npc,shift,frac)
#endif
        !else
        !   call zoomshiftmap(map_far(:,:,j),map_3D(:,:,j),nc,npc,shift,newdefl,frac)
        !endif
        write(*,*)'delta mean,min,max:' , sum(map_3D(:,:,j)/nc/nc), minval(map_3D(:,:,j)),maxval(map_3D(:,:,j)) 

!##############
#ifdef calshear
#ifndef  SemiBorn


          write(*,*)'shear p mean,min,max:' , sum(shear(:,:)%p/nc/nc), minval(shear(:,:)%p),maxval(shear(:,:)%p) 
          write(*,*)'shear c mean,min,max:' , sum(shear(:,:)%c/nc/nc), minval(shear(:,:)%c),maxval(shear(:,:)%c) 

          tmp_map = shear%p
          call zoomshiftmap(tmp_map,newshear(:,:,j)%p,zoom_map, nc,npc,shift,CorrBornDefl,frac)

          tmp_map = shear%c
          call zoomshiftmap(tmp_map,newshear(:,:,j)%c,zoom_map, nc,npc,shift,CorrBornDefl,frac)
          !call zoomshiftmap(shear%p,newshear(:,:,j)%p,zoom_map, nc,npc,shift,CorrBornDefl,frac)
          !call zoomshiftmap(shear%c,newshear(:,:,j)%c,zoom_map, nc,npc,shift,CorrBornDefl,frac)

          tmp_map = defl%x
          call zoomshiftmap(tmp_map,tmp_pix_map,zoom_map, nc,npc,shift,CorrBornDefl,frac)
          newdefl(:,:,j)%x = tmp_pix_map*(box/nc)/chi(z_write(j),omegam,h)/h
          write(*,*)'newdefl x mean,min,max:' , sum(newdefl(:,:,j)%x/nc/nc), minval(newdefl(:,:,j)%x),maxval(newdefl(:,:,j)%x) 

          tmp_map = defl%y
          call zoomshiftmap(tmp_map,tmp_pix_map,zoom_map, nc,npc,shift,CorrBornDefl,frac)
          newdefl(:,:,j)%y = tmp_pix_map*(box/nc)/chi(z_write(j),omegam,h)/h
 
          !call zoomshiftmap(defl%x,newdefl(:,:,j)%x,zoom_map, nc,npc,shift,CorrBornDefl,frac)
          !call zoomshiftmap(defl%y,newdefl(:,:,j)%y,zoom_map, nc,npc,shift,CorrBornDefl,frac)
          
          ! Scale the angles
          !newdefl(:,:,j)%x=newdefl(:,:,j)%x*(box/nc)/chi(z_write(j),omegam,h)/h
          !newdefl(:,:,j)%y=newdefl(:,:,j)%y*(box/nc)/chi(z_write(j),omegam,h)/h
         
          write(*,*)'newshear p mean,min,max:' , sum(newshear(:,:,j)%p/nc/nc), minval(newshear(:,:,j)%p),maxval(newshear(:,:,j)%p) 
          write(*,*)'newshear c mean,min,max:' , sum(newshear(:,:,j)%c/nc/nc), minval(newshear(:,:,j)%c),maxval(newshear(:,:,j)%c) 

          !CorrBornDefl(:,:)%x=newdefl(:,:,1)%x
          !CorrBornDefl(:,:)%y=newdefl(:,:,1)%y

#endif
#endif
          write(*,*) 'Randomly Shifted and Zoomed on the  maps'        
#endif
! ifndef halo_only


        icount=icount+1
!        if(mod(j,2) .eq. 0) icount=icount+1
!        write(*,*) '******TEST: incrementing rotation every two steps...'
     enddo !!j=1,nslice


#ifdef halo_only      
  enddo !ir=1,nr  
#else

     !*****************
     !loop over sources
     !*****************
     !do i=nslice,18,-1
     do i=nslice,1,-1

       write(*,*) '******************************************'
       write(*,*) 'Processing sources at z =  ', z_write_s(i) 

!############
! Get Kernel:
!############
!! Full Kernel=3/(2)*(H0/c)**2*omegam*chi(z(i))*(1+z(i))*(1-chi(z(i))/chi(zs(j)))*dchi 
!! Partial Kernel=3/(2)*(H0/c)**2*omegam*chi(z(i))*(1+z)*dchi
!! Note the unit of chi.f90 is Mpc, not h^-1*Mpc 
!! H0/c=100/3*e5=3.E3 (h Mpc^-1)
!! Also,  dchi = (box/nc)  
!! And the deflection angle must further be scaled by dChi/Chi_lens 

!###################
#ifdef full_geometry
!###################

        ! Multiply by weight and geometry, integrate up to the lens redshift
        !write(*,*) 'Integrating up to z=',z_write(i)

        cumul_kappa(:,:)=0.0
#ifdef calshear
        cumul_gamma1(:,:)=0.0
        cumul_gamma2(:,:)=0.0
        cumul_deflx(:,:)=0.0
        cumul_defly(:,:)=0.0
#endif

        !****************
        !loop over lenses
        !****************
        do j = nslice,i,-1
 

           ! Get the good box size: Those of the lenses we sum over!
           if(z_write(j) > z_switch) box=sbox
           if(z_write(j) <= z_switch) box=lbox


           !******************************
           ! Randomly shift the origin and  
           ! zoom the slices by 'frac = 
           !(theta_source/theta_lens)' to fit the field of view
           !******************************

           if(chi(z_write(j),omegam,h) .eq.0) then 
              frac=angle/(box / 0.00001) 
           else
              frac=angle/(box / chi(z_write(j),omegam,h) / h) 
           endif

           !write(*,*) 'z = ', z_write(j)
           !write(*,*) 'opening angle = ', angle
           !write(*,*)  'lens angle = ', box/(chi(z_write(j),omegam,h)*h)
           !write(*,*) 'wrong lens angle = ' ,box / chi(z_write(j),omegam,h)*h
           !write(*,*) 'wrong frac = ', frac
           !write(*,*) 'good frac = ' , angle/(box / chi(z_write(j),omegam,h)/h)
           !pause

           !Get the Kernel in!!!
           kernel = (3./2)*omegam*(box/nc)/(3.E3)**2*(chi(z_write(j),omegam,h)*h)*(1+z_write(j))*(1-chi(z_write(j),omegam,h)/chi(z_write_s(i),omegam,h))

           
           cumul_kappa(:,:) =cumul_kappa(:,:) +map_3D(:,:,j)*kernel


!##############
#ifdef calshear

#ifdef SemiBorn
           !In Post Born Approx, kappa is still calculated at pixel location.
           !Only gamma's and defl's are calulated at light ray location.

           write(*,*)'shear p mean,min,max:' , sum(shear(:,:)%p/nc/nc), minval(shear(:,:)%p),maxval(shear(:,:)%p) 
           write(*,*)'shear c mean,min,max:' , sum(shear(:,:)%c/nc/nc), minval(shear(:,:)%c),maxval(shear(:,:)%c) 

           call zoomshiftmap(shear%p,newshear(:,:,j)%p,zoom_map,nc,npc,shift,CorrBornDefl,frac)
           call zoomshiftmap(shear%c,newshear(:,:,j)%c,zoom_map,nc,npc,shift,CorrBornDefl,frac)
           call zoomshiftmap(defl%x,newdefl(:,:,j)%x,zoom_map,nc,npc,shift,CorrBornDefl,frac)
           call zoomshiftmap(defl%y,newdefl(:,:,j)%y,zoom_map,nc,npc,shift,CorrBornDefl,frac)
           
           ! Scale the angles
           newdefl(:,:,j)%x=newdefl(:,:,j)%x*(box/nc)/chi(z_write(j),omegam,h)/h
           newdefl(:,:,j)%y=newdefl(:,:,j)%y*(box/nc)/chi(z_write(j),omegam,h)/h
         
           write(*,*)'newshear p mean,min,max:' , sum(newshear(:,:,j)%p/nc/nc), minval(newshear(:,:,j)%p),maxval(newshear(:,:,j)%p) 
           write(*,*)'newshear c mean,min,max:' , sum(newshear(:,:,j)%c/nc/nc), minval(newshear(:,:,j)%c),maxval(newshear(:,:,j)%c)

           CorrBornDefl(:,:)%x=newdefl(:,:,j)%x*kernel
           CorrBornDefl(:,:)%y=newdefl(:,:,j)%y*kernel

           write(*,*) 'Randomly Shifted and Zoomed on the  maps'
#endif

           cumul_gamma1(:,:)=cumul_gamma1(:,:)+newshear(:,:,j)%p*kernel
           cumul_gamma2(:,:)=cumul_gamma2(:,:)+newshear(:,:,j)%c*kernel
           cumul_deflx(:,:)=cumul_deflx(:,:)+newdefl(:,:,j)%x*kernel
           cumul_defly(:,:)=cumul_defly(:,:)+newdefl(:,:,j)%y*kernel
#endif
           
          
       enddo ! End loop over lenses

!#####
#endif
!#####
! ended ifdef full geometry...



!##############
! Write to file
!##############

!##############
#ifdef z_slices
!##############

#ifdef full_geometry
!-----------------
#ifdef write_kappa
        write (fn,'(f5.3,"kappa_weight.dat")') z_write(i)
        !open(10,file=Lens_output_path//"run_"//digitn(ir,2)//"/"//fn,form='binary')
        !open(10,file=Lens_output_path//trim(LOS_str)//"/"//fn, access = output_access, form = output_form, status = 'replace')
        open(10,file=Lens_output_path//trim(fn)//'_LOS'//trim(LOS_str), access = output_access, form = output_form, status = 'replace')
        write(10) cumul_kappa(:,:)
        close(10)
        write(*,*)'Wrote kappa Map for z=',z_write(i)
#endif
!----------------
#ifdef calshear
        write (fn,'(f5.3,"gamma1_weight.dat")') z_write(i)
        !open(11,file=Lens_output_path//"run_"//digitn(ir,2)//"/"//fn,form='binary')
        !open(11,file=Lens_output_path//trim(LOS_str)//"/"//fn, access = output_access, form = output_form, status = 'replace')
        open(11,file=Lens_output_path//trim(fn)//'_LOS'//trim(LOS_str), access = output_access, form = output_form, status = 'replace')
        write(11) cumul_gamma1(:,:)
        close(11)

        write (fn,'(f5.3,"gamma2_weight.dat")') z_write(i)
        !open(10,file=Lens_output_path//"run_"//digitn(ir,2)//"/"//fn,form='binary')\
        !open(12,file=Lens_output_path//trim(LOS_str)//"/"//fn, access = output_access, form = output_form, status = 'replace')
        open(12,file=Lens_output_path//trim(fn)//'_LOS'//trim(LOS_str), access = output_access, form = output_form, status = 'replace')
        write(12) cumul_gamma2(:,:)
        close(12)
        write(*,*)'Wrote shear Maps for z=',z_write(i)
#endif
#ifdef write_defl
        write (fn,'(f5.3,"deflx_weight.dat")') z_write(i)
        !open(10,file=Lens_output_path//"run_"//digitn(ir,2)//"/"//fn,form='binary')
        !open(13,file=Lens_output_path//trim(LOS_str)//"/"//fn, access = output_access, form = output_form, status = 'replace')
        open(13,file=Lens_output_path//trim(fn)//'_LOS'//trim(LOS_str), access = output_access, form = output_form, status = 'replace')
        write(13) cumul_deflx(:,:)
        close(13)
        write (fn,'(f5.3,"defly_weight.dat")') z_write(i)
        !open(10,file=Lens_output_path//"run_"//digitn(ir,2)//"/"//fn,form='binary')\
        !open(14,file=Lens_output_path//trim(LOS_str)//"/"//fn, access = output_access, form = output_form, status = 'replace')
        open(14,file=Lens_output_path//trim(fn)//'_LOS'//trim(LOS_str), access = output_access, form = output_form, status = 'replace')
        write(14) cumul_defly(:,:)
        close(14)
        write(*,*)'Wrote shear Maps for z=',z_write(i)
#endif
!----------------
#ifdef power_spectrum

     power = 0.0

     !-------------
     ! No zero-pad:
     !-------------

     call ps2_deconvolve(cumul_kappa,cumul_kappa,power,npc)
     !!call ps2(cumul_gamma1,cumul_gamma1,power,npc)
     !!call ps2(cumul_gamma2,cumul_gamma2,power,npc)
     !write (fn,'(f5.3,"l2cl_kappa_ngp_new.dat_LOS")') z_write(i)
     write (fn,'(f5.3,"l2cl_kappa_ngp.dat_LOS")') z_write(i)
     !write (fn,'(f5.3,"l2cl_kappa.dat_LOS")') z_write(i)
     !!write (fn,'(f5.3,"l2cl_gamma1.dat_LOS")') z_write(i)
     !!write (fn,'(f5.3,"l2cl_gamma2.dat_LOS")') z_write(i)

    
     !----------
     ! Zero-pad:
     !----------
     !zero-pad by inserting cumul_kappa at the center of the zeropad_map array, which is (2*npc, 2*npc):
     ! So the array starts at npc/2 + 1 up to 3*npc/2 in both dimensions.
     
     !zeropad_map = 0
     !zeropad_map(npc/2+1:3*npc/2, npc/2+1:3*npc/2) = cumul_kappa
     !call ps2_deconvolve(zeropad_map,zeropad_map,power,npc*2)
     !write (fn,'(f5.3,"l2cl_kappa.dat_zeropad_LOS")') z_write(i)

     !--- Then, for both: 


     open(30,file=Lens_output_path//trim(fn)//trim(LOS_str))
     do i1=1,npc/2   ! leave the npc/2 for no zeropad, remove the '/2' for zeropad
        write(30,*) i1*(2*pi/angle - 1),power(i1)
        !write(30,*) i*2*pi/angle,power(i)
     enddo
     close(30)
     write(*,*) 'Got l2Cl for z=',z_write(i)
#endif
!-----------------
#endif 
!End #ifdef full_geometry

!----------------------
#ifdef write_delta_maps
        write(*,*) 'Wrinting delta maps'
        write (fn,'(f5.3,"delta.dat_bicubic")') z_write(i)
        !open(10,file=Lens_output_path//"run_"//digitn(ir,2)//"/"//fn,form='binary')
        !open(15,file=Lens_output_path//trim(LOS_str)//"/"//fn, access = output_access, form = output_form, status = 'replace')
        open(15,file=Lens_output_path//trim(fn)//'_LOS'//trim(LOS_str), access = output_access, form = output_form, status = 'replace')
        write(15) map_3D(:,:,i)
        close(15)
        write(*,*)'Wrote delta Map for z=',z_write(i)
#endif
!-------------------

!#####
#endif
!##### 
!end of '#ifdef z_slice'
 

     enddo ! End loop over sources

     write(*,*) 'Done integral over redshift slices'


!###############   
!#ifdef integrate      
!#ifdef full_geometry
     !open(20,file=Lens_output_path//fn_kappa,form='binary')
     !open(20,file=Lens_output_path//trim(LOS_str)//"/"//fn_kappa, access = output_access, form = output_form, status = 'replace')
!     open(20,file=Lens_output_path//trim(fn_kappa)//'_LOS'//trim(LOS_str), access = output_access, form = output_form, status = 'replace')
!     write(20)cumul_kappa(:,:)
!     close(20)
!     write(*,*)'Wrote integrated kappa Maps'

!#ifdef calshear
     !open(20,file=Lens_output_path//trim(LOS_str)//"/"//fn_gamma1, access = output_access, form = output_form, status = 'replace')
!     open(20,file=Lens_output_path//trim(fn_gamma1)//'_LOS'//trim(LOS_str), access = output_access, form = output_form, status = 'replace')
!     write(20) cumul_gamma1(:,:)
!     close(20)
     !finalmap=0
     !do i=1,nslice
     !   if(z_write(i).le.zs)finalmap(:,:)=finalmap(:,:)+newshear(:,:,i)%c
     !enddo
     !open(20,file=Lens_output_path//fn_gamma2,form='binary')
     !open(20,file=Lens_output_path//trim(LOS_str)//"/"//fn_gamma2, access = output_access, form = output_form, status = 'replace')
!     open(20,file=Lens_output_path//trim(fn_gamma2)//'_LOS'//trim(LOS_str), access = output_access, form = output_form, status = 'replace')
!     write(20) cumul_gamma2(:,:)
!     close(20)     
!     write(*,*)'Wrote integrated shear Maps'
!#endif

!#ifdef write_defl
     !open(20,file=Lens_output_path//fn_defl,form='binary')
     !open(20,file=Lens_output_path//trim(LOS_str)//"/"//fn_defl, access = output_access, form = output_form, status = 'replace')
!     open(20,file=Lens_output_path//frim(fn_defl)//'_LOS'//trim(LOS_str), access = output_access, form = output_form, status = 'replace')
!     write(20)cumul_deflx(:,:),cumul_defly(:,:)
!     close(20)

!     write(*,*)'Wrote integrated defl Maps'
!#endif

!####################
!#ifdef power_spectrum
!-----------------
!#ifdef write_kappa 
!     power = 0.0
!     call ps2(cumul_kappa,cumul_kappa,power,npc)
     !open(30,file=Lens_output_path//'l2cl_kappa.'//digitn(ir,2))
!     open(30,file=Lens_output_path//trim(LOS_str)//'/l2cl_kappa.'//digitn(ir,2))
!     do i=1,npc/2
!        write(30,*) i*2*pi/angle,power(i)
!     enddo
!     close(30)
!#endif
!-----------------
!#ifdef calshear
!     power = 0.0
!     call ps2(cumul_gamma1,cumul_gamma1,power,npc)
     !open(30,file=Lens_output_path//'l2cl_kappa.'//digitn(ir,2))
!     open(30,file=Lens_output_path//trim(LOS_str)//'/l2cl_gamma1.'//digitn(ir,2))
!     do i=1,npc/2
!        write(30,*) i*2*pi/angle,power(i)
!     enddo
!     close(30)
!     power = 0.0
!     call ps2(cumul_gamma2,cumul_gamma2,power,npc)
     !open(30,file=Lens_output_path//'l2cl_kappa.'//digitn(ir,2))
!     open(30,file=Lens_output_path//trim(LOS_str)//'/l2cl_gamma2.'//digitn(ir,2))
!     do i=1,npc/2
!        write(30,*) i*2*pi/angle,power(i)
!     enddo
!     close(30)
!#endif
!-----------------
!     write(*,*)'Wrote Power Spectrum'
!#endif  
!####################   
!#else
!     write(*,*) 'Cannot integrate without full geometry'
!#endif
!end of full geometry
!#endif
!end of integrate

  enddo !! ir=1,nr


#ifdef mix_nbody_runs
  open(150,file=Lens_output_path//trim(LOS_str)//'/count_run.dat', access = 'sequential')
!  do i=1,num_nbody_runs
!     write(150,*) i,count_run(i)
!  enddo
  do i=1,nslice
     write(150,*) rand_LOS(i)
  enddo
  close(150)
#endif

#endif 
!ifndef halo_only

end program SimulLens

!--------------------

  subroutine zoomshiftmap(map1,map2,map3, nc,npc,shift,defl,frac)
    use Lensing
    implicit none
    
    integer nc,npc 
    real,dimension(nc,nc) :: map1
    real,dimension(2*nc,2*nc) :: map3 
    type(vec2D) shift
    type(vec2D), dimension(npc,npc) :: defl
    real frac

    real, dimension(npc,npc) :: map2 
    integer i1,j1
   
    !write(*,*) 'Inside zoomshiftmap'
    map3 = 0.0

    !$ call omp_set_num_threads(24)
    !$omp parallel do default(shared) private(j1,i1)
    do j1=1,nc
       do i1=1,nc
          map3(i1,j1)=map1(i1,j1)
          map3(i1+nc,j1)=map1(i1,j1)
          map3(i1,j1+nc)=map1(i1,j1)
          map3(i1+nc,j1+nc)=map1(i1,j1)
       enddo
    enddo
    !$omp end parallel do        

    call zoom(map3(:,:),map2(:,:),nc,npc,shift,defl,frac)
    !call zoom_bicubic(map3(:,:),map2(:,:),nc,npc,shift,defl,frac)
  
    !write(*,*) 'Done zoom'
   
    return
  end subroutine zoomshiftmap 

!--------------------

  subroutine zoomshiftmap_nodefl(map1,map2,map3,nc,npc,shift,frac)
    use Lensing
    implicit none
    
    integer nc,npc 
    real,dimension(nc,nc) :: map1
    real,dimension(2*nc,2*nc) :: map3 
    type(vec2D) shift
    real frac

    real, dimension(npc,npc) :: map2 
    integer i1,j1
   
    map3 = 0.0

    !write(*,*) 'Inside zoomshiftmap'

    !$ call omp_set_num_threads(24)
    !$omp parallel do default(shared) private(j1,i1)
    do j1=1,nc
       do i1=1,nc
          map3(i1,j1)=map1(i1,j1)
          map3(i1+nc,j1)=map1(i1,j1)
          map3(i1,j1+nc)=map1(i1,j1)
          map3(i1+nc,j1+nc)=map1(i1,j1)
       enddo
    enddo
    !$omp end parallel do        
        
    call zoom_nodefl(map3(:,:),map2(:,:),nc,npc,shift,frac)
    !call zoom_bicubic(map3(:,:),map2(:,:),nc,npc,shift,defl,frac)
  
    !write(*,*) 'Done zoom_nodefl'
   
    return
  end subroutine zoomshiftmap_nodefl

!-------------

  subroutine zoom(tmp,map2,n1,n2,shift,defl,frac)
  !! tmp doubles the original slice by the periodic condition
    use Lensing 
    implicit none

    integer n1,n2,i,j,i1,j1,i2,j2,ib,ip,jp,jb
    type(vec2D) shift
    type(vec2D), dimension(n2,n2) :: defl 
    real map2(n2,n2),x,y,frac,tmp(2*n1,2*n1),w,frac1,s,w1,w2

    !write(*,*) 'Inside zoom subroutine'

    map2=0
    frac1=n1*frac/n2

    !write(*,*) '***************'
    !write(*,*) 'frac = ',frac,frac1 
    !write(*,*) '***************'
    !pause

    !$ call omp_set_num_threads(24)
    !$omp parallel do default(shared) private(j2,i2,jb,ib,w2,w1,ip,jp)
    do j2=1,n2
       do i2=1,n2
          jb=shift%y*n1+(j2-1)*frac1+1!+defl(i2,j2)%y*n2
          ib=shift%x*n1+(i2-1)*frac1+1!+defl(i2,j2)%x*n2
          w2=shift%y*n1+(j2-1)*frac1+1-jb!+defl(i2,j2)%y*n2
          w1=shift%x*n1+(i2-1)*frac1+1-ib!+defl(i2,j2)%x*n2
          jb=modulo(jb-1,2*n1)+1
          ib=modulo(ib-1,2*n1)+1
          jp=modulo(jb,2*n1)+1
          ip=modulo(ib,2*n1)+1
          map2(i2,j2)=tmp(ip,jp)*w1*w2+tmp(ip,jb)*w1*(1-w2)+&
                          tmp(ib,jp)*(1-w1)*w2+tmp(ib,jb)*(1-w1)*(1-w2)
       enddo
    enddo
    !$omp end parallel do        

    return
  end subroutine zoom

!-------------

  subroutine zoom_nodefl(tmp,map2,n1,n2,shift,frac)
  !! tmp doubles the original slice by the periodic condition
    use Lensing 
    implicit none

    integer n1,n2,i,j,i1,j1,i2,j2,ib,ip,jp,jb
    type(vec2D) shift
    real map2(n2,n2),x,y,frac,tmp(2*n1,2*n1),w,frac1,s,w1,w2

    !write(*,*) 'Inside zoom subroutine'

    map2=0
    frac1=n1*frac/n2

    !write(*,*) '***************'
    !write(*,*) 'frac = ',frac,frac1 
    !write(*,*) '***************'
    !pause

    !$ call omp_set_num_threads(24)
    !$omp parallel do default(shared) private(j2,i2,jb,ib,w2,w1,ip,jp)
    do j2=1,n2
       do i2=1,n2
          jb=shift%y*n1+(j2-1)*frac1+1
          ib=shift%x*n1+(i2-1)*frac1+1
          w2=shift%y*n1+(j2-1)*frac1+1-jb
          w1=shift%x*n1+(i2-1)*frac1+1-ib
          jb=modulo(jb-1,2*n1)+1
          ib=modulo(ib-1,2*n1)+1
          jp=modulo(jb,2*n1)+1
          ip=modulo(ib,2*n1)+1
          map2(i2,j2)=tmp(ip,jp)*w1*w2+tmp(ip,jb)*w1*(1-w2)+&
                          tmp(ib,jp)*(1-w1)*w2+tmp(ib,jb)*(1-w1)*(1-w2)
       enddo
    enddo
    !$omp end parallel do        

    return
  end subroutine zoom_nodefl

  !***********************
  subroutine zoom_bicubic(tmp,map2,n1,n2,shift,defl,frac)
  !! tmp doubles the original slice by the periodic condition
    use Lensing 
    use omp_lib
    implicit none

    integer n1,n2,i,j,i1,j1,i2,j2,ib,ip,jp,jb,nthreads
    type(vec2D) shift
    type(vec2D), dimension(n2,n2) :: defl 
    real map2(n2,n2),x,y,frac,tmp(2*n1,2*n1),w,frac1,s,w1,w2
    real ddx_map1(2*n1,2*n1), ddy_map1(2*n1,2*n1), ddxy_map1(2*n1,2*n1)
    real A(16,16), coeff(4,4)
    real X_array(16), alpha(16)

    write(*,*)'Using bicubic interpolation'
    !write(*,*)'defl mean,min,max:' , sum(defl%x/n2/n2), minval(defl%x),maxval(defl%x) 
    !write(*,*)'defl mean,min,max:' , sum(defl%y/n2/n2), minval(defl%y),maxval(defl%y) 

    map2=0
    ddx_map1  = 0
    ddy_map1  = 0
    ddxy_map1 = 0


    frac1=n1*frac/n2

    !write(*,*) '***************'
    !write(*,*) 'frac = ',frac 
    !write(*,*) '***************'
    !pause

    !OMP PARALLEL?!?!?

#ifdef OPENMP

    !call OMP_SET_NUM_THREADS(8)
    !nthreads =  omp_get_num_threads()
    nthreads =  1
    write(*,*)'Num OMP threads = ', nthreads
#endif    

    !1- Compute d/dx, d/dy, d/dxdy on the grid

    !$omp parallel do default(shared) private(i1,j1)
    do j1 = 1,2*n1
       do i1 = 1,2*n1
          ddx_map1(i1,j1) = (tmp(modulo(i1,2*n1)+1,j1)-tmp(modulo(i1-2,2*n1)+1,j1))/2.0
          ddy_map1(i1,j1) = (tmp(i1,modulo(j1,2*n1)+1)-tmp(i1,modulo(j1-2,2*n1)+1))/2.0
          ddxy_map1(i1,j1) =(tmp(modulo(i1,2*n1)+1,modulo(j1,2*n1)+1)- &
               tmp(modulo(i1-2,2*n1)+1,modulo(j1-2,2*n1)+1))/sqrt(8.0)

          !write(*,*) (modulo(i1,2*n1)+1, modulo(i1-2,2*n1) +1)/2
          !pause
       enddo
    enddo
    !$omp end parallel do

    !2- load the bicubic matrix

    A(1,:) =  (/  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
    A(2,:) =  (/  0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
    A(3,:) =  (/ -3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
    A(4,:) =  (/  2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
    A(5,:) =  (/  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 /)
    A(6,:) =  (/  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 /)
    A(7,:) =  (/  0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0 /)
    A(8,:) =  (/  0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0 /)
    A(9,:) =  (/ -3, 0, 3, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0 /)
    A(10,:) = (/  0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0,-2, 0,-1, 0 /)
    A(11,:) = (/  9,-9,-9, 9, 6, 3,-6,-3, 6,-6, 3,-3, 4, 2, 2, 1 /)
    A(12,:) = (/ -6, 6, 6,-6,-3,-3, 3, 3,-4, 4,-2, 2,-2,-2,-1,-1 /)
    A(13,:) = (/  2, 0,-2, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0 /)
    A(14,:) = (/  0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 1, 0, 1, 0 /)
    A(15,:) = (/ -6, 6, 6,-6,-4,-2, 4, 2,-3, 3,-3, 3,-2,-1,-2,-1 /)
    A(16,:) = (/  4,-4,-4, 4, 2, 2,-2,-2, 2,-2, 2,-2, 1, 1, 1, 1 /)

    !3- Loop over pixels, get the 4 pixel corners

    !$omp parallel do default(shared) private(i1,j1,i2,j2,jb,jp,ib,ip,w1,w2,X_array,alpha,coeff)
    do j2=1,n2
       do i2=1,n2
          jb=shift%y*n1+(j2-1)*frac1+1+defl(i2,j2)%y*n2
          ib=shift%x*n1+(i2-1)*frac1+1+defl(i2,j2)%x*n2
          jb=modulo(jb-1,2*n1)+1
          ib=modulo(ib-1,2*n1)+1
          jp=modulo(jb,2*n1)+1
          ip=modulo(ib,2*n1)+1
 
          !write(*,*) ib,ip,jb,jp
          !pause

          ! The interpolation formula is between 0 and 1, so we use the weigths!
          ! (1-w) = distance to lower point
          w2=shift%y*n1+(j2-1)*frac1+1-jb+defl(i2,j2)%y*n2 
          w1=shift%x*n1+(i2-1)*frac1+1-ib+defl(i2,j2)%x*n2 



          !4- load the function and derivatives at the 4 corners
          !   f(x,y), dfdx(x,y), dfdy(x,y), dfdxy(x,y)
          
          X_array(1:4)   =   (/ tmp(ib,jb),tmp(ip,jb),tmp(ib,jp),tmp(ip,jp)/)
          X_array(5:8)   =   (/ ddx_map1(ib,jb),ddx_map1(ip,jb),ddx_map1(ib,jp),ddx_map1(ip,jp)  /)
          X_array(9:12)  =   (/ ddy_map1(ib,jb),ddy_map1(ip,jb),ddy_map1(ib,jp),ddy_map1(ip,jp)  /)
          X_array(13:16) =   (/ ddxy_map1(ib,jb),ddxy_map1(ip,jb),ddxy_map1(ib,jp),ddxy_map1(ip,jp)  /)

          !5- Get the bicubic coefficients alpha 
          !   i1 loop over lines, j1 loops over columns

          alpha = 0

          do i1 = 1,16
             do j1 = 1,16
                alpha(i1) = alpha(i1) + A(i1,j1)*X_array(j1)
            enddo
          enddo

          coeff(1,:)= alpha(1:4)
          coeff(2,:)= alpha(5:8)
          coeff(3,:)= alpha(9:12)
          coeff(4,:)= alpha(13:16)

          !write(*,*) 'coeff = ', coeff

          !6- Get the values!

          do i1 = 1,4
             do j1 = 1,4
                !map2(i2,j2) = map2(i2,j2) + coeff(i1,j1)*(w1)**(i1-1)*(w2)**(j1-1)
                !map2(i2,j2) = map2(i2,j2) + coeff(i1,j1)*(1-w1)**(i1-1)*(1-w2)**(j1-1)
                !map2(i2,j2) = map2(i2,j2) + coeff(i1,j1)*(1-w1)**(j1-1)*(1-w2)**(i1-1) !upside down?
                map2(i2,j2) = map2(i2,j2) + coeff(i1,j1)*(w1)**(j1-1)*(w2)**(i1-1) !upside down?
                !map2(i2,j2) = map2(i2,j2) + coeff(i1,j1)*(1-w1)**(j1-1)*(1-w2)**(i1-1) !upside down?
              enddo
          enddo


#ifdef debug_bicubic          
          write(*,*) 'Starting map(',i2,j2,')'
          write(*,*) 'Nearest neighbours   : ', tmp(ib,jp),tmp(ip,jp)
          write(*,*) '                     : ', tmp(ib,jb),tmp(ip,jb)
          write(*,*) '1-w1 = ' ,(1-w1)
          write(*,*) '1-w2 = ' ,(1-w2)
          write(*,*) 'linear interpolation : ', tmp(ip,jp)*w1*w2+tmp(ip,jb)*w1*(1-w2)+&
                          tmp(ib,jp)*(1-w1)*w2+tmp(ib,jb)*(1-w1)*(1-w2)         
          write(*,*) 'bicubic interpolation : ', map2(i2,j2) 
          pause 
#endif


       enddo
    enddo
    !$omp end parallel do

    return
  end subroutine zoom_bicubic
!-------------

subroutine init_var
  !use Lensing
  implicit none
  include 'Lens.fh'

  integer(4) :: k

!  newshear(:,:,:)%p  =0.
!  newshear(:,:,:)%c  =0.
!  newdefl(:,:,:)%x  =0.
!  newdefl(:,:,:)%y  =0.

!  CorrBornDefl(:,:)%x=0.
!  CorrBornDefl(:,:)%y=0.

!  do k = 1,nslice
!    map_3D(:,:,k) =0.
!  enddo
#ifdef write_delta_maps
  write(*,*) '****************'
  write(*,*) 'write delta maps'
  write(*,*) '****************'
#endif
#ifdef write_kappa
  write(*,*) '****************'
  write(*,*) 'write kappa maps'
  write(*,*) '****************'
#endif
#ifdef full_geometry
  write(*,*) '*************'
  write(*,*) 'full geometry'
  write(*,*) '*************'
#endif
#ifdef calshear
  write(*,*) '*****************'
  write(*,*) 'calculating shear'
  write(*,*) '*****************'
#endif

  write(*,*) '*****************'
  write(*,*) 'Running SimulLens'
  write(*,*) '*****************'

  write(*,*) 'nc_far :',nc
  write(*,*) 'nc_close:',nc
  write(*,*) 'npc:',npc
  write(*,*) 'max source redshift:',zs
  !write(*,*) 'Done init_var'

  return
end subroutine init_var

!--------------
subroutine shift_maps(shift, method, patch,j,newLOS_str)
  use StringOpt
  use Lensing
  implicit none
  include 'Lens.fh'


    type(vec2D) shift
    integer method, patch, j, i3, j2
    character(len=7) newLOS_str

!------------------------
! Randomly shift  the map
!------------------------

!method =  1: Write shift to a file
!method = -1: Read shift from file
!method =  0: No shift

        if(method.eq.1)then
           call random_number(shift%x)
           call random_number(shift%y)

           if(patch==0)then

           !---------------------------
           ! Replace the old shift file
           if(j.eq.1)then
              open(100,file=Lens_output_path//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(npc),4)//"/random_shift")
              write(100,*) shift%x, shift%y
              close(100)
           else
              open(100,file=Lens_output_path//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(npc),4)//"/random_shift", position = 'append')
              write(100,*) shift%x, shift%y
              close(100)
           endif
           !----------------------------

           else  !patch = 1

              !-----------------------------------
              ! must read the other file's random shift and change it first!
              open(101,file='/scratch/jharno/SimulLens_2048_to_1024/LOS'//trim(newLOS_str)//'/zs_'//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(npc),4)//"/random_shift", status='old')
              
              do i3=1,j
                 read(101,*) shift%x, shift%y
              enddo
              close(101)

              if (shift%x .le. 0.5) then
                  shift%x = shift%x + 0.5
              else 
                   shift%x = shift%x - 0.5
              endif
              if (shift%y .le. 0.5) then 
                  shift%y = shift%y + 0.5
              else
                  shift%y = shift%y - 0.5
              endif

              open(100,file=Lens_output_path//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(npc),4)//"/random_shift", position = 'append')
              write(100,*) shift%x, shift%y
              close(100)
              !--------------------------------------
             
           endif


        elseif(method.eq.-1)then
           open(100,file=Lens_output_path//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(npc),4)//"/random_shift", status='old',position='asis')
           ! Sloppy way to get the right entry number
           do j2= 1,j          
              read(100,*) shift%x, shift%y
           enddo
           close(100)
        elseif(method.eq.0)then
           shift%x = 0;
           shift%y = 0;
        else
           write(*,*) 'Must specify a shift type. (No shift,  read from or write to file)'
        endif     
           
        write(*,*)'Random shift ', shift%x, shift%y
!--------------
! Done shifting
!--------------        

return
end subroutine shift_maps

!  subroutine read_halo(peak,xv,l,v_disp,radius_scale,halo_mass,halo_mass_pp, halo_mass1,var,halo_pid,fn)
!    implicit none

    !real(4), dimension(20) :: halo_input_buffer
    !integer nh_total, nploc(nn), ip,ii,fu2,file_status
    !character (len=4) :: rank_s
    !character(len=180) :: halofn,fn
    !integer(kind=8), dimension(10,10000) :: halo_pid 
    !integer(kind=8), parameter :: np_max = hc**3 
    !real, dimension(np_max) :: v_disp, radius_scale, halo_mass, halo_mass_pp, halo_mass1
    !real, dimension(3,np_max) :: peak, l, var
    !real, dimension(6,np_max) :: xv
    


    !ip = 0
    !nh_total = 0
    !do ii=1,nn
       
    !   write(*,*) 'Reading Node ',ii
    !   write(rank_s,'(i4)') ii-1
    !   rank_s=adjustl(rank_s)
    !   !write(zstring,'(f5.3)') z3dps       ! Need (f6.3) for Z > 10
    !   halofn=fn(1:len_trim(fn))//'halo'//rank_s(1:len_trim(rank_s))//".dat"        
    !   
    !   fu2=20+ii
    !   open(unit=fu2,file=halofn,form='binary',status='old',iostat=file_status)
    !   write(*,*) 'Opened', halofn
    !   
    !   read(fu2) nploc(ii)
    !   !write(*,*) 'nploc(ii) =',nploc(ii)
    !   do
    !      nh_total=nh_total+1
    !      !write(*,*)'Reading halo', nh_total
!#ifdef HALOPID
 !         read(fu2,end=112,err=113) halo_input_buffer, halo_pid(:,nh_total)           
!#else
 !         read(fu2,end=112,err=113) halo_input_buffer
!#endif              
    !      peak(:,nh_total)=halo_input_buffer(1:3)
    !      xv(:,nh_total)=halo_input_buffer(4:9)
    !      l(:,nh_total)=halo_input_buffer(10:12)
    !      v_disp(nh_total)=halo_input_buffer(13)
    !      radius_scale(nh_total)=halo_input_buffer(14)
    !      halo_mass(nh_total)=halo_input_buffer(15)
    !      halo_mass_pp(nh_total)=halo_input_buffer(16)
    !      halo_mass1(nh_total)=halo_input_buffer(17)
    !      var(:,nh_total)=halo_input_buffer(18:20)
    !      
          !write(*,*)'peak=',peak(:,nh_total)
          !write(*,*)'xv=',xv(:,nh_total)
          !write(*,*)'l=',cubepm_l(:,nh_total)
          !write(*,*)'radius_scale=',cubepm_radius_scale(nh_total)
          !write(*,*)'halo_mass=',halo_mass(nh_total)
          !write(*,*)'halo_mass_pp=',halo_mass_pp(nh_total)
          !write(*,*)'halo_mass1=',halo_mass1(nh_total)
          !write(*,*)'var=',var(:,nh_total)
!#ifdef HALOPID
          !write(*,*)'PID=',halo_pid(:,nh_total)   
!#endif
          
!113       continue
!       enddo
!112    close(fu2)
       
       
       
!       ip=ip+nploc(ii)
!       write(*,*) 'np cumulative = ', ip,', np local = ', nploc(ii)
       
!    enddo
    
!    return
!  end subroutine read_halo
  
