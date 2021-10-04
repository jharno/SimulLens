! Lensing code written by TingTing Lu, modifid by Joachim Harnois-Deraps
! to accomodate more than one box size, and to compute the geometry with source
! redshifts not equal to lens redshifts -- i.e. sources are behind the box, 
! while lenses are at the mid plane.
!
! Organized to read 2D slabs written by cubep3m.
!
! This code also zoom on the halo catalogue and produces a 'zoomed catalogue',
! for exact matching with lenses.
!
!
! To be run on the TCS or GPC cluster at SciNet

program SimulLens
  use StringOpt
  use Lensing
  implicit none
  include 'Lens.fh'

  integer, dimension(num_nbody_runs) :: sublist_run,count_run
  real, dimension(nc,nc) :: rho_pxy_far!,rho_pxz,rho_pyz
  real, dimension(nc,nc) :: rho_pxy_close
  real, dimension(npc,npc,nslice) :: map
  real, dimension(nc,nc,nslice_far) :: map_far
  real, dimension(nc,nc,nslice_close) :: map_close
  real, dimension(npc,npc) :: cumul_kappa,cumul_gamma1,cumul_gamma2,cumul_deflx,cumul_defly 
  real, dimension(npc) :: power
  real, dimension(nslice) :: z_write,z_write_s
  real, dimension(nc,nc) :: phi !not up to date with two boxes
  type(vec2D), dimension(nc,nc) :: defl !not up to date with two boxes
  type(vec2D), dimension(npc,npc,nslice) :: newdefl
  type(vec2D), dimension(npc,npc) :: CorrBornDefl
  type(Spin2), dimension(nc,nc) :: shear  !not up to date with two boxes
  type(Spin2), dimension(npc,npc,nslice) :: newshear
  

  type(vec2D) shift
  integer i,j,fu,i1,j1,j2,i3,icount,kx,ky,ir,index_nbody_run,file_status , newLOS,patch!,redshift, z
  real lense_weight,frac,ncr,angle,dang,chi,pi,random_index_run,box,kernel ,ScaleFactor
  real(kind=8) rhomean
  character(len=180) :: fn,fn1,fp,my_status, test_str
  character(len=7) z_string, index_str, newLOS_str,LOS_str

  !-------- Multi-plateform code:
  ! in/output_form =  TCS: 'unformatted' ; GCP: 'binary'
  ! output_access  =  TCS: 'stream'      ; GPC: 'sequential'

  ! GCP:
  character(*), parameter :: output_form=   'binary'
  character(*), parameter :: output_access= 'sequential' 
  character(*), parameter :: input_form=    'binary' 

  ! TCS:
  !character(*), parameter :: output_form=   'unformatted'
  !character(*), parameter :: output_access= 'stream' 
  !character(*), parameter :: input_form=    'unformatted' 

  !character(len=max_path) :: seedfile
  integer(4) :: seedsize
  integer(4), allocatable, dimension(:) :: iseed,old
  integer k_seed, clock
  integer, dimension(8) :: values
  real  rand
  integer, dimension(nslice) :: rand_LOS

  real, parameter :: z_switch = 2.51 ! 2.51!1.0 !2.4 = CMBLens ! 1.0 = CFHT Lens, 2.51 = RCSLenS

!--------------------
#ifdef HALO_CATALOGUE
  
!not up to date with two boxes

  real(4), dimension(28) :: halo_input_buffer, halo_output_buffer
  integer nh_total, nploc(nn), ip, fu2,fu3,ii,jj,kk,ll, nh_final, nh_final_local
  !integer(kind=8), dimension(10,10000) :: halo_pid
  integer(kind=8),dimension(10) :: halo_pid
  character (len=4) :: rank_s
  character (len=180) :: halofn, new_halofn
  integer(kind=8), parameter :: np_max = hc**3 
  !real, dimension(np_max) :: v_disp, radius_scale, halo_mass, halo_mass_pp, halo_mass1
  real, dimension(3,np_max) :: var, l!, peak
  real, dimension(6,np_max) :: xv
  real, dimension(nc+2,nc,nc) :: d
  real, dimension(3) :: euler, I_ij_vector_in, I_ij_vector_out
  real :: mag, chi_halo, dummy1, dummy2  

#endif
!--------------------

! Select one of these
  integer, parameter ::  write_shift = 1
  integer, parameter ::  read_shift = 0
  integer, parameter ::  no_shift = 0
  integer, parameter ::  method = -1
!method =  1: Write shift to a file
!method = -1: Read shift from file
!method =  0: No shift



  box=lbox
  ncr=nc
  pi=acos(-1.)

  map(:,:,:)  =0;
  newshear(:,:,:)%p  =0.
  newshear(:,:,:)%c  =0.
  newdefl(:,:,:)%x  =0.
  newdefl(:,:,:)%y  =0.

  CorrBornDefl(:,:)%x=0.
  CorrBornDefl(:,:)%y=0.

  write(*,*) '*****************'
  write(*,*) 'Running SimulLens'
  write(*,*) '*****************'

  write(*,*) 'nc_far :',nc
  write(*,*) 'nc_close:',nc
  write(*,*) 'npc:',npc
  write(*,*) 'source redshift:',zs
 
 !Field of view and resolution:
  !angle=5./180.*pi  ! For CMBLens
  !angle = sqrt(12.84)/180*pi  ! For UBC Lens
  angle = sqrt(Area)/180*pi  ! For UBC Lens
  !angle=box/chi(zs,omegam,h)/h
  dang=angle/npc

  write(*,*) 'Field of view :',(angle*180./pi)**2,'sq.degrees'
  write(*,*) 'Resolution    :',(dang*180./pi)*60,'arc min.'

  !read z-lens
  open(11,file=fn_z)  
  do i=1,nslice
     read(11,*) z_write(i)
  enddo
  close(11)

  !read z_sources
  open(11,file=fn_z_s)  
  do i=1,nslice
     read(11,*) z_write_s(i)
  enddo
  close(11)

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

     icount=0
     
     !Loop over the slices
     write(*,*) '****************'
     write(*,*) 'sbox = ', sbox
     write(*,*) '****************'

!     write(*,*) '*************** Doing only a few slices for testing *****************'
!     do j = 8,nslice
     
     do j=1,nslice

        if(j==nslice_far+1)then  
           write(*,*) '****************'
           write(*,*) 'switching to lbox = ', lbox
           write(*,*) '****************'
        endif

        i=mod(icount,3)

        ! Get the Redshift
        write(z_string,'(f7.3)') z_write(j)
        z_string=adjustl(z_string)

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
        call GetArg(1,LOS_str)  
        ! Enforcing to read always LOS51 ?
        !if(z_write(j) <= z_switch)fn=proj_path//"51"//"/"//z_string(1:len_trim(z_string)) 
        !if(z_write(j) > z_switch)fn=proj_path2//"51"//"/"//z_string(1:len_trim(z_string)) 
        if(z_write(j) <= z_switch)fn=proj_path//trim(LOS_str)//"/"//z_string(1:len_trim(z_string)) 
        if(z_write(j) > z_switch)fn=proj_path2//trim(LOS_str)//"/"//z_string(1:len_trim(z_string)) 


#endif
!########################################

        fn=adjustl(fn)

        !if(z_write(j) < 0.050)then

        !use always the same orientation, say xy for now...
        !if(i.eq.0)fp=fn(1:len_trim(fn))//'proj_xy_test.dat'
        !if(i.eq.1)fp=fn(1:len_trim(fn))//'proj_yz_test.dat'
        !if(i.eq.2)fp=fn(1:len_trim(fn))//'proj_xz_test.dat'
        !fp=fn(1:len_trim(fn))//'proj_xy.dat'         

        !else

        if(i.eq.0)fp=fn(1:len_trim(fn))//'proj_half_xy.dat'
        if(i.eq.1)fp=fn(1:len_trim(fn))//'proj_half_yz.dat'
        if(i.eq.2)fp=fn(1:len_trim(fn))//'proj_half_xz.dat'


        !if(i.eq.0)fp=fn(1:len_trim(fn))//'proj_xy.dat'
        !if(i.eq.1)fp=fn(1:len_trim(fn))//'proj_yz.dat'
        !if(i.eq.2)fp=fn(1:len_trim(fn))//'proj_xz.dat'

        !endif

        !if(i.eq.0)fp=fn(1:len_trim(fn))//'proj_xy.dat_subscale'
        !if(i.eq.1)fp=fn(1:len_trim(fn))//'proj_yz.dat_subscale'
        !if(i.eq.2)fp=fn(1:len_trim(fn))//'proj_xz.dat_subscale'
        !if(i.eq.0)fp=fn(1:len_trim(fn))//'proj_xy-'//LOS//'.dat'
        !if(i.eq.1)fp=fn(1:len_trim(fn))//'proj_yz-'//LOS//'.dat'
        !if(i.eq.2)fp=fn(1:len_trim(fn))//'proj_xz-'//LOS//'.dat'
        fu=10+i

        patch=0 !0  ! -1 means don't keep track of shift. 
                     ! 0 means write shift to file
        open(unit=fu,file=fp,form=input_form,status='old',iostat=file_status)
        if (file_status.gt.0) then 
            close(fu)
            write(*,*) 'Failed to open file ', fp
            write(*,*) 'file_status',file_status
            write(*,*) 'Trying to open another file, modifying random shifting...'
            write(*,*) 'Currently set for CFHT sims...'
            call GetArg(2,newLOS_str)
            !read(newLOS_str,*) newLOS
            write(*,*) 'Switching ', LOS, ' for ',newLOS_str
            !newLOS = newLOS-1
            !write(newLOS_str,'(f5.3)') newLOS
            !write(*,*) LOS ,newLOS,newLOS_str

            if(z_write(j) <= z_switch)then
               !fn = '/scratch/p/pen/haoran/ACT_lensing/LOS'//trim(newLOS_str)//'/'//z_string(1:len_trim(z_string))
               fn = '/scratch/p/pen/jharno/Lensing_2048/LOS'//trim(newLOS_str)//'/'//z_string(1:len_trim(z_string))
            endif
            if(z_write(j) >= z_switch)then
               !fn = '/scratch/p/pen/haoran/ACT_lensing/LOS'//trim(newLOS_str)//'/'//z_string(1:len_trim(z_string))
               fn = '/scratch/p/pen/jharno/NewLensing_2048/LOS'//trim(newLOS_str)//'/'//z_string(1:len_trim(z_string))
            endif

            if(i.eq.0)fp=fn(1:len_trim(fn))//'proj_xy.dat'
            if(i.eq.1)fp=fn(1:len_trim(fn))//'proj_yz.dat'
            if(i.eq.2)fp=fn(1:len_trim(fn))//'proj_xz.dat'

            open(unit=fu,file=fp,form=input_form,status='old',iostat=file_status)
            write(*,*) 'Opened', fp
            patch=1 

        else
           write(*,*) 'Opened', fp
        endif        

        !Read the files
        read(fu) ScaleFactor
        write(*,*) 'Got Scalefactor  = ', ScaleFactor


        if(z_write(j) <= z_switch) then
            read(fu) rho_pxy_close
            write(*,*) 'Read', fp
            close(fu)
            !pause
            map_close(:,:,j-nslice_far)=rho_pxy_close
        else
            read(fu) rho_pxy_far
            !write(*,*) 'Read', fp
            close(fu)
            !pause
            map_far(:,:,j)=rho_pxy_far
        endif 
        write(*,*) 'Got projection!'

#ifdef test_read
        write (fn,'(f5.3,"test.dat")') z_write(j)
        open(10,file=Lens_output_path//trim(LOS_str)//"/"//fn, access = output_access, form = output_form, status = 'replace')
        !open(10,file=Lens_output_path//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(npc),4)//"/"//fn, access = output_access, form = output_form, status = 'replace')
        write(10) map_close(:,:,j)
        close(10)
        write(*,*)'Wrote test Map for z=',z_write(j)
        stop
#endif

!------------------------
! Randomly shift the map
!------------------------
        if(write_shift.eq.1)then
           call random_number(shift%x)
           call random_number(shift%y)           

           if(patch==0)then
        !      test_str=Lens_output_path//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(npc),4)
         !   write(*,*) 'Test string' ,trim(test_str)
           !---------------------------
           ! Replace the old shift file
           if(j.eq.1)then
              open(100,file=Lens_output_path//trim(LOS_str)//"/random_shift")!access = output_access, form = output_form, status = 'replace')
              !open(100,file=Lens_output_path//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(npc),4)//"/random_shift")!access = output_access, form = output_form, status = 'replace')
              write(100,*) shift%x, shift%y
              close(100)
           else
              open(100,file=Lens_output_path//trim(LOS_str)//"/random_shift", position = 'append')
              !open(100,file=Lens_output_path//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(npc),4)//"/random_shift", position = 'append')
              write(100,*) shift%x, shift%y
              close(100)
           endif
           !----------------------------

           elseif(patch .eq. 1)then

              !-----------------------------------
              ! must read the other file's random shift and change it first!
              open(101,file='/scratch/p/pen/jharno/SimulLens_2048_to_1024/LOS'//trim(newLOS_str)//"/random_shift", status='old')
              !open(101,file='/scratch/p/pen/jharno/SimulLens_2048_to_1024/LOS'//trim(newLOS_str)//'/zs_'//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(npc),4)//"/random_shift", status='old')
              
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

              open(100,file=Lens_output_path//trim(LOS_str)//"/random_shift", position = 'append')
              !open(100,file=Lens_output_path//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(npc),4)//"/random_shift", position = 'append')
              write(100,*) shift%x, shift%y
              close(100)
              !--------------------------------------
           else !patch = -1
             write(*,*) 'Not keeping track of random shift'
           endif


        elseif(read_shift.eq.1)then
           open(100,file=Lens_output_path//trim(LOS_str)//"/random_shift", status='old',position='asis')!, access = output_access, form = output_form)
           !open(100,file=Lens_output_path//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(npc),4)//"/random_shift", status='old',position='asis')!, access = output_access, form = output_form)
           ! Sloppy way to get the right entry number
           !if(j .le. 2) then
           !   shift%x = 0;
           !   shift%y = 0;
           !else
           !   do j2= 1,(j-1)/2 
              do j2= 1,j          
                 read(100,*) shift%x, shift%y
              enddo
           !endif
           close(100)
        elseif(no_shift.eq.1)then
           shift%x = 0;
           shift%y = 0;
        else
           write(*,*) 'Must specify a shift type. (No shift,  read from or write to file)'
        endif     
           
        write(*,*)'Random shift ', shift%x, shift%y
!--------------
! Done shifting
!--------------        


!********************
#ifdef HALO_CATALOGUE
!********************

        ip = 0
        nh_total = 0
        nh_final = 0
        !do ii=1,nn
           ii = 1
           nh_final_local=0

           ! The following is for parallel catalogs

           !write(*,*) 'Reading Node ',ii           
           !write(rank_s,'(i4)') ii-1
           !rank_s=adjustl(rank_s)
           !write(zstring,'(f5.3)') z3dps       ! Need (f6.3) for Z > 10
           !halofn=fn(1:len_trim(fn))//'halo'//rank_s(1:len_trim(rank_s))//'.dat'        
           !new_halofn = Lens_output_path//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(npc),4)//"/"//z_string(1:len_trim(z_string))//'zoomed_halo'//rank_s(1:len_trim(rank_s))//".dat"
           
           ! This is for merged catalogs

           halofn=fn(1:len_trim(fn))//'halo.dat'        
           new_halofn = Lens_output_path//trim(LOS_str)//"/"//z_string(1:len_trim(z_string))//'zoomed_halo.dat'
           !new_halofn = Lens_output_path//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(npc),4)//"/"//z_string(1:len_trim(z_string))//'zoomed_halo.dat'

           fu2=20+ii
           fu3=30+ii

           !open(unit=fu2,file=halofn,form=input_form,status='old',iostat=file_status) ! for read by node
           open(unit=fu2,file=halofn,form='formatted',status='old',iostat=file_status) ! for reading merged catalogs
           write(*,*) 'Opened', halofn
           open(unit=fu3,file=new_halofn, access = output_access, form = output_form, status = 'replace',iostat=file_status)
           write(*,*) 'Opened', new_halofn
          

           ! Merged halos do not have a header
           !read(fu2) nploc(ii)
           !write(fu3) nploc(ii)
           !write(*,*) 'nploc(ii) =',nploc(ii)
           do
              !nh_total=nh_total+1
              !write(*,*)'Reading halo', nh_total + 1
#ifdef HALOPID
              read(fu2,'(28f20.4)',end=112,err=113) halo_input_buffer, halo_pid!(:,nh_total)           ! for reading merged catalogues
              !read(fu2,end=112,err=113) halo_input_buffer, halo_pid!(:,nh_total)           ! for read by node
#else
              read(fu2,'(28f20.4)', end=112,err=113) halo_input_buffer ! for reading merged catalogues
              !read(fu2, end=112,err=113) halo_input_buffer ! for read by nodes

#endif              
              if(maxval(halo_input_buffer)==0) then
#ifdef debug_halo
                 !write(*,*) 'No halos in  node', ii 
                 !write(*,*) 'nh_total =' , nh_total
#endif
                 !cycle
                 exit
              else
                 nh_total=nh_total+1
              endif

              !write(*,*) 'For testing, keep only 2000 halos'
              !if(nh_total>2000)exit
              



              !********************************
              ! 1- Rotate the halo catalogues *
              !********************************
              
              !***************
              if(i .eq. 0)then
                 !no rotation
                 halo_output_buffer = halo_input_buffer
#ifdef debug_halo_intense
                 write(*,*)'halo_input_buffer(4:6)', halo_input_buffer(4:6)
#endif                 
              endif

              !***************
              if(i .eq. 1)then                 
                 !rotate about y axis by pi/2, counter-clockwise

#ifdef debug_halo_intense
                 write(*,*)'halo_input_buffer(4:6)', halo_input_buffer(4:6)
#endif

                 call Rotate(halo_input_buffer(1:3),0.,-pi/2,0.,halo_output_buffer(1:3) )
                 call Rotate(halo_input_buffer(4:6),0.,-pi/2,0.,halo_output_buffer(4:6) )
                 call Rotate(halo_input_buffer(7:9),0.,-pi/2,0.,halo_output_buffer(7:9) )
                 call Rotate(halo_input_buffer(10:12),0.,-pi/2,0.,halo_output_buffer(10:12) )
                 call Rotate(halo_input_buffer(13:15),0.,-pi/2,0.,halo_output_buffer(13:15) )
                 call Rotate(halo_input_buffer(20:22),0.,-pi/2,0.,halo_output_buffer(20:22) )
                 halo_output_buffer(20:22) = abs(halo_output_buffer(20:22))

                 ! For the I_ij matrix, I need to rotate the three rows, and keep only the relevant ones
                 I_ij_vector_in(:) = halo_input_buffer(23:25)
                 call Rotate(I_ij_vector_in,0.,-pi/2,-pi/2.,halo_output_buffer(23:25))            !I_xx,I_xy, I_xz
                 halo_output_buffer(23) = abs(halo_output_buffer(23))                                !I_xx = positive definite

                 I_ij_vector_in(1) =  halo_input_buffer(24)
                 I_ij_vector_in(2) =  halo_input_buffer(26)
                 I_ij_vector_in(3) =  halo_input_buffer(27)
                 call Rotate(I_ij_vector_in,0.,-pi/2,-pi/2.,I_ij_vector_out)
                 !halo_output_buffer(24) = I_ij_vector_out(1)
                 halo_output_buffer(26) = I_ij_vector_out(2)
                 halo_output_buffer(27) = I_ij_vector_out(3)
                 halo_output_buffer(26) = abs(halo_output_buffer(26))                                !I_yy = positive definite


                 I_ij_vector_in(1) =  halo_input_buffer(25)
                 I_ij_vector_in(2) =  halo_input_buffer(27)
                 I_ij_vector_in(3) =  halo_input_buffer(28)
                 call Rotate(I_ij_vector_in,0.,-pi/2,-pi/2.,I_ij_vector_out)
                 !halo_output_buffer(25) = I_ij_vector_out(1)
                 !halo_output_buffer(27) = I_ij_vector_out(2)
                 halo_output_buffer(28) = I_ij_vector_out(3)
                 halo_output_buffer(28) = abs(halo_output_buffer(28))                                !I_zz = positive definite



#ifdef debug_halo_intense
                 write(*,*)'halo_output_buffer(4:6)', halo_output_buffer(4:6)
#endif
                                  

                 !halo_output_buffer(1) = -halo_input_buffer(3)
                 !halo_output_buffer(2) =  halo_input_buffer(2)
                 !halo_output_buffer(3) =  halo_input_buffer(1)
                 
                 !halo_output_buffer(4) = -halo_input_buffer(6)
                 !halo_output_buffer(5) =  halo_input_buffer(5)
                 !halo_output_buffer(6) =  halo_input_buffer(4)
                 
                 !halo_output_buffer(7) = -halo_input_buffer(9)
                 !halo_output_buffer(8) =  halo_input_buffer(8)
                 !halo_output_buffer(9) =  halo_input_buffer(7)

                 !halo_output_buffer(10) = -halo_input_buffer(12)
                 !halo_output_buffer(11) =  halo_input_buffer(11)
                 !halo_output_buffer(12) =  halo_input_buffer(10)

                 halo_output_buffer(16:19) =  halo_input_buffer(16:19)

                 

                 !halo_output_buffer(18) = -halo_input_buffer(20)
                 !halo_output_buffer(19) =  halo_input_buffer(19)
                 !halo_output_buffer(20) =  halo_input_buffer(18)
              endif


              !***************
              if(i .eq. 2)then                 
                 !rotate about x axis by pi/2, clockwise

#ifdef debug_halo_intense
                 write(*,*)'halo_input_buffer(4:6)', halo_input_buffer(4:6)
#endif

                 call Rotate(halo_input_buffer(1:3),pi/2.,-pi/2,-pi/2.,halo_output_buffer(1:3) )     !x_peak
                 call Rotate(halo_input_buffer(4:6),pi/2.,-pi/2,-pi/2.,halo_output_buffer(4:6) )     !x_cm
                 call Rotate(halo_input_buffer(7:9),pi/2.,-pi/2,-pi/2.,halo_output_buffer(7:9) )     !v_cm
                 call Rotate(halo_input_buffer(10:12),pi/2.,-pi/2,-pi/2.,halo_output_buffer(10:12) ) !l_cm
                 call Rotate(halo_input_buffer(13:15),pi/2.,-pi/2,-pi/2.,halo_output_buffer(13:15) ) !v2_disp
                 call Rotate(halo_input_buffer(20:22),pi/2.,-pi/2,-pi/2.,halo_output_buffer(20:22) ) !var2
                 halo_output_buffer(20:22) = abs(halo_output_buffer(20:22))                          !these need to remain positive

                 ! For the I_ij matrix, I need to rotate the three rows, and keep only the relevant ones
                 I_ij_vector_in(:) = halo_input_buffer(23:25)
                 call Rotate(I_ij_vector_in,pi/2.,-pi/2,-pi/2.,halo_output_buffer(23:25))            !I_xx,I_xy, I_xz
                 halo_output_buffer(23) = abs(halo_output_buffer(23))                                !I_xx = positive definite

                 I_ij_vector_in(1) =  halo_input_buffer(24)
                 I_ij_vector_in(2) =  halo_input_buffer(26)
                 I_ij_vector_in(3) =  halo_input_buffer(27)
                 call Rotate(I_ij_vector_in,pi/2.,-pi/2,-pi/2.,I_ij_vector_out)
                 !halo_output_buffer(24) = I_ij_vector_out(1)
                 halo_output_buffer(26) = I_ij_vector_out(2)
                 halo_output_buffer(27) = I_ij_vector_out(3)
                 halo_output_buffer(26) = abs(halo_output_buffer(26))                                !I_yy = positive definite


                 I_ij_vector_in(1) =  halo_input_buffer(25)
                 I_ij_vector_in(2) =  halo_input_buffer(27)
                 I_ij_vector_in(3) =  halo_input_buffer(28)
                 call Rotate(I_ij_vector_in,pi/2.,-pi/2,-pi/2.,I_ij_vector_out)
                 !halo_output_buffer(25) = I_ij_vector_out(1)
                 !halo_output_buffer(27) = I_ij_vector_out(2)
                 halo_output_buffer(28) = I_ij_vector_out(3)
                 halo_output_buffer(28) = abs(halo_output_buffer(28))                                !I_zz = positive definite




#ifdef debug_halo_intense
                 write(*,*)'halo_output_buffer(4:6)', halo_output_buffer(4:6)
#endif

                 !halo_output_buffer(1) =  halo_input_buffer(1)
                 !halo_output_buffer(2) = -halo_input_buffer(3)
                 !halo_output_buffer(3) =  halo_input_buffer(2)
                 
                 !halo_output_buffer(4) =  halo_input_buffer(4)
                 !halo_output_buffer(5) = -halo_input_buffer(6)
                 !halo_output_buffer(6) =  halo_input_buffer(5)
                 
                 !halo_output_buffer(7) =  halo_input_buffer(7)
                 !halo_output_buffer(8) = -halo_input_buffer(9)
                 !halo_output_buffer(9) =  halo_input_buffer(8)

                 !halo_output_buffer(10) =  halo_input_buffer(10)
                 !halo_output_buffer(11) = -halo_input_buffer(12)
                 !halo_output_buffer(12) =  halo_input_buffer(11)

                 halo_output_buffer(16:19) =  halo_input_buffer(16:19)

                 !halo_output_buffer(18) =  halo_input_buffer(18)
                 !halo_output_buffer(19) = -halo_input_buffer(20)
                 !halo_output_buffer(20) =  halo_input_buffer(19)
              endif
  
              ! Apply random shift as well, then use periodic BC:             

              do jj = 1,6

                 if(jj==1 .or. jj==4)halo_output_buffer(jj) = halo_output_buffer(jj) + shift%x
                 if(jj==2 .or. jj==5)halo_output_buffer(jj) = halo_output_buffer(jj) + shift%y

                 if(halo_output_buffer(jj) < 0.0) then
                    halo_output_buffer(jj) = halo_output_buffer(jj) + nc
                 elseif(halo_output_buffer(jj) >nc) then
                    halo_output_buffer(jj) = halo_output_buffer(jj) - nc
                 endif
              enddo

#ifdef debug_halo_intense
              write(*,*)'halo_output_buffer(4:6)', halo_output_buffer(4:6)
#endif

              !*********************************
              ! 2- Find magnification and zoom *
              !*********************************

              ! Get chi(halo)
              chi_halo = chi(z_write(j),omegam,h)*h + (halo_output_buffer(6) - nc/2)*(box/nc)

#ifdef debug_halo_intense
              write(*,*) 'chi(z_s) =',chi(zs,omegam,h)*h, 'Mpc/h'
              write(*,*) 'chi(z_i) =',chi(z_write(j),omegam,h)*h, 'Mpc/h'
              write(*,*) 'chi_halo =',chi_halo, 'Mpc/h'
#endif

              if(chi_halo>chi(zs,omegam,h)*h) then
                 !write(*,*)'Halo behind source'
                 cycle
              endif

  
              !*** Conical geometry ***!
              ! Interpolate npc pixels in a continuous conical volume, 
              ! Find distance to center of lens, then zoom on that volume
              !mag = ((angle/2)*chi_halo*(nc/box))

              !*** Blocks geometry ***
              ! Interpolate npc pixels in a volume that matches the projected
              ! densities. Use the lens distance for all halos.

              mag =  ((angle/2)*chi(z_write(j),omegam,h)*(nc/box))

#ifdef debug_halo_intense
              !write(*,*) 'mag = ',mag*box/nc, 'Mpc/h'
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
             
 
              halo_output_buffer(1)  = hc/(mag)*(halo_output_buffer(1) - hc) + hc
              halo_output_buffer(2)  = hc/(mag)*(halo_output_buffer(2) - hc) + hc
              halo_output_buffer(4)  = hc/(mag)*(halo_output_buffer(4) - hc) + hc
              halo_output_buffer(5)  = hc/(mag)*(halo_output_buffer(5) - hc) + hc

              halo_output_buffer(7)  = hc/(mag)*halo_output_buffer(7)
              halo_output_buffer(8)  = hc/(mag)*halo_output_buffer(8)
              halo_output_buffer(10) = hc/(mag)*halo_output_buffer(10)
              halo_output_buffer(11) = hc/(mag)*halo_output_buffer(11)
              halo_output_buffer(13) = hc/(mag)*halo_output_buffer(13) ! not sure though... v_disp is more complicated...
              halo_output_buffer(14) = hc/(mag)*halo_output_buffer(14)

              halo_output_buffer(16) = (hc/(mag))**2*halo_output_buffer(16)

              halo_output_buffer(20) = (hc/(mag))**2*halo_output_buffer(20)
              halo_output_buffer(21) = (hc/(mag))**2*halo_output_buffer(21)

              !*** Blocks geometry ***
              ! Interpolate npc pixels in a volume that matches the projected
              ! densities

#ifdef debug_halo_intense
              write(*,*)'x_after mag=',halo_output_buffer(4) 
              write(*,*)'y_after mag=',halo_output_buffer(5) 
              !write(*,*)'var_x_after=',halo_output_buffer(20)
              !write(*,*)'var_y_after=',halo_output_buffer(21)              
#endif

              ! if out of box, don't write!              
              if((halo_output_buffer(4) < nc) .and. (halo_output_buffer(4).ge.0) .and. &
                   (halo_output_buffer(5) < nc) .and. (halo_output_buffer(5).ge.0))then

                 nh_final_local = nh_final_local+1

#ifdef HALOPID
                 !write(*,*) 'Writing to zoomed_halo file'
                 write(fu3) halo_output_buffer, halo_pid!(:,nh_total)           
#else
                 !write(*,*) 'Writing to zoomed_halo file'
                 write(fu3) halo_output_buffer
#endif                               
              endif



              !peak(:,nh_total)=halo_input_buffer(1:3)
              !xv(:,nh_total)=halo_input_buffer(4:9)
              !l(:,nh_total)=halo_input_buffer(10:12)
              !v_disp(nh_total)=halo_input_buffer(13)
              !radius_scale(nh_total)=halo_input_buffer(14)
              !halo_mass(nh_total)=halo_input_buffer(15)
              !halo_mass_pp(nh_total)=halo_input_buffer(16)
              !halo_mass1(nh_total)=halo_input_buffer(17)
              !var(:,nh_total)=halo_input_buffer(18:20)
              
#ifdef debug_halo_intense
              write(*,*)'peak=', halo_output_buffer(1:3)
              write(*,*)'xv=',halo_output_buffer(4:9)
              write(*,*)'l=',halo_output_buffer(10:12)
              write(*,*)'v_disp=',halo_output_buffer(13:15)
              write(*,*)'radius_scale=',halo_output_buffer(16)
              write(*,*)'halo_mass=',halo_output_buffer(17)
              write(*,*)'halo_mass_pp=',halo_output_buffer(18)
              write(*,*)'halo_mass1=',halo_output_buffer(19)
              write(*,*)'var=',halo_output_buffer(20:22)              
              write(*,*)'I_ij=', halo_output_buffer(23:28)             
#ifdef HALOPID
              write(*,*)'PID=',halo_pid!(:,nh_total)   
#endif              

#endif
             
              !stop

113           continue
           enddo
112        close(fu2)
           close(fu3)
           write(*,*)'closed original and new halo files'

           nh_final = nh_final + nh_final_local
           ip=ip+nploc(ii)
           write(*,*) 'np cumulative = ', ip,', np local = ', nploc(ii)
           write(*,*) 'nh final cumulative = ', nh_final ,'nh final local = ', nh_final_local
        !enddo ! end loop over nodes

        !call cic(xv,d,nc,ip)

#endif
!********************
! done HALO_CATALOGUE      
!********************

#ifdef cubepm
        !! substract the mean to get the overdensity
        write(*,*) 'substracting mean'
        !write(*,*) 'test:', map_close(1,1,1), map_close(1,1,2), map_close(1,1,3), sum(map_close(1,1,1:3))
        if(z_write(j) <= z_switch) then
           rhomean=sum(real(map_close(:,:,j-nslice_far),kind=8))/nc/nc 
           map_close(:,:,j-nslice_far)=map_close(:,:,j-nslice_far)-rhomean  
        else
           rhomean=sum(real(map_far(:,:,j),kind=8))/nc/nc
           map_far(:,:,j)=map_far(:,:,j)-rhomean !nc 
        endif
        write(*,*) 'Mean before subtraction= ',rhomean

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
        if(z_write(j) <= z_switch) then
           write(*,*)'before kappa2defl delta mean,min,max:' , sum(map_close(:,:,j-nslice_far)/nc/nc), minval(map_close(:,:,j-nslice_far)),maxval(map_close(:,:,j-nslice_far)) 
           call kappa_to_defl(map_close(:,:,j-nslice_far),defl,shear,phi,nc) !Using discrete differentiation
           write(*,*)'after kappa2defl delta mean,min,max:' , sum(map_close(:,:,j-nslice_far)/nc/nc), minval(map_close(:,:,j-nslice_far)),maxval(map_close(:,:,j-nslice_far)) 
           !call kappa2gamma(map_close(:,:,j-nslice_far),defl,shear,phi,nc)  !Using continuous differentiation
        else
           write(*,*)'before kappa2defl delta mean,min,max:' , sum(map_far(:,:,j)/nc/nc), minval(map_far(:,:,j)),maxval(map_far(:,:,j)) 
           call kappa_to_defl(map_far(:,:,j),defl,shear,phi,nc) !Using discrete differentiation
           !call kappa2gamma(map_far(:,:,j),defl,shear,phi,nc)  !Using continuous differentiation
           write(*,*)'after kappa2defl delta mean,min,max:' , sum(map_far(:,:,j)/nc/nc), minval(map_far(:,:,j)),maxval(map_far(:,:,j)) 
        endif
        write(*,*) 'Calculated shear and deflection'
        !write(*,*)' Unzoomed shear p mean,min,max:' , sum(shear(:,:)%p/nc/nc), minval(shear(:,:)%p),maxval(shear(:,:)%p) 
        !write(*,*)' Unzoomed shear c mean,min,max:' , sum(shear(:,:)%c/nc/nc), minval(shear(:,:)%c),maxval(shear(:,:)%c) 

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
       
        if(z_write(j) <= z_switch) then
           call zoomshiftmap(map_close(:,:,j-nslice_far),map(:,:,j),nc,npc,shift,newdefl,frac)
        else
           call zoomshiftmap(map_far(:,:,j),map(:,:,j),nc,npc,shift,newdefl,frac)
        endif
        write(*,*)'delta mean,min,max:' , sum(map(:,:,j)/nc/nc), minval(map(:,:,j)),maxval(map(:,:,j)) 

!##############
#ifdef calshear
#ifndef  SemiBorn


          write(*,*)'shear p mean,min,max:' , sum(shear(:,:)%p/nc/nc), minval(shear(:,:)%p),maxval(shear(:,:)%p) 
          write(*,*)'shear c mean,min,max:' , sum(shear(:,:)%c/nc/nc), minval(shear(:,:)%c),maxval(shear(:,:)%c) 

          call zoomshiftmap(shear%p,newshear(:,:,j)%p,nc,npc,shift,CorrBornDefl,frac)
          call zoomshiftmap(shear%c,newshear(:,:,j)%c,nc,npc,shift,CorrBornDefl,frac)
          call zoomshiftmap(defl%x,newdefl(:,:,j)%x,nc,npc,shift,CorrBornDefl,frac)
          call zoomshiftmap(defl%y,newdefl(:,:,j)%y,nc,npc,shift,CorrBornDefl,frac)
          
          ! Scale the angles
          newdefl(:,:,j)%x=newdefl(:,:,j)%x*(box/nc)/chi(z_write(j),omegam,h)/h
          newdefl(:,:,j)%y=newdefl(:,:,j)%y*(box/nc)/chi(z_write(j),omegam,h)/h
         
          write(*,*)'newshear p mean,min,max:' , sum(newshear(:,:,j)%p/nc/nc), minval(newshear(:,:,j)%p),maxval(newshear(:,:,j)%p) 
          write(*,*)'newshear c mean,min,max:' , sum(newshear(:,:,j)%c/nc/nc), minval(newshear(:,:,j)%c),maxval(newshear(:,:,j)%c) 

          !CorrBornDefl(:,:)%x=newdefl(:,:,1)%x
          !CorrBornDefl(:,:)%y=newdefl(:,:,1)%y

#endif
#endif
          write(*,*) 'Randomly Shifted and Zoomed on the  maps'        

        icount=icount+1
!        if(mod(j,2) .eq. 0) icount=icount+1
!        write(*,*) '******TEST: incrementing rotation every two steps...'
     enddo !!j=1,nslice

 


     !*****************
     !loop over sources
     !*****************
     do i=nslice,1,-1

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
        cumul_gamma1(:,:)=0.0
        cumul_gamma2(:,:)=0.0
        cumul_deflx(:,:)=0.0
        cumul_defly(:,:)=0.0

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

           
           cumul_kappa(:,:) =cumul_kappa(:,:) +map(:,:,j)*kernel


!##############
#ifdef calshear

#ifdef SemiBorn
           !In Post Born Approx, kappa is still calculated at pixel location.
           !Only gamma's and defl's are calulated at light ray location.

           write(*,*)'shear p mean,min,max:' , sum(shear(:,:)%p/nc/nc), minval(shear(:,:)%p),maxval(shear(:,:)%p) 
           write(*,*)'shear c mean,min,max:' , sum(shear(:,:)%c/nc/nc), minval(shear(:,:)%c),maxval(shear(:,:)%c) 

           call zoomshiftmap(shear%p,newshear(:,:,j)%p,nc,npc,shift,CorrBornDefl,frac)
           call zoomshiftmap(shear%c,newshear(:,:,j)%c,nc,npc,shift,CorrBornDefl,frac)
           call zoomshiftmap(defl%x,newdefl(:,:,j)%x,nc,npc,shift,CorrBornDefl,frac)
           call zoomshiftmap(defl%y,newdefl(:,:,j)%y,nc,npc,shift,CorrBornDefl,frac)
           
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
           
          
       enddo ! Loop over lenses

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
        open(10,file=Lens_output_path//trim(LOS_str)//"/"//fn, access = output_access, form = output_form, status = 'replace')
        !open(10,file=Lens_output_path//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(npc),4)//"/"//fn, access = output_access, form = output_form, status = 'replace')
        write(10) cumul_kappa(:,:)
        close(10)
        write(*,*)'Wrote kappa Map for z=',z_write(i)
#endif
!----------------
#ifdef calshear
        write (fn,'(f5.3,"gamma1_weight.dat")') z_write(i)
        !open(11,file=Lens_output_path//"run_"//digitn(ir,2)//"/"//fn,form='binary')
        open(11,file=Lens_output_path//trim(LOS_str)//"/"//fn, access = output_access, form = output_form, status = 'replace')
        !open(11,file=Lens_output_path//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(npc),4)//"/"//fn, access = output_access, form = output_form, status = 'replace')
        write(11) cumul_gamma1(:,:)
        close(11)
        write (fn,'(f5.3,"gamma2_weight.dat")') z_write(i)
        !open(10,file=Lens_output_path//"run_"//digitn(ir,2)//"/"//fn,form='binary')\
        open(12,file=Lens_output_path//trim(LOS_str)//"/"//fn, access = output_access, form = output_form, status = 'replace')
        !open(12,file=Lens_output_path//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(npc),4)//"/"//fn, access = output_access, form = output_form, status = 'replace')
        write(12) cumul_gamma2(:,:)
        close(12)
        write(*,*)'Wrote shear Maps for z=',z_write(i)

        write (fn,'(f5.3,"deflx_weight.dat")') z_write(i)
        !open(10,file=Lens_output_path//"run_"//digitn(ir,2)//"/"//fn,form='binary')
        open(13,file=Lens_output_path//trim(LOS_str)//"/"//fn, access = output_access, form = output_form, status = 'replace')
        !open(13,file=Lens_output_path//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(npc),4)//"/"//fn, access = output_access, form = output_form, status = 'replace')
        write(13) cumul_deflx(:,:)
        close(13)
        write (fn,'(f5.3,"defly_weight.dat")') z_write(i)
        !open(10,file=Lens_output_path//"run_"//digitn(ir,2)//"/"//fn,form='binary')\
        open(14,file=Lens_output_path//trim(LOS_str)//"/"//fn, access = output_access, form = output_form, status = 'replace')
        !open(14,file=Lens_output_path//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(npc),4)//"/"//fn, access = output_access, form = output_form, status = 'replace')
        write(14) cumul_defly(:,:)
        close(14)
        write(*,*)'Wrote shear Maps for z=',z_write(i)
#endif
!----------------
#endif

!----------------------
#ifdef write_delta_maps
        write(*,*) 'Wrinting delta maps'
        write (fn,'(f5.3,"delta.dat_bicubic")') z_write(i)
        !open(10,file=Lens_output_path//"run_"//digitn(ir,2)//"/"//fn,form='binary')
        open(15,file=Lens_output_path//trim(LOS_str)//"/"//fn, access = output_access, form = output_form, status = 'replace')
        !open(15,file=Lens_output_path//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(npc),4)//"/"//fn, access = output_access, form = output_form, status = 'replace')
        write(15) map(:,:,i)
        close(15)
        write(*,*)'Wrote delta Map for z=',z_write(i)
#endif
!-------------------

!#####
#endif
!##### 
!end of 'ifdef z_slice'
 

     enddo ! i=1,nslice

     write(*,*) 'Done integral over redshift slices'


!###############   
#ifdef integrate      
#ifdef full_geometry
     !open(20,file=Lens_output_path//fn_kappa,form='binary')
     open(20,file=Lens_output_path//trim(LOS_str)//"/"//fn_kappa, access = output_access, form = output_form, status = 'replace')
     !open(20,file=Lens_output_path//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(npc),4)//"/"//fn_kappa, access = output_access, form = output_form, status = 'replace')
     write(20)cumul_kappa(:,:)
     close(20)
     write(*,*)'Wrote integrated kappa Maps'

#ifdef calshear
     open(20,file=Lens_output_path//trim(LOS_str)//"/"//fn_gamma1, access = output_access, form = output_form, status = 'replace')
     !open(20,file=Lens_output_path//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(npc),4)//"/"//fn_gamma1, access = output_access, form = output_form, status = 'replace')
     write(20) cumul_gamma1(:,:)
     close(20)
     !finalmap=0
     !do i=1,nslice
     !   if(z_write(i).le.zs)finalmap(:,:)=finalmap(:,:)+newshear(:,:,i)%c
     !enddo
     !open(20,file=Lens_output_path//fn_gamma2,form='binary')
     open(20,file=Lens_output_path//trim(LOS_str)//"/"//fn_gamma2, access = output_access, form = output_form, status = 'replace')
     !open(20,file=Lens_output_path//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(npc),4)//"/"//fn_gamma2, access = output_access, form = output_form, status = 'replace')
     write(20) cumul_gamma2(:,:)
     close(20)     
     !open(20,file=Lens_output_path//fn_defl,form='binary')
     open(20,file=Lens_output_path//trim(LOS_str)//"/"//fn_defl, access = output_access, form = output_form, status = 'replace')
     !open(20,file=Lens_output_path//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(npc),4)//"/"//fn_defl, access = output_access, form = output_form, status = 'replace')
     write(20)cumul_deflx(:,:),cumul_defly(:,:)
     close(20)

     write(*,*)'Wrote Integrated shear Maps'
#endif

!####################
#ifdef power_spectrum
!-----------------
#ifdef write_kappa 
     power = 0.0
     call ps2(cumul_kappa,cumul_kappa,power,npc)
     !open(30,file=Lens_output_path//'l2cl_kappa.'//digitn(ir,2))
     open(30,file=Lens_output_path//trim(LOS_str)//'/l2cl_kappa.'//digitn(ir,2))
     !open(30,file=Lens_output_path//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(npc),4)//'/l2cl_kappa.'//digitn(ir,2))
     do i=1,npc/2
        write(30,*) i*2*pi/angle,power(i)
     enddo
     close(30)
#endif
!-----------------
#ifdef calshear
     power = 0.0
     call ps2(cumul_gamma1,cumul_gamma1,power,npc)
     !open(30,file=Lens_output_path//'l2cl_kappa.'//digitn(ir,2))
     open(30,file=Lens_output_path//trim(LOS_str)//'/l2cl_gamma1.'//digitn(ir,2))
     !open(30,file=Lens_output_path//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(npc),4)//'/l2cl_gamma1.'//digitn(ir,2))
     do i=1,npc/2
        write(30,*) i*2*pi/angle,power(i)
     enddo
     close(30)
     power = 0.0
     call ps2(cumul_gamma2,cumul_gamma2,power,npc)
     !open(30,file=Lens_output_path//'l2cl_kappa.'//digitn(ir,2))
     open(30,file=Lens_output_path//trim(LOS_str)//'/l2cl_gamma2.'//digitn(ir,2))
     !open(30,file=Lens_output_path//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(npc),4)//'/l2cl_gamma2.'//digitn(ir,2))
     do i=1,npc/2
        write(30,*) i*2*pi/angle,power(i)
     enddo
     close(30)
#endif
!-----------------
     write(*,*)'Wrote Power Spectrum'
#endif  
!####################   
#else
     write(*,*) 'Cannot integrate without full geometry'
#endif
!end of full geometry
#endif
!end of integrate

  enddo !! ir=1,nr

#ifdef mix_nbody_runs
  open(150,file=Lens_output_path//trim(LOS_str)//'/count_run.dat', access = 'sequential')
  !open(150,file=Lens_output_path//"zs_"//digitn(int(zs),1)//"_"//digitn(int(lbox),3)//"_mix_"//digitn(int(sbox),3)//"Mpc_"//digitn(int(npc),4)//'/count_run.dat', access = 'sequential')
!  do i=1,num_nbody_runs
!     write(150,*) i,count_run(i)
!  enddo
  do i=1,nslice
     write(150,*) rand_LOS(i)
  enddo
  close(150)
#endif

end program SimulLens

  subroutine zoomshiftmap(map1,map2,nc,npc,shift,defl,frac)
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
   
    write(*,*) 'Calling zoomshiftmap'

    do j1=1,nc
       do i1=1,nc
          map3(i1,j1)=map1(i1,j1)
          map3(i1+nc,j1)=map1(i1,j1)
          map3(i1,j1+nc)=map1(i1,j1)
          map3(i1+nc,j1+nc)=map1(i1,j1)
       enddo
    enddo
        
    call zoom(map3(:,:),map2(:,:),nc,npc,shift,defl,frac)
    !call zoom_bicubic(map3(:,:),map2(:,:),nc,npc,shift,defl,frac)
     
    return
  end subroutine zoomshiftmap 

  subroutine zoom(tmp,map2,n1,n2,shift,defl,frac)
  !! tmp doubles the original slice by the periodic condition
    use Lensing 
    implicit none

    integer n1,n2,i,j,i1,j1,i2,j2,ib,ip,jp,jb
    type(vec2D) shift
    type(vec2D), dimension(n2,n2) :: defl 
    real map2(n2,n2),x,y,frac,tmp(2*n1,2*n1),w,frac1,s,w1,w2

    map2=0
    frac1=n1*frac/n2

    !write(*,*) '***************'
    !write(*,*) 'frac = ',frac 
    !write(*,*) '***************'
    !pause

    do j2=1,n2
       do i2=1,n2
          jb=shift%y*n1+(j2-1)*frac1+1+defl(i2,j2)%y*n2
          ib=shift%x*n1+(i2-1)*frac1+1+defl(i2,j2)%x*n2
          w2=shift%y*n1+(j2-1)*frac1+1-jb+defl(i2,j2)%y*n2
          w1=shift%x*n1+(i2-1)*frac1+1-ib+defl(i2,j2)%x*n2
          jb=modulo(jb-1,2*n1)+1
          ib=modulo(ib-1,2*n1)+1
          jp=modulo(jb,2*n1)+1
          ip=modulo(ib,2*n1)+1
          map2(i2,j2)=tmp(ip,jp)*w1*w2+tmp(ip,jb)*w1*(1-w2)+&
                          tmp(ib,jp)*(1-w1)*w2+tmp(ib,jb)*(1-w1)*(1-w2)
       enddo
    enddo

    return
  end subroutine zoom

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
    nthreads =  omp_get_num_threads()
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
  
