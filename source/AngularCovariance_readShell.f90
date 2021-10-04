
!! modified on Hy Trac's + Ting Ting Lu's code.
program cicpow
  implicit none
#include 'par.fh'

  !! np should be set to nc (1:1) or hc (1:2)
  integer, parameter :: nt=1
  integer(kind=8), parameter :: np=hc**3
  real, parameter :: mp=int(nc,kind=8)**3/np !nc**3/np

  !! cubepm
!  integer, parameter :: nodes_dim=2
  real, parameter  :: ncc=nc/nodes_dim !! number of cells / cubic 
  real, parameter  :: rnc = nc

  !! Dark matter arrays
  real, dimension(6,np) :: xv
  integer, dimension(np) :: ll
  integer, dimension(nc,nc) :: hoc
  integer, dimension(2,nc,nc,nt) :: htoc

  !! Power spectrum arrays
  real, dimension(nc) :: ps
  real, dimension(nc) :: PS_Ave,  CrossCorrAve, N_CrossCorrAve !
  real, dimension(nc+2,nc,nc) :: d, CrossCorr, N_CrossCorr!, psShell_Inner, psShell_Outer, x_Hat
 
#ifdef FFT_TEST
 real, dimension(nc+2,nc,nc) :: Sync, Gauss, DeltaFunction, N_Hat, Hat
 real, dimension(nc) :: Xsi, HatAve, GaussAve, DeltaAve, FFT_Hat
#endif

  !real, dimension(0:hc,-hc+1:hc,-hc+1:hc) :: Hat
  integer, parameter :: MSL=100

  complex, dimension(nc/2 + 1,nc,nc) ::  d_COMPLEX, psV_COMPLEX, CrossCorr_COMPLEX, N_CrossCorr_COMPLEX 
  !complex, dimension(nc/2 + 1,nc,nc) ::  psShell_Inner_COMPLEX, psShell_Outer_COMPLEX, N_psShell_Inner_COMPLEX,N_psShell_Outer_COMPLEX

  equivalence (d, d_COMPLEX) 
  equivalence (CrossCorr, CrossCorr_COMPLEX) 
  equivalence (N_CrossCorr, N_CrossCorr_COMPLEX) 

  integer nploc(nn)
  real pi

  !! variables in CubePM
  integer cubepm_nts,cubepm_cur_checkpoint,cubepm_cur_projection,cubepm_cur_halofind,i1,j1,k1, node_coords(3)
  real  cubepm_a,cubepm_t,cubepm_tau,cubepm_dt_f_acc,cubepm_dt_pp_acc,cubepm_dt_c_acc, cubepm_mass_p
  
  !! other parameters
  real MAX_Outer, MAX_Inner
  character*3 Max_Inner_char, Max_Outer_char


  common /rarr/ d,xv,ps 
  common /iarr/ ll,htoc,hoc

  pi=acos(-1.)

#ifndef RANDOM
  call readdm 
  call cic
#endif

  call powerspectrum
  call writeps

#ifdef AngularCovariance
  call CrossCorrelationOfShells
  call writeCov
#endif

contains

  subroutine readdm
    implicit none
    character*100 fn

    character*8 t1,t2
    integer i,j 
    integer(kind=8) ip
    real HubbleScale
    real Conversion
    character (len=4) :: rank_s
    character (len=MSL) :: ofile,zstring
    character (len=MSL) :: ifile

    !! Read particle data file
    ip=0
    do i=1,nn
       write(*,*) 'Reading Node ',i
       write(rank_s,'(i4)') i-1
       rank_s=adjustl(rank_s)
       write(zstring,'(f5.3)') z3dps
       ifile=trim(zstring)//'xv'//rank_s(1:len_trim(rank_s))//".dat"
       ofile=proj_path//Version//'/out/RUN-'//Run//'/'//trim(ifile)
       write(*,*) 'opening ',ofile
       open (unit=12,file=ofile,form='binary')
#ifdef pmfast
       read(12) nploc(i)
#endif
#ifdef cubepm
       read(12) nploc(i), cubepm_a,cubepm_t,cubepm_tau,cubepm_nts,cubepm_dt_f_acc,cubepm_dt_pp_acc,cubepm_dt_c_acc,cubepm_cur_checkpoint, cubepm_cur_projection,cubepm_cur_halofind,cubepm_mass_p
#endif
       read(12) xv(1:6,ip+1:ip+nploc(i))
       close(12)
#ifdef pmfast       
       xv(3,ip+1:ip+nploc(i))=xv(3,ip+1:ip+nploc(i))+(i-1)*nc/nn
#endif
#ifdef cubepm
       do k1=1,nodes_dim
          do j1=1,nodes_dim
             do i1=1,nodes_dim
                if (i-1 == (i1-1)+(j1-1)*nodes_dim+(k1-1)*nodes_dim**2)  &
                    node_coords(:)=(/(i1-1),(j1-1),(k1-1)/)
             enddo
          enddo
       enddo
       xv(1,ip+1:ip+nploc(i))=modulo(xv(1,ip+1:ip+nploc(i))+node_coords(1)*ncc,rnc)
       xv(2,ip+1:ip+nploc(i))=modulo(xv(2,ip+1:ip+nploc(i))+node_coords(2)*ncc,rnc)
       xv(3,ip+1:ip+nploc(i))=modulo(xv(3,ip+1:ip+nploc(i))+node_coords(3)*ncc,rnc)
#endif
#ifdef Kaiser
       !Red Shift Distortion: x_z -> x_z +  v_z/H(Z)
       !Using H0=0.23335, where I have truncated a factor of 10**-17 to avoid NaN

       !Convert seconds into simulation time units: 
       Conversion = 1/(5.8317*cubepm_a**2) !Added a factor of 10**17 to avoid NaN. 
                                           !Already compensated in the definition of
       !write(*,*) 'Conversion = ', Conversion !IT WORKS
       !write(*,*) 'Test : Conversion*H^-1= ', Conversion*1.5151*(10**17)

       !Compute Hubble Scale in matter era:

       HubbleScale = 0.23335*cubepm_a**(-3/2)

       xv(3,:)=xv(3,:) + xv(6,:)*Conversion/HubbleScale

       if(i==nn) then
          write(*,*) 'Included Kaiser Effect'
          write(*,*) 'Test : Scale Factor = ', cubepm_a 
          write(*,*) 'Test : Hubble Scale = ', HubbleScale 
          write(*,*) 'Test : in Simulation time units: = ', HubbleScale*Conversion 
       endif
#endif

       ip=ip+nploc(i)
       !write(*,*) 'np cumulative = ', ip,', np local = ', nploc(i)
    enddo
 
    !write(*,*) 'xv', maxval(xv(1:3,:)),minval(xv(1:3,:)),maxval(xv(4:6,:)),minval(xv(4:6,:))
    write(*,*) '*************'
    write(*,*) '*Done readdm*'
    write(*,*) '*************'
    return
  end subroutine readdm

  !************************
  subroutine writeps
    implicit none
    integer k
    character*150 fn1


    fn1=PS_dir//'z_'//Z//'/cicpow_'//RedShift//'_Run'//Run//'-'//Version//'.dat'

#ifdef Kaiser
    fn1=PS_dir//'z_'//Z//'/cicpow_'//RedShift//'_Run'//Run//'-'//Version//'-RSD.dat'
#endif

#ifdef RANDOM
    fn1=PS_dir//'z_'//Z//'/cicpow_'//RedShift//'_Run'//Run//'-'//Version//'-RDM.dat'
#endif

#ifdef GAUSS
    fn1=PS_dir//'z_'//Z//'/cicpow_'//RedShift//'_Run'//Run//'-'//Version//'-GAUSS.dat'
#ifdef TEST
    fn1=Test_dir//'cicpow_'//RedShift//'_Run'//Run//'-'//Version//'-GAUSS.dat'
#endif
#endif

    !! Output power spectrum
    !! First column is physical k
    !! Second column is \Delta^2

    open(11,file=fn1)
    do k=1,hc !hc+1
       write(11,*) 2*pi/lbox*(k), ps(k) !, exp(-2*ps(1,k)**2/sig**2)/box**3*ps(1,k)**3/2/pi**2*AA**2
    enddo
    close(11)

    write(*,*) 'Wrote ', fn1

    write(*,*) '**************'
    write(*,*) '*done writeps*'
    write(*,*) '**************'

    return
  end subroutine writeps

  !************************

  subroutine writeCov
    implicit none
    integer k
    character*300, fn2

    !fn2=Test_dir//'CrossCorr_'//RedShift//'_Run'//Run//'-'//Version//'.dat'

    fn2=CrossCorr_dir//'z_'//Z//'/CrossCorr_'//RedShift//'_'//MAX_Outer_char//'_'//MAX_Inner_char//'_Run'//Run//'-'//Version//'.dat'
#ifdef Kaiser
    fn2=CrossCorr_dir//'z_'//Z//'/CrossCorr_'//RedShift//'_'//MAX_Outer_char//'_'//MAX_Inner_char//'_Run'//Run//'-'//Version//'-RSD.dat'
#endif

#ifdef RANDOM
    fn2=CrossCorr_dir//'z_'//Z//'/CrossCorr_'//RedShift//'_'//MAX_Outer_char//'_'//MAX_Inner_char//'_Run'//Run//'-'//Version//'-RDM.dat'
#endif

#ifdef GAUSS
    fn2=CrossCorr_dir//'z_'//Z//'/CrossCorr_'//RedShift//'_'//MAX_Outer_char//'_'//MAX_Inner_char//'_Run'//Run//'-'//Version//'-GAUSS.dat'
#ifdef TEST
    fn2=Test_dir//'CrossCorr_'//RedShift//'_'//MAX_Outer_char//'_'//MAX_Inner_char//'_Run'//Run//'-'//Version//'-GAUSS.dat'
#endif
#endif



#ifdef debug   
    write(*,*) 'TEST',MAX_Inner,MAX_Outer,2*MAX_Inner*MAX_Outer,2.0*MAX_Inner*MAX_Outer
    write(*,*) 'TEST',((MAX_Inner**2+MAX_Outer**2-(10-1.)**2)/(2*MAX_Inner*MAX_Outer))
    write(*,*) 'TEST',acos((MAX_Inner**2+MAX_Outer**2-(10-1.)**2)/(2*MAX_Inner*MAX_Outer))
#endif

    open(22,file=fn2)
    do k=int(MAX_Outer-MAX_Inner)+1,int(MAX_Outer+MAX_Inner)+1!1,hc+1
       write(22,*) (k-1),acos((MAX_Inner**2+MAX_Outer**2-(k-1)**2)/(2*MAX_Inner*MAX_Outer)),N_CrossCorrAve(k),CrossCorrAve(k),CrossCorrAve(k)/N_CrossCorrAve(k) 
    enddo
    close(22)

    write(*,*) 'Wrote ', fn2

#ifdef FFT_TEST

    integer r
    character*100 fn8, fn9

    fn8=dir_work//'FFT_XSI.dat'
    fn9=dir_work//'FFT_HAT.dat'


    open(20,file=fn8)
    do r=2,nc! r=1 Asks for FFT_Hat(0) which is a NaN
       write(20,*) r-1, Xsi(r)   
    enddo
    close(20)

    open(21,file=fn9)
    do k=2,hc+1 
       write(21,*) k-1, HatAve(k), FFT_Hat(k),  GaussAve(k), DeltaAve(k) 
    enddo
    close(21)

    write(*,*) 'Wrote ', fn8
    write(*,*) 'Wrote ', fn9


#endif


    write(*,*) '**************'
    write(*,*) '*done writeCov*'
    write(*,*) '**************'

    return
  end subroutine writeCov



  !************************************
  !*** cic computes the density 'd' ***
  !************************************

  subroutine cic
    implicit none
    integer, parameter :: kpt=nc/nt
    integer, parameter :: npt=np/nt

    integer it,i,j,k,ip,OutBound
    real toe

    !! Construct chaining lists in parallel
    !$omp parallel do default(shared) private(it,ip,j,k)
    OutBound = 0
    do it=1,nt
       htoc(:,:,:,it)=0
       do ip=1+(it-1)*npt,min(np,it*npt)
          j=floor(xv(2,ip))+1
          k=floor(xv(3,ip))+1
          if((j > nc)) then
#ifdef debug
             write (*,*) '#### PROBLEM!!! (j = floor(xv(2,ip))+1) =',j
             write (*,*) 'Enforcing Periodic BC Manually'
#endif
             j = j-nc
             OutBound = OutBound +1
             !pause
          endif
          if((k > nc)) then
#ifdef debug
             write (*,*) '#### PROBLEM!!! (k = floor(xv(2,ip))+1) =',k
             write (*,*) 'Enforcing Periodic BC Manually'
#endif
             k = k-nc
             OutBound = OutBound +1
             !pause
          endif


          if (htoc(1,j,k,it) .eq. 0) then
             ll(ip)=0
             htoc(:,j,k,it)=ip
          else
             ll(htoc(2,j,k,it))=ip
             ll(ip)=0
             htoc(2,j,k,it)=ip
          endif
       enddo
    enddo
    !$omp end parallel do

    write(*,*) 'Enforced BC with ',OutBound, 'particles'

    !! Merge chaining lists
    !$omp parallel do default(shared) private(it,j,k,toe)
    do k=1,nc
       do j=1,nc
          hoc(j,k)=0
          do it=1,nt
             if (hoc(j,k) .eq. 0) then
                hoc(j,k)=htoc(1,j,k,it)
                toe=htoc(2,j,k,it)
             else
                if (htoc(1,j,k,it) .ne. 0) then
                   ll(toe)=htoc(1,j,k,it)
                   toe=htoc(2,j,k,it)
                endif
             endif
          enddo
       enddo
    enddo
    !$omp end parallel do


    !! Initialize density field
    !$omp parallel do default(shared) private(it,k)
    do it=1,nt
       do k=1+(it-1)*kpt,min(nc,it*kpt)
          d(:,:,k)=0
       enddo
    enddo
    !$omp end parallel do


    !! Add particle density to density field
!    do ko=1,4
    !$omp parallel do default(shared) private(j,k,ip)
       do k=1,nc
          do j=1,nc
             ip=hoc(j,k)
             call cicmass(ip)
          enddo
       enddo
!   enddo
    !$omp end parallel do


  write(*,*) '**********'
  write(*,*) '*done cic*'
  write(*,*) '**********'


    return
  end subroutine cic



  subroutine cicmass(ip)
    implicit none
    real, parameter :: ncr=nc

    integer ip

    integer i1,i2,j1,j2,k1,k2
    real x,y,z,dx1,dx2,dy1,dy2,dz1,dz2

    do
       if (ip .eq. 0) exit

       x=modulo(xv(1,ip)-0.5+ncr,ncr)
       y=modulo(xv(2,ip)-0.5+ncr,ncr)
       z=modulo(xv(3,ip)-0.5+ncr,ncr)

       i1=floor(x)+1
       i2=mod(i1,nc)+1
       dx1=i1-x
       dx2=1-dx1
       j1=floor(y)+1
       j2=mod(j1,nc)+1
       dy1=j1-y
       dy2=1-dy1
       k1=floor(z)+1
       k2=mod(k1,nc)+1
       dz1=k1-z
       dz2=1-dz1

       dz1=mp*dz1
       dz2=mp*dz2
       d(i1,j1,k1)=d(i1,j1,k1)+dx1*dy1*dz1
       d(i2,j1,k1)=d(i2,j1,k1)+dx2*dy1*dz1
       d(i1,j2,k1)=d(i1,j2,k1)+dx1*dy2*dz1
       d(i2,j2,k1)=d(i2,j2,k1)+dx2*dy2*dz1
       d(i1,j1,k2)=d(i1,j1,k2)+dx1*dy1*dz2
       d(i2,j1,k2)=d(i2,j1,k2)+dx2*dy1*dz2
       d(i1,j2,k2)=d(i1,j2,k2)+dx1*dy2*dz2
       d(i2,j2,k2)=d(i2,j2,k2)+dx2*dy2*dz2
          
       ip=ll(ip)
    enddo
   
    return
  end subroutine cicmass


!!--------------------------------------------------------------!!

  subroutine powerspectrum
    implicit none

    !real pst(2,nc,nt)
    integer i

#ifdef GAUSS
    real, dimension(nc) :: kgauss
    real, dimension(hc) :: ps_gauss
    character(len=MSL) :: psGaussFile

#endif 

    !write(*,*) 'd before fft', maxval(d),minval(d) 

    !! Convert density to overdensity
    d = d-1

#ifdef RANDOM   
    
    write(*,*) '******************************'
    write(*,*) 'Generating Flat Random Density'
    write(*,*) '******************************'
   
    call random_number(d) 
    d = d-0.5

#endif

    call ps3_r2c(d,ps,nc)

    write(*,*) 'd after fft', maxval(d),minval(d) 
    write(*,*) 'Delta^2(k) after fft', maxval(ps),minval(ps) 

#ifdef GAUSS

    write(*,*) '********************************'
    write(*,*) 'Generating Gaussian Random Field'
    write(*,*) '********************************'

    do i = 1,nc
       kgauss(i) = i*2*pi/lbox
    enddo

#ifdef debug
    write(*,*)'Filled kgauss'
#endif

    call GaussRandomField_3d_r2c(d, lbox, nc, kgauss, ps, nc)

    call ps3_r2c(d,ps_gauss,nc)

    !Get Average PowerSpectrum from 43 simulations from File
    psGaussFile = dir_work//'/PowerSpectrum/z_'//Z//'/psGaussAve_Z'//Z//'.txt'

    write(*,*) 'writing', psGaussFile
    open(unit=25,file=psGaussFile)
 
   
    !***** To write/read the power spectrum to/from a file
    do i=1,hc 
       write(25,*) kgauss(i), ps_gauss(i)
    enddo
    close(25)


#endif



    write(*,*) '********************'
    write(*,*) '*done powerspectrum*'
    write(*,*) '********************'



    return
  end subroutine powerspectrum

!!--------------------------------------------------------------!!

!!--------------------------------------------------------------!!

  subroutine CrossCorrelationOfShells
    implicit none


    integer i,j,k,ii, kx, ky, kz, N_offshell, N_onshell, N_onshell_kx, thread
    !real MAX_Inner, Max_Outer!, MIN_Inner, MIN_Outer, kr, Thickness_Inner, Thickness_Outer, sigma, r0, k0
    real kr, MIN_Inner, MIN_Outer, sigma, r0, k0
    complex(8) psShellCumul_COMPLEX, psShellMean_COMPLEX, N_psShellCumul_COMPLEX, N_psShellMean_COMPLEX

    real, dimension(nc+2,nc,nc) ::  psShell_Inner, psShell_Outer, N_psShell_Inner, N_psShell_Outer
    complex, dimension(nc/2 + 1,nc,nc) ::  psShell_Inner_COMPLEX, psShell_Outer_COMPLEX, N_psShell_Inner_COMPLEX,N_psShell_Outer_COMPLEX
    complex(8), dimension(nc/2 + 1,nc,nc) ::  psShell_COMPLEX, N_psShell_COMPLEX

    equivalence (psShell_Inner,psShell_Inner_COMPLEX)
    equivalence (psShell_Outer,psShell_Outer_COMPLEX)
    equivalence (N_psShell_Inner,N_psShell_Inner_COMPLEX)
    equivalence (N_psShell_Outer,N_psShell_Outer_COMPLEX)

    real, dimension(hc) :: KFromFile
    real, dimension(hc) :: PSFromFile

    character(len=MSL) :: psFile, par_shell!, Max_Inner_char, Max_Outer_char

    !Get Average PowerSpectrum from 43 simulations from File
    !***** To write/read the power spectrum to/from a file

    psFile = dir_work//'/PowerSpectrum/z_'//Z//'/psAve_Z'//Z//'.txt'
    write(*,*) 'opening', psFile
    open(unit=24,file=psFile)
  
    do i=1,hc 
       read(24,*) KFromFile(i), PSFromFile(i)
    enddo
    close(24)

    par_shell = dir_work//'../CubePMLens/par_shell.fh'
    write(*,*) 'opening', par_shell    
    open(unit=25,file=par_shell)
    
    read(25,*) MAX_Inner 
    read(25,*) MAX_Outer 
    read(25,*) MAX_Inner_char 
    read(25,*) MAX_Outer_char

    close(25)
    write(*,*) 'Read',MAX_Inner,MAX_Outer,MAX_Inner_char,MAX_Outer_char 

    ! Initialize variables
    N_offshell = 0 
    N_onshell = 0
    N_onshell_kx = 0
    N_psShellCumul_COMPLEX = 0.
    psShellCumul_COMPLEX = 0.
    thread = 0


    call sfft3_r2c(d,nc,1)
    write(*,*) 'Called First FFT'



    psV_COMPLEX = lbox**3*abs(d_COMPLEX/nc**3)**2
    
#ifdef debug
    write(*,*) 'Got Vectorial Power Spectrum, max Re = ',maxval(real(psV_COMPLEX)), 'min Re= ', minval(real(psV_COMPLEX))
    write(*,*) 'Got Vectorial Power Spectrum, max Im = ',maxval(aimag(psV_COMPLEX)), 'min Im= ', minval(aimag(psV_COMPLEX))
#endif

#ifdef LegendreTest


    !$omp parallel do default(shared) private(i,j,k,kz)
    do k=1,nc

       if(k .lt. hc+2) then
          kz = k-1 ! ->[0,hc]
       else
          kz = k-1-nc ! ->[-hc+1,-1]
       endif
        
       do j=1,nc
          
          do i=1,hc+1

             psV_COMPLEX(i,j,k) = kz**2

          enddo
       enddo

       !write(*,*) 'psV_COMPLEX(4,4,k) =', psV_COMPLEX(4,4,k) 

    enddo
    !$omp end parallel do

    write(*,*) '*******************************'
    write(*,*) 'Overwriting with Legendre Field'
    write(*,*) '*******************************'

    write(*,*) 'Check that these numbers agree with P(k) = k_z^2:'

    write(*,*) 'psV(1,1,1)=', psV_COMPLEX(1,1,1)
    write(*,*) 'psV(2,2,2)=', psV_COMPLEX(2,2,2)
    write(*,*) 'psV(128,128,128)=', psV_COMPLEX(128,128,128)



#endif



    !***************

    ! Define a shell, in fraction of hc (radius of biggest sphere in the k-box)
    ! and set the power spectrum to 0 everywhere else. 
    ! The Shell must be smaller than hc


    write(*,*) '******************************'
    write(*,*) 'Starting Loop over Inner Shell'
    write(*,*) '******************************'



    MIN_Inner = MAX_Inner -  Thickness_Inner


    ! Loop over all k-cells, get magnitude kr  

    !$omp  parallel do default(shared) private(i,j,k,kx,ky,kz,kr)
    !$ thread=omp_get_thread_num() +  1
    !write(*,*) 'thread number =',thread
    do k=1,nc
       if(k .lt. hc+2) then
          kz = k-1 ! ->[0,hc]
       else
          kz = k-1-nc ! ->[-hc+1,-1]
       endif
       do j=1,nc
          if(j .lt. hc+2) then
             ky = j-1 !-> [0,hc]
          else
             ky = j-1-nc ! -> [-hc+1,-1]
          endif

          !do kx=1,hc+1 ! -> [1,hc+1]
          do i=1,hc+1 ! -> [1,hc+1]
             
             kx = i-1 ! -> [0,hc]

             !Spherical Coordinates
             kr= ( (kx)**2 + (ky)**2 + (kz)**2 )**(0.5)

           

             if (kr < MIN_Inner .or. kr > MAX_Inner) then
                N_offshell = N_offshell +1
                psShell_COMPLEX(i,j,k) = 0 ! Set power spectrum to zero
                N_psShell_COMPLEX(i,j,k) = 0
             else
                N_onshell = N_onshell +1 ! Count # cells on shell
                psShell_COMPLEX(i,j,k) = psV_COMPLEX(i,j,k) ! Assign power spectrum
                psShellCumul_COMPLEX = psShellCumul_COMPLEX + psShell_COMPLEX(i,j,k)
                N_psShell_COMPLEX(i,j,k) = 1
                N_psShellCumul_COMPLEX = N_psShellCumul_COMPLEX + N_psShell_COMPLEX(i,j,k)

                if(kx == 0) N_onshell_kx = N_onshell_kx+1

             endif
             
 
          enddo
       enddo
    enddo
    !$omp  end parallel do

    ! **********************************

    ! Compute Mean powerspectrum on the Shell
    psShellMean_COMPLEX = psShellCumul_COMPLEX/N_onshell
    N_psShellMean_COMPLEX = N_psShellCumul_COMPLEX/N_onshell

    ! **********************************


    write(*,*) 'MAX = ', MAX_Inner, 'MIN = ',MIN_Inner

#ifdef debug

    write(*,*) 'N k-cells offshell : ', N_offshell
    write(*,*) 'N k-cells onshell : ', N_onshell
    write(*,*) 'N k-cells onshell_kx : ', N_onshell_kx
    write(*,*) '*******************************'
    write(*,*) 'Total Power on Shell = ', psShellCumul_COMPLEX
    write(*,*) 'Mean Power on Shell = ', psShellMean_COMPLEX
    write(*,*) 'Mean Norm on Shell = ', N_psShellMean_COMPLEX ! Make sure this equals to 1!!!!
    write(*,*) '*******************************'
    write(*,*) 'Delta(k = k-shell) = ',ps(MAX_Inner)
    write(*,*) 'P(k = k-shell-1) = Delta*2*pi**2/k**3 = ',ps(MAX_Inner-1)*2*pi**2/((MAX_Inner-1)*2*pi/lbox)**3
    write(*,*) 'P(k = k-shell) = ',ps(MAX_Inner)*2*pi**2/(MAX_Inner*2*pi/lbox)**3
    write(*,*) 'P(k = k-shell+1) = ',ps(MAX_Inner+1)*2*pi**2/((MAX_Inner+1)*2*pi/lbox)**3
    write(*,*) '<P>(k = k-shell) = ',PSFromFile(MAX_Inner)
    write(*,*) '<P>(k = k-shell-1) = ',PSFromFile(MAX_Inner-1)
    write(*,*) '*******************************'

#endif

    ! **********************************

    !Substract mean from power spectrum on shell:

    write(*,*) 'sum psShell (before mean substraction) =' , sum(psShell_COMPLEX)


    ! Loop again over all k-cells, get magnitude kr  

    !$omp  parallel do default(shared) private(i,j,k,kx,ky,kz,kr)
    do k=1,nc
       if(k .lt. hc+2) then
          kz = k-1 ! ->[0,hc]
       else
          kz = k-1-nc ! ->[-hc+1,-1]
       endif
       do j=1,nc
          if(j .lt. hc+2) then
             ky = j-1 !-> [0,hc]
          else
             ky = j-1-nc ! -> [-hc+1,-1]
          endif

          !do kx=1,hc+1 ! -> [1,hc+1]
          do i=1,hc+1 ! 
             
             kx = i-1 ! -> [0,hc]
             
             !kr=( (kx)**2 + (ky-hc)**2 + (kz-hc)**2)**(0.5)
             kr=( (kx)**2 + (ky)**2 + (kz)**2 )**(0.5)

             if (kr < MIN_Inner .or. kr > MAX_Inner) then
                psShell_COMPLEX(i,j,k) = psShell_COMPLEX(i,j,k)
                N_psShell_COMPLEX(i,j,k) = N_psShell_COMPLEX(i,j,k)
             else
                !Un-comment the next line to include exact mean subtraction
                !psShell_COMPLEX(i,j,k) = psShell_COMPLEX(i,j,k) - psShellMean_COMPLEX                

                !Un-comment the next line to include average over 43 simulations subtraction
                psShell_COMPLEX(i,j,k) = psShell_COMPLEX(i,j,k) - 0.5*(PSFromFile(MAX_Inner)+PSFromFile(MAX_Inner-1))         

                !Un-comment the next line to verify code:  Norm should be exactly == 0
                !N_psShell_COMPLEX(i,j,k) = N_psShell_COMPLEX(i,j,k) - N_psShellMean_COMPLEX 

             endif

          enddo
       enddo
    enddo
    !$omp end parallel do


#ifdef debug
    write(*,*) 'sum psShell (after mean substraction) =' , sum(psShell_COMPLEX)
#endif

    psShell_Inner_COMPLEX = psShell_COMPLEX
    N_psShell_Inner_COMPLEX = N_psShell_COMPLEX

    !****************************


    ! Re-Initialize variables
    N_offshell = 0 
    N_onshell = 0
    N_onshell_kx = 0
    N_psShellCumul_COMPLEX = 0.
    psShellCumul_COMPLEX = 0.



    write(*,*) '******************************'
    write(*,*) 'Starting Loop over Outer Shell'
    write(*,*) '******************************'



    MIN_Outer = MAX_Outer -  Thickness_Outer



    ! Loop over all k-cells, get magnitude kr  

    !$omp  parallel do default(shared) private(i,j,k,kx,ky,kz,kr)
    do k=1,nc
       if(k .lt. hc+2) then
          kz = k-1 ! ->[0,hc]
       else
          kz = k-1-nc ! ->[-hc+1,-1]
       endif
       do j=1,nc
          if(j .lt. hc+2) then
             ky = j-1 !-> [0,hc]
          else
             ky = j-1-nc ! -> [-hc+1,-1]
          endif

          !do kx=1,hc+1 ! -> [1,hc+1]
          do i=1,hc+1
             
             kx = i-1 ! -> [0,hc]

             !Spherical Coordinates
             kr= ( (kx)**2 + (ky)**2 + (kz)**2 )**(0.5)

             

             if (kr < MIN_Outer .or. kr > MAX_Outer) then
                N_offshell = N_offshell +1
                psShell_COMPLEX(i,j,k) = 0 ! Set power spectrum to zero
                N_psShell_COMPLEX(i,j,k) = 0
             else
                N_onshell = N_onshell +1 ! Count # cells on shell
                psShell_COMPLEX(i,j,k) = psV_COMPLEX(i,j,k) ! Assign power spectrum
                psShellCumul_COMPLEX = psShellCumul_COMPLEX + psShell_COMPLEX(i,j,k)
                N_psShell_COMPLEX(i,j,k) = 1
                N_psShellCumul_COMPLEX = N_psShellCumul_COMPLEX + N_psShell_COMPLEX(i,j,k)

                if(kx == 0) N_onshell_kx = N_onshell_kx+1

             endif

 
          enddo
       enddo
    enddo
    !$omp end parallel do

    ! **********************************

    ! Compute Mean powerspectrum on the Shell
    psShellMean_COMPLEX = psShellCumul_COMPLEX/N_onshell
    N_psShellMean_COMPLEX = N_psShellCumul_COMPLEX/N_onshell

    ! **********************************


    write(*,*) 'MAX = ', MAX_Outer, 'MIN = ',MIN_Outer

#ifdef debug

    write(*,*) 'N k-cells offshell : ', N_offshell
    write(*,*) 'N k-cells onshell : ', N_onshell
    write(*,*) 'N k-cells onshell_kx : ', N_onshell_kx
    write(*,*) '*******************************'
    write(*,*) 'Total Power on Shell = ', psShellCumul_COMPLEX
    write(*,*) 'Mean Power on Shell = ', psShellMean_COMPLEX
    write(*,*) 'Mean Norm on Shell = ', N_psShellMean_COMPLEX ! Make sure this equals to 1!!!!
    write(*,*) '*******************************'
    write(*,*) 'Delta(k = k-shell) = ',ps(MAX_Outer)
    write(*,*) 'P(k = k-shell) = Delta*2*pi**2/k**3 = ',ps(MAX_Outer)*2*pi**2/(MAX_Inner*2*pi/lbox)**3
    write(*,*) '*******************************'

#endif

    ! **********************************

    !Substract mean from power spectrum on shell:

    write(*,*) 'sum psShell (before mean substraction) =' , sum(psShell_COMPLEX)


    ! Loop again over all k-cells, get magnitude kr  

    !$omp  parallel do default(shared) private(i,j,k,kx,ky,kz,kr)
    do k=1,nc
       if(k .lt. hc+2) then
          kz = k-1 ! ->[0,hc]
       else
          kz = k-1-nc ! ->[-hc+1,-1]
       endif
       do j=1,nc
          if(j .lt. hc+2) then
             ky = j-1 !-> [0,hc]
          else
             ky = j-1-nc ! -> [-hc+1,-1]
          endif

          !do kx=1,hc+1 ! -> [1,hc+1]
          do i=1,hc+1
             
             kx = i-1 ! -> [0,hc]
             
             !kr=( (kx)**2 + (ky-hc)**2 + (kz-hc)**2)**(0.5)
             kr=( (kx)**2 + (ky)**2 + (kz)**2 )**(0.5)

             if (kr < MIN_Outer .or. kr > MAX_Outer) then
                psShell_COMPLEX(i,j,k) = psShell_COMPLEX(i,j,k)
                N_psShell_COMPLEX(i,j,k) = N_psShell_COMPLEX(i,j,k)
             else
                !Un-comment the next line to include mean subtraction
                !psShell_COMPLEX(i,j,k) = psShell_COMPLEX(i,j,k) - psShellMean_COMPLEX

                !Un-comment the next line to include average over 43 simulations subtraction
                psShell_COMPLEX(i,j,k) = psShell_COMPLEX(i,j,k) - 0.5*(PSFromFile(MAX_Outer)+PSFromFile(MAX_Outer-1)) 

                !Un-comment the next line to verify code:  Norm should be exactly == 0
                !N_psShell_COMPLEX(i,j,k) = N_psShell_COMPLEX(i,j,k) - N_psShellMean_COMPLEX 

             endif

          enddo
       enddo
    enddo
    !$omp end parallel do

#ifdef debug
    write(*,*) 'sum psShell (after mean substraction) =' , sum(psShell_COMPLEX)
#endif

    psShell_Outer_COMPLEX = psShell_COMPLEX
    N_psShell_Outer_COMPLEX = N_psShell_COMPLEX




    !******************************************************
    write(*,*) '*** Cross-Correlation of Normalization ***'
    !******************************************************

    call GetCrossCorr(N_psShell_Inner, N_psShell_Outer, N_CrossCorr, nc)

#ifdef debug
    write(*,*) 'Vectorial Cross-Correlation Re: ', maxval(real(N_CrossCorr_COMPLEX)),' > N_CrossCorr > ',minval(real(N_CrossCorr_COMPLEX))
    write(*,*) 'Vectorial Cross-Correlation Im: ', maxval(aimag(N_CrossCorr_COMPLEX)),' > N_CrossCorr > ',minval(aimag(N_CrossCorr_COMPLEX))
#endif

    call Get_KAverage(N_CrossCorr,N_CrossCorrAve,nc,1.,1.)

#ifdef debug
    write(*,*) 'Averaged Cross-Correlation', maxval(N_CrossCorrAve),' > CrossCorrAve > ',minval(N_CrossCorrAve) 
#endif



    !**********************************************
    write(*,*) '*** Cross-Correlation of Shell ***'
    !**********************************************

#ifdef debug
    write(*,*) 'sum psShell_Inner (after mean substraction) =' , sum(psShell_Inner_COMPLEX)
    write(*,*) 'sum psShell_Outer (after mean substraction) =' , sum(psShell_Outer_COMPLEX)
#endif

    call GetCrossCorr(psShell_Inner, psShell_Outer, CrossCorr, nc)

#ifdef debug
    write(*,*) 'Vectorial Cross-Correlation', maxval(CrossCorr),' > CrossCorr > ',minval(CrossCorr)
    write(*,*) 'Vectorial Cross-Correlation Re: ', maxval(real(CrossCorr_COMPLEX)),' > CrossCorr > ',minval(real(CrossCorr_COMPLEX))
    write(*,*) 'Vectorial Cross-Correlation Im: ', maxval(aimag(CrossCorr_COMPLEX)),' > CrossCorr > ',minval(aimag(CrossCorr_COMPLEX))

    write(*,*) 'TEST : SUM over CrossCorrelation Re = ', sum(real(CrossCorr_COMPLEX))
    write(*,*) 'TEST : SUM over CrossCorrelation Im = ', sum(aimag(CrossCorr_COMPLEX))
#endif


!! The following prints the Values of CrossCorr above a certain threshold, along with the coordinates
!! Typically, these appear 

!    do k=1,nc
!       if(k .lt. hc+2) then
!          kz = k-1 ! ->[0,hc]
!       else
!          kz = k-1-nc ! ->[-hc+1,-1]
!       endif
!       do j=1,nc
!          if(j .lt. hc+2) then
!             ky = j-1 !-> [0,hc]
!          else
!             ky = j-1-nc ! -> [-hc+1,-1]
!          endif
!          
!          !do kx=1,hc+1 ! -> [1,hc+1]
!          do i=1,hc+1
!             
!             kx = i-1 ! -> [0,hc]
!
!             kr=( (kx)**2 + (ky)**2 + (kz)**2 )**(0.5)
!             
!             if(real(CrossCorr_COMPLEX(i,j,k)) > 3) then
!                write(*,*) 'k = (',kx,ky,kz,'), kr = ',kr  
!                write(*,*) 'CrossCorr(k) = ', CrossCorr_COMPLEX(i,j,k)
!                pause
!             endif
!
!          enddo
!       enddo
!    enddo


    call Get_KAverage(CrossCorr,CrossCorrAve,nc,1.,1.)

#ifdef debug
    write(*,*) 'Averaged Cross-Correlation', maxval(CrossCorrAve),' > CrossCorrAve > ',minval(CrossCorrAve) 
#endif


 



    write(*,*) '***********************'
    write(*,*) '*done CrossCorrelation*'
    write(*,*) '***********************'



    return
  end subroutine CrossCorrelationOfShells

!*********************************

!! 1- FFT-Inverse a Power Spectrum to.get the Correlation Function
!! 2- AbsoluteValue-Square the Correlation Function
!! 3- FFT to get a Cross-Correlation(delta_k)

subroutine GetCrossCorr(ps1,ps2,Cross,nc)
!subroutine GetCrossCorr(ps1_COMPLEX,ps2_COMPLEX,Cross_COMPLEX,nc)
  implicit none

  integer nc
  real, dimension(nc+2,nc,nc) :: ps1, ps2, Cross  

#ifdef debug
  write(*,*) 'P1:sum over density in Fourier Space =' , sum(ps1),'nc = ', nc
#endif

  call sfft3_r2c(ps1,nc,-1) 

#ifdef debug
  write(*,*) 'Xi_1(0,0,0) = ', ps1(1,1,1)   ! NOTE: Xi(0,0,0) ~ 2*sum/nc**3 (this double counts the x=0 slab)
  write(*,*) 'P2:sum over density in Fourier Space =' , sum(ps2),'nc = ', nc
#endif

  call sfft3_r2c(ps2,nc,-1) 
#ifdef debug
  write(*,*) 'Xi_2(0,0,0) = ', ps2(1,1,1)   ! NOTE: Xi(0,0,0) ~ 2*sum/nc**3 (this double counts the x=0 slab)
#.le.endif


  !! Now ps1 and ps2 are a "real" function, going from 1 to nc 
  !! in all 3 dimensions. 

  Cross = ps1*ps2


#ifdef debug
  write(*,*) 'Xi_1(1,0,0) = ', ps1(2,1,1)
  write(*,*) 'Xi_1(0,1,0) = ', ps1(1,2,1)
  write(*,*) 'Xi_1(0,0,1) = ', ps1(1,1,2)
  write(*,*) 'Xi_1(10,0,0) = ', ps1(11,1,1)
  write(*,*) 'Xi_1(0,10,0) = ', ps1(1,11,1)
  write(*,*) 'Xi_1(0,0,10) = ', ps1(1,1,11)

  write(*,*) 'Xi_1(10,10,10) = ', ps1(11,11,11)
  write(*,*) 'Xi_1(11,10,10) = ', ps1(12,11,11)
  write(*,*) 'Xi_1(10,11,10) = ', ps1(11,12,11)
  write(*,*) 'Xi_1(10,10,11) = ', ps1(11,11,12)
  write(*,*) 'Xi_1(21,20,20) = ', ps1(22,21,21)
  write(*,*) 'Xi_1(20,21,20) = ', ps1(21,22,21)
  write(*,*) 'Xi_1(20,20,21) = ', ps1(21,21,22)

  write(*,*) 'Xi_1*Xi_2(0,0,0) = ', Cross(1,1,1)
  write(*,*) 'Xi_1*Xi_2(1,0,0) = ', Cross(2,1,1)
  write(*,*) 'Xi_1*Xi_2(0,1,0) = ', Cross(1,2,1)
  write(*,*) 'Xi_1*Xi_2(0,0,1) = ', Cross(1,1,2)
  write(*,*) 'Xi_1*Xi_2(10,0,0) = ', Cross(11,1,1)
  write(*,*) 'Xi_1*Xi_2(0,10,0) = ', Cross(1,11,1)
  write(*,*) 'Xi_1*Xi_2(0,0,10) = ', Cross(1,1,11)

  write(*,*) 'Cross before FFT : ', maxval(Cross),' > Cross > ' ,minval(Cross)
#endif
 
  call sfft3_r2c(Cross,nc,1) 

  !Now Cross is in k-space. This subroutine defined Cross as a real, 
  !but even x-array elements represent imaginary parts


#ifdef debug
  write(*,*) 'Cross after FFT : ', maxval(Cross),' > Cross > ' ,minval(Cross)
#endif


  return
end subroutine GetCrossCorr

!*******************************


!======================================================
 Subroutine Get_KAverage(map,Ave,n,normp,normk)
!======================================================
! Calculates the Power Spectrum and Dumps it into a file


   Integer                            :: n
   Real                               :: normp,normk
   Real, Dimension(1:n+2,1:n,1:n)  :: map
   Character(Len=120)              :: file

   !Real*4, Dimension(1:n+2,1:n,1:n):: fft
   Real, Dimension(3,1:n)        :: pst
   Real, Dimension(n)        :: Ave

   Real  :: w1,w2,kz,kx,ky,kr,pow
   Integer :: i,j,k,hn,k1,k2

   If (mod(n,2) .NE. 0) Then
      Write(*,*) ' n must be even'
      Stop
   End If


   !! == Dump power spectra
   !Open(unit=40,file=file,status='replace')
   hn  = n/2
   pst = 0.0

   !$omp  parallel do default(shared) private(i,j,k,kx,ky,kz,kr, k1,k2,pow,w1,w2)
   Do k = 1,n
      If (k .Lt. hn+2) Then
         kz = k-1
      Else
         kz = k-1-n
      Endif
      Do j = 1,n
         If (j .Lt. hn+2) Then
            ky = j-1
         Else
            ky = j-1-n
         Endif
         Do i = 1,n+2,2
            kx = (i-1)/2
            kr = Sqrt(kx**2+ky**2+kz**2)

            If (kr .Ne. 0.) Then
               k1  = Ceiling(kr)
               k2  = k1+1
               w1  = k1-kr
               w2  = 1-w1
               pow = map(i,j,k)
               pst(1,k1)=pst(1,k1)+w1*pow
               pst(2,k1)=pst(2,k1)+w1*pow**2
               pst(3,k1)=pst(3,k1)+w1 ! Count the number of elements
               pst(1,k2)=pst(1,k2)+w2*pow
               pst(2,k2)=pst(2,k2)+w2*pow**2
               pst(3,k2)=pst(3,k2)+w2 ! Count the number of elements
               !if (kr<2) write(*,*) 'kr =',kr, 'k1 =',k1,'k2=',k2
               !if (kr<2) write(*,*) 'w1 =',w1,'kw=',w2, 'pst(k1) =',pst(1,k1),'pst(k2) =',pst(1,k2) 
            Else
               !write(*,*) 'kr =',kr, 'pst(kr)=', pst(1,1) 
            Endif
         Enddo
      Enddo
   End Do
   !$omp end parallel do

   Ave(1) = map(1,1,1)
   Do k = 2,n
      If (pst(3,k) .Eq. 0) Then
         Ave(k) = 0
      Else 
         Ave(k) = pst(1,k)/pst(3,k)
      Endif
   Enddo
   !write(*,*) 'Ave(1) =', Ave(1), 'Ave(2) =', Ave(2)
   !Close(40)
   !Write(*,*) '3D power spectra written in ',Trim(file)
   .le.Write(*,*) '3D Average Done '

   return

End Subroutine Get_KAverage


!*********************************************


!*********************************

! Usefull after a Fourier Transform the density, this routine rearranges indices to run from [1,hc+1],[1,nc],[1,nc]
! Corresponding to physical wave nubers of [1,hc+1],[-hc+1,hc] and [-hc+1,hc]

! d must be in 3-D Fourrier space
 
! NOTE : the density will have Real/Imaginary components arranged as follow:
! Re(rho) from kx -> [1,hc+1], then Im(rho) from kx -> [hc+2,nc+2], 

! Ex: After Rearrange(d,nc), a slice through the kx axis (ky = kz = 0) is given by:
! do k = 1:hc
!   d(k,hc,hc)
! enddo

subroutine Rearrange(d,nc)
  implicit none
  
  integer nc,hc,FB
  integer, parameter :: nt=4
  
  real, dimension(nc+2,nc,nc) :: d
  real, dimension(nc+2,nc,nc, nt) :: rhot
  
  integer it,i,j,k,kpt
  real kx,ky,kz,pi
  


  pi=acos(-1.)
  hc=nc/2
  kpt=nc/nt

  !$omp parallel do default(shared) &
  !$omp& private(it,i,j,k,kr,kx,ky,kz,k1,k2,w1,w2,pow)
  do it=1,nt
     rhot(:,:,:,it)=0
     !weightt(:,it)=0
     do k=1+(it-1)*kpt,min(nc,it*kpt)
        if (k .lt. hc+2) then
            kz=k-1 ! kz -> [0,hc]
        else
            kz=k-1-nc ! kz -> [-hc+1, -1] 
        endif

        ! Now, kz runs in the range [-hc+1,hc]
        ! I will shift values of kz so that they run from [1,nc]
        kz = kz + hc
        
        do j=1,nc
           if (j .lt. hc+2) then
              ky=j-1 ! ky -> [0,hc]
           else
              ky=j-1-nc ! ky -> [-hc+1, -1] 
           endif

           ! Now, ky runs in the range [-hc+1,hc]
           ! I will shift values of ky so that they run from [1,nc]
           ky = ky + hc


          ! kx is particular: (i odd = real, i even = imaginary)
          ! so I will bring Re(rho) from kx -> [1,hc+1], then Im(rho) from kx -> [hc+2,nc+2]


           do i=1,nc+1,2 
              kx=(i+1)/2 

              ! kx runs from [1,hc+1]. 

              ! Filling Re(rho) from [1,hc+1]
              rhot(kx,ky,kz,it)=d(i,j,k)

              ! Filling Im(rho) from [hc+2,nc+2]
              rhot(hc+kx+1,ky,kz,it)=d(i+1,j,k)
                     
              !write(*,*) 'Test : kx = ', kx, 'ky = ', ky, 'kz = ',kz
              !pause


           enddo
        enddo
     enddo
  enddo
  !$omp end parallel do
  

  !! Merge density from threads

  d=0
  do it=1,nt  
     d=d+rhot(:,:,:,it)
  enddo

  !write(*,*) 'Density : ', maxval(d), ' > d > ',minval(d)
  
 
  return
end subroutine Rearrange

!*****************************************



!! 1- FFT-Inverse a Power Spectrum to.get the Correlation Function
!! 2- AbsoluteValue-Square the Correlation Function
!! 3- FFT to get a Covariance(delta_k)

subroutine GetCov(PS,Cov,nc)
  implicit none

  integer nc
  real, dimension(nc+2,nc,nc) :: PS, Cov  
  integer i,j,k

  write(*,*) 'PS before backward FFT : ', maxval(PS),' > PS > ' ,minval(PS)

  !do ii = 1,nc
     !write(*,*) 'PPShell(', ii, ', hc, hc) = ', d(ii,hc,hc)  
  !enddo

  call sfft3_r2c(PS,nc,-1) 


  !! Now PS is a "mostly-real" function, going from 1 to nc 
  !! in all 3 dimensions.Let's absolute-Square the power spectrum

  Cov = abs(PS)**2

#ifdef debug
  write(*,*) 'PS after backward FFT : ', maxval(PS),' > PS > ' ,minval(PS)

  write(*,*) 'Xi^2(0,0,0) = ', Cov(1,1,1)
  write(*,*) 'Xi^2(1,0,0) = ', Cov(2,1,1)
  write(*,*) 'Xi^2(0,1,0) = ', Cov(1,2,1)
  write(*,*) 'Xi^2(0,0,1) = ', Cov(1,1,2)

  write(*,*) 'Cov before FFT : ', maxval(Cov),' > Cov > ' ,minval(Cov)
#endif

  call sfft3_r2c(Cov,nc,1) 

#.le.ifdef debug
  write(*,*) 'Cov after FFT : ', maxval(Cov),' > Cov > ' ,minval(Cov)
#endif

  !Now Cov is in k-space, orderer as the fftw wishes it.

  return
end subroutine GetCov

!*******************************




end program cicpow
