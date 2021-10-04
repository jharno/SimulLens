! To be compiled and run on lobster6, otherwise the calshear part will
! dump core and stop.
!

program FullCMBLens
  use StringOpt
  use Lensing
  implicit none
#include 'CMBLens.fh'


  real, dimension(nc,nc) :: rho_pxy,rho_pxz,rho_pyz
  real, dimension(nc,nc,nslice) :: map1
  real, dimension(nc,nc) :: cumul_kappa,cumul_gamma1,cumul_gamma2,cumul_deflx,cumul_defly
  real, dimension(nslice) :: z_write,box
  real, dimension(nc,nc) :: phi
  type(vec2D), dimension(nc,nc) :: defl
  type(vec2D), dimension(nc,nc,nslice) :: newdefl
  type(Spin2), dimension(nc,nc) :: shear
  type(Spin2), dimension(nc,nc,nslice) :: newshear
  type(vec2D), dimension(nc,nc) :: CorrBornDefl

  type(vec2D) shift
  integer i,j,fu,i1,j1,icount,kx,ky,ir,index_nbody_run,file_status,ScaleFactor!,redshift, z
  real lense_weight,frac,ncr,angle,dang,chi,pi,random_index_run!,box
  real(kind=8) rhomean
  character(len=180) :: fn,fn1,fp
  character(len=7) z_string

  box=lbox
  ncr=nc
  pi=acos(-1.)

  write(*,*) '*****************'
  write(*,*) 'Running SimulLens'
  write(*,*) '*****************'

  write(*,*) 'nc :',nc
  write(*,*) 'npc = nc by construction'
  write(*,*) 'source redshift:',zs

 !Field of view and resolution:
  !angle=sqrt(12.84)/180*pi
  !angle=box/chi(zs,omegam,h)/h
  !dang=angle/nc

  !write(*,*) 'Field of view :',(angle*180./pi)**2,'sq.degrees'
  !write(*,*) 'Resolution    :',(dang*180./pi)*60,'arc min.'
  open(11,file=fn_z)
  do i=1,nslice
     read(11,*) z_write(i)
  enddo
  close(11)
  open(11,file=fn_box)
  do i=1,nslice
     read(11,*) box(i)
  enddo
  close(11)


  do ir=1,nr
     write(*,*) ir
     call random_seed()
     !call random_number(shift%x)
     !call random_number(shift%y)

     icount=0

     !Loop over the slices
     do j=1,nslice
        i=mod(icount,3)

        ! Get the Redshift
        write(z_string,'(f7.3)') z_write(j)
        z_string=adjustl(z_string)

        ! for mixed size of simulations
        !if(z_write(j) > 1.0) box=sbox !!!! To modify... I do not need this in fact...
        !if(z_write(j) <= 1.0) box=lbox


        ! Get the projection
        !if(z_write(j) <=1.0)fn=proj_path//z_string(1:len_trim(z_string)) ! simple structure                
        !if(z_write(j) >1.0)fn=proj_path2//z_string(1:len_trim(z_string)) ! simple structure        
        fn=proj_path//z_string(1:len_trim(z_string)) 

        fn=adjustl(fn)
        if(i.eq.0)fp=fn(1:len_trim(fn))//'proj_xy.dat.bin'
        if(i.eq.1)fp=fn(1:len_trim(fn))//'proj_yz.dat.bin'
        if(i.eq.2)fp=fn(1:len_trim(fn))//'proj_xz.dat.bin'
        fu=10+i

        open(unit=fu,file=fp,form='binary',status='old',iostat=file_status)
        write(*,*) 'Opened', fp
        !if(file_status.gt.0)write(*,*) 'file_status',file_status


        !Read the files
        read(fu) ScaleFactor
        read(fu) rho_pxy
        !write(*,*) 'Read', fp
        close(fu)
        !pause
        map1(:,:,j)=rho_pxy
        write(*,*) 'Got projection!'

#ifdef test_read
        write (fn,'(f5.3,"test.dat")') z_write(j)
        open(10,file=Lens_output_path//"/"//fn,form='binary')
        write(10) map1(:,:,j)
        close(10)
        write(*,*)'Wrote test Map for z=',z_write(j)
        stop
#endif


        !call random_number(shift%x)
        !call random_number(shift%y)

        !! substract the mean to get the overdensity
        rhomean=sum(real(map1(:,:,j),kind=8))/nc/nc
        map1(:,:,j)=map1(:,:,j)-rhomean !nc 
 
        call kappa_to_defl(map1(:,:,j),defl,shear,phi,nc) !Using discrete differentiation
        !call kappa2gamma(map1(:,:,j),defl,shear,phi,nc)  !Using continuous differentiation 
        write(*,*) 'Calculated shear and deflection'

        !frac=angle/(box / chi(z_write(j),omegam,h) / h)

        !call zoomshiftmap(map1(:,:,j),map1(:,:,j),nc,nc,shift,newdefl,frac)
        !write(*,*) 'Randomly Shifted the delta maps'

        if(z_write(j)<10) write (fn,'(f5.3,"delta.dat")') z_write(j)
        if(z_write(j)>10 .and. z_write(j)<100) write (fn,'(f6.3,"delta.dat")') z_write(j)
        if(z_write(j)>100)write (fn,'(f7.3,"delta.dat")') z_write(j)
        !open(10,file=Lens_output_path//"run_"//digitn(ir,2)//"/"//fn,form='binary')
        open(10,file=Lens_output_path//"/"//fn,form='binary')
        write(10) map1(:,:,j)
        close(10)
        write(*,*)'Wrote kappa Map for z=',z_write(j)


        !call zoomshiftmap(shear%p,newshear(:,:,j)%p,nc,nc,shift,CorrBornDefl,frac)
        !call zoomshiftmap(shear%c,newshear(:,:,j)%c,nc,nc,shift,CorrBornDefl,frac)
        !call zoomshiftmap(defl%x,newdefl(:,:,j)%x,nc,nc,shift,CorrBornDefl,frac)
        !call zoomshiftmap(defl%y,newdefl(:,:,j)%y,nc,nc,shift,CorrBornDefl,frac)


        ! multiply the deflection angle by the basic units d_theta_lens
        write(*,*) 'Theta factor = ', (box(j)/nc)/chi(z_write(j),omegam,h)/h
        newdefl(:,:,j)%x=defl(:,:)%x*(box(j)/nc)/chi(z_write(j),omegam,h)/h
        newdefl(:,:,j)%y=defl(:,:)%y*(box(j)/nc)/chi(z_write(j),omegam,h)/h
        newshear(:,:,j)%p = shear(:,:)%p
        newshear(:,:,j)%c = shear(:,:)%c
        !write(*,*) 'Randomly Shifted the shear and deflection maps'

        icount=icount+1
     enddo !!j=1,nslice


     !Start the integral over redshift slices 

     do i=nslice,1,-1

        ! Multiply by weight and geometry, intgrate up to the lens
        ! redshift
        !write(*,*) 'Integrating up to z=',z_write(i)

        cumul_kappa(:,:)=0.0
        cumul_gamma1(:,:)=0.0
        cumul_gamma2(:,:)=0.0
        cumul_deflx(:,:)=0.0
        cumul_defly(:,:)=0.0

        do j = nslice,i,-1

           ! Get the good box size: Those of the lenses we sum over!
           !if(z_write(j) > 1.0) box=sbox
           !if(z_write(j) <= 1.0) box=lbox

        !   write(*,*) 'Summing contribution of z=',z_write(j)                 
           cumul_kappa(:,:) =cumul_kappa(:,:)+map1(:,:,j)*(3./2)*omegam*(box(j)/nc)/(3.E3)**2*(chi(z_write(j),omegam,h)*h)*(1+z_write(j))*(1-chi(z_write(j),omegam,h)/chi(z_write(i),omegam,h))
           cumul_gamma1(:,:)=cumul_gamma1(:,:)+newshear(:,:,j)%p*(3./2)*omegam*(box(j)/nc)/(3.E3)**2*(chi(z_write(j),omegam,h)*h)*(1+z_write(j))*(1-chi(z_write(j),omegam,h)/chi(z_write(i),omegam,h))
           cumul_gamma2(:,:)=cumul_gamma2(:,:)+newshear(:,:,j)%c*(3./2)*omegam*(box(j)/nc)/(3.E3)**2*(chi(z_write(j),omegam,h)*h)*(1+z_write(j))*(1-chi(z_write(j),omegam,h)/chi(z_write(i),omegam,h))
           cumul_deflx(:,:)=cumul_deflx(:,:)+newdefl(:,:,j)%x*(3./2)*omegam*(box(j)/nc)/(3.E3)**2*(chi(z_write(j),omegam,h)*h)*(1+z_write(j))*(1-chi(z_write(j),omegam,h)/chi(z_write(i),omegam,h))
           cumul_defly(:,:)=cumul_defly(:,:)+newdefl(:,:,j)%y*(3./2)*omegam*(box(j)/nc)/(3.E3)**2*(chi(z_write(j),omegam,h)*h)*(1+z_write(j))*(1-chi(z_write(j),omegam,h)/chi(z_write(i),omegam,h))
        enddo

        !Write to files


        if(z_write(i)<10) write (fn,'(f5.3,"kappa_weight.dat")') z_write(i)
        if(z_write(i)>10 .and. z_write(i)<100) write (fn,'(f6.3,"kappa_weight.dat")') z_write(i)
        if(z_write(i)>100)write (fn,'(f7.3,"kappa_weight.dat")') z_write(i)
        !open(10,file=Lens_output_path//"run_"//digitn(ir,2)//"/"//fn,form='binary')
        open(10,file=Lens_output_path//"/"//fn,form='binary')
        write(10) cumul_kappa(:,:)
        close(10)
        write(*,*)'Wrote kappa Map for z=',z_write(i)

        if(z_write(i)<10) write (fn,'(f5.3,"gamma1_weight.dat")') z_write(i)
        if(z_write(i)>10 .and. z_write(i)<100) write (fn,'(f6.3,"gamma1_weight.dat")') z_write(i)
        if(z_write(i)>100)write (fn,'(f7.3,"gamma1_weight.dat")') z_write(i)
        !open(10,file=Lens_output_path//"run_"//digitn(ir,2)//"/"//fn,form='binary')
        open(10,file=Lens_output_path//"/"//fn,form='binary')
        write(10) cumul_gamma1(:,:)
        close(10)
        if(z_write(i)<10) write (fn,'(f5.3,"gamma2_weight.dat")') z_write(i)
        if(z_write(i)>10 .and. z_write(i)<100) write (fn,'(f6.3,"gamma2_weight.dat")') z_write(i)
        if(z_write(i)>100)write (fn,'(f7.3,"gamma2_weight.dat")') z_write(i)
        !open(10,file=Lens_output_path//"run_"//digitn(ir,2)//"/"//fn,form='binary')\
        open(10,file=Lens_output_path//"/"//fn,form='binary')
        write(10) cumul_gamma2(:,:)
        close(10)
        write(*,*)'Wrote shear Maps for z=',z_write(i)

        if(z_write(i)<10) write (fn,'(f5.3,"deflx_weight.dat")') z_write(i)
        if(z_write(i)>10 .and. z_write(i)<100) write (fn,'(f6.3,"deflx_weight.dat")') z_write(i)
        if(z_write(i)>100)write (fn,'(f7.3,"deflx_weight.dat")') z_write(i)
       !open(10,file=Lens_output_path//"run_"//digitn(ir,2)//"/"//fn,form='binary')
        open(10,file=Lens_output_path//"/"//fn,form='binary')
        write(10) cumul_deflx(:,:)
        close(10)
        if(z_write(i)<10) write (fn,'(f5.3,"defly_weight.dat")') z_write(i)
        if(z_write(i)>10 .and. z_write(i)<100) write (fn,'(f6.3,"defly_weight.dat")') z_write(i)
        if(z_write(i)>100)write (fn,'(f7.3,"defly_weight.dat")') z_write(i)
        !open(10,file=Lens_output_path//"run_"//digitn(ir,2)//"/"//fn,form='binary')\
        open(10,file=Lens_output_path//"/"//fn,form='binary')
        write(10) cumul_defly(:,:)
        close(10)
        write(*,*)'Wrote defl Maps for z=',z_write(i)

     enddo

     write(*,*) 'Done integral over redshift slices'

     open(20,file=Lens_output_path//"/"//fn_kappa,form='binary')
     write(20) cumul_kappa(:,:)
     close(20)
     write(*,*)'Wrote integrated kappa Maps'

     open(20,file=Lens_output_path//"/"//fn_gamma1,form='binary')
     write(20) cumul_gamma1(:,:)
     close(20)
     !finalmap=0
     !do i=1,nslice
     !   if(z_write(i).le.zs)finalmap(:,:)=finalmap(:,:)+newshear(:,:,i)%c
     !enddo
     !open(20,file=Lens_output_path//fn_gamma2,form='binary')
     open(20,file=Lens_output_path//"/"//fn_gamma2,form='binary')
     write(20) cumul_gamma2(:,:)
     close(20)
     !open(20,file=Lens_output_path//fn_defl,form='binary')
     open(20,file=Lens_output_path//"/"//fn_defl,form='binary')
     write(20)cumul_deflx(:,:),cumul_defly(:,:)
     close(20)

     write(*,*)'Wrote Integrated shear Maps'

  enddo !! ir-1,nr

end program FullCMBLens 

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

    !write(*,*) 'Calling zoomshiftmap'

    do j1=1,nc
       do i1=1,nc
          map3(i1,j1)=map1(i1,j1)
          map3(i1+nc,j1)=map1(i1,j1)
          map3(i1,j1+nc)=map1(i1,j1)
          map3(i1+nc,j1+nc)=map1(i1,j1)
       enddo
    enddo

    call zoom(map3(:,:),map2(:,:),nc,npc,shift,defl,frac)

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






