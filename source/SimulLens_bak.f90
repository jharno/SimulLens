
program SimulLens
  use StringOperation
  use Lensing
  implicit none
#include 'par.fh'

  real*8, dimension(nc,nc) :: rho_pxy,rho_pxz,rho_pyz
  real, dimension(npc,npc,nslice) :: map
  real, dimension(nc,nc,nslice) :: map1
  real, dimension(npc,npc) :: finalmap 
  real, dimension(npc) :: power
  real, dimension(nslice) :: z_write
  real, dimension(nc,nc) :: phi
  type(vec2D), dimension(nc,nc) :: defl
  type(vec2D), dimension(npc,npc,nslice) :: newdefl
  type(vec2D), dimension(npc,npc) :: CorrBornDefl
  type(Spin2), dimension(nc,nc) :: shear  
  type(Spin2), dimension(npc,npc,nslice) :: newshear
  

  type(vec2D) shift
  integer i,j,fu,i1,j1,icount,kx,ky,ir
  real redshift,lense_weight,z,frac,ncr,angle,dang,chi,pi
  real(kind=8) rhomean
  character(len=180) :: fn,fn1,fp

  ncr=nc
  pi=acos(-1.)
  angle=lbox/chi(zs,omegam,h)/h
  dang=angle/npc
  open(11,file=fn_z)  
  do i=1,nslice
     !read(unit=11,fmt='(f7.3)') z_write(i)
     read(11,*) z_write(i)
  enddo
  close(11)

  do ir=1,nr
     write(*,*) ir
     call random_seed()
     icount=0
   
     do j=1,nslice
        i=mod(icount,3)
        if (i.eq.0) write (fn,'(f7.3,"proj_xy.dat")') z_write(j)
        if (i.eq.1) write (fn,'(f7.3,"proj_yz.dat")') z_write(j)
        if (i.eq.2) write (fn,'(f7.3,"proj_xz.dat")') z_write(j)
        fn=adjustl(fn)
        write(*,*) fn
        fp=PROJPATH//fn
        fu=10+i
        open(unit=fu,file=fp,form='binary')
        read(fu) redshift
        read(fu) rho_pxy
        close(fu)
        map1(:,:,j)=rho_pxy

        !! overdensity from PMFAST, density from CUBEPM
        rhomean=sum(real(map1(:,:,j),kind=8))/nc/nc 
        map1(:,:,j)=map1(:,:,j)-nc 
        lense_weight=(chi(z_write(j),omegam,h)*h)*(1+z_write(j))*(1-chi(z_write(j),omegam,h)/chi(zs,omegam,h))

        map1(:,:,j)=map1(:,:,j)*lense_weight*lbox/nc/(3.E3)**2
        map1(:,:,j)=3./2*omegam*map1(:,:,j)

        call kappa_to_defl(map1(:,:,j),defl,shear,phi,nc)

        !! zoom the slices to fit the field of view
        frac=angle*chi(z_write(j),omegam,h)*h/lbox 
        call random_number(shift%x)
        call random_number(shift%y)

        call zoomshiftmap(map1(:,:,j),map(:,:,j),nc,npc,shift,newdefl,frac)

        call zoomshiftmap(shear%p,newshear(:,:,j)%p,nc,npc,shift,CorrBornDefl,frac)
        call zoomshiftmap(shear%c,newshear(:,:,j)%c,nc,npc,shift,CorrBornDefl,frac)
        call zoomshiftmap(defl%x,newdefl(:,:,j)%x,nc,npc,shift,CorrBornDefl,frac)
        call zoomshiftmap(defl%y,newdefl(:,:,j)%y,nc,npc,shift,CorrBornDefl,frac)
        !newdefl%x=newdefl%x*(lbox/nc)/chi(z_write(j),omegam,h)/h
        !newdefl%y=newdefl%y*(lbox/nc)/chi(z_write(j),omegam,h)/h
        newdefl(:,:,j)%x=newdefl(:,:,j)%x*(lbox/nc)/chi(z_write(j),omegam,h)/h
        newdefl(:,:,j)%y=newdefl(:,:,j)%y*(lbox/nc)/chi(z_write(j),omegam,h)/h

        icount=icount+1

#ifdef SemiBorn
        CorrBornDefl(:,:)%x=newdefl(:,:,j)%x
        CorrBornDefl(:,:)%y=newdefl(:,:,j)%y
#endif
     enddo

     finalmap=0
     do i=1,nslice
        finalmap(:,:)=finalmap(:,:)+map(:,:,i)
     enddo

     open(20,file=dir_work//fn_kappa,form='binary')
     write(20)finalmap
     close(20)

     open(20,file=dir_work//fn_defl,form='binary')
     write(20)sum(newdefl(:,:,:)%x,dim=3),sum(newdefl(:,:,:)%y,dim=3)
     close(20)

     !! lense_weight=distance(z)*(1+z)*(1-chi(z(i))/chi(zs))
     !! Note the unit of chi is Mpc, not h^-1*Mpc
     !! 3.E3=H0/c=100/3*e5 (h Mpc^-1)

     call ps2(finalmap,finalmap,power,npc)
     open(30,file=dir_work//'l2cl_kappa.'//digitn(ir,2))
     do i=1,npc/2
        write(30,*) i*2*pi/angle,power(i)
     enddo
     close(30)

  enddo

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


