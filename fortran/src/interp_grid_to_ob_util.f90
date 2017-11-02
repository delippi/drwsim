module interp_util


!******************************************************************
!    Below are various global variables needed throughout
!******************************************************************

  use kinds, only: r_kind,i_kind,i_byte
  implicit none
  
  real(r_kind),parameter:: r1_5=1.5_r_kind
  integer(i_kind),parameter :: nlon=1728_i_kind	
  integer(i_kind),parameter :: nlat=1440_i_kind
  integer(i_kind),parameter :: nsig=63_i_kind  
!  integer(i_kind),parameter :: ntime=2_i_kind  
  integer(i_kind) :: nxtilde,nytilde
  real(r_kind)    :: rlambda0,pihalf,sign_pole,atilde_x,atilde_y,btilde_x,btilde_y
  real(r_kind)    :: rlon_min_dd,rlon_max_dd,rlat_min_dd,rlat_max_dd
  real(r_kind)    :: rlon_min_ll,rlon_max_ll,rlat_min_ll,rlat_max_ll
  real(r_kind)    :: btilde_xinv,btilde_yinv
  
  integer(i_kind) :: istart,jstart,lon1
   
  !arrays
  integer(i_kind),allocatable :: i0_tilde(:,:),j0_tilde(:,:)
  integer(i_byte),allocatable :: ip_tilde(:,:),jp_tilde(:,:)
  real(r_kind),allocatable    :: xtilde0(:,:),ytilde0(:,:)
  real(r_kind),allocatable    ::cos_beta_ref(:,:),sin_beta_ref(:,:)
  real(r_kind),allocatable    :: glon(:,:),glat(:,:),glon_an(:,:),glat_an(:,:)


  type :: ob
     real(r_kind) :: omf       ! (O-F)
     real(r_kind) :: omf2      ! (O-F)^2
     real(r_kind) :: ob        ! radar observation
     real(r_kind) :: ges       ! Forecast
  end type ob

 CONTAINS


subroutine gridmod_extract
 use constants, only: ione,one,deg2rad
 implicit none
 
 integer(i_kind) :: i,k


  istart=ione
  jstart=ione
  lon1=nlon

! Main driving subroutine which obtains all grid variables necessary for obtaining the
!  x,y locations of the obs on the model grid

! Set up constants
    
       allocate(glon_an(nlon,nlat)) ! unlike rest of gsi these arrays are lon,lat
       allocate(glat_an(nlon,nlat))

       rlon_min_ll=one
       rlat_min_ll=one
       rlon_max_ll=nlon
       rlat_max_ll=nlat
       rlat_min_dd=rlat_min_ll+r1_5
       rlat_max_dd=rlat_max_ll-r1_5
       rlon_min_dd=rlon_min_ll+r1_5
       rlon_max_dd=rlon_max_ll-r1_5
    
       if(1==1) then
          write(6,*)' in grid, rlat_min_dd=',rlat_min_dd
          write(6,*)' in grid, rlat_max_dd=',rlat_max_dd
          write(6,*)' in grid, rlon_min_dd=',rlon_min_dd
          write(6,*)' in grid, rlon_max_dd=',rlon_max_dd
          write(6,*)' in grid, rlat_min_ll=',rlat_min_ll
          write(6,*)' in grid, rlat_max_ll=',rlat_max_ll
          write(6,*)' in grid, rlon_min_ll=',rlon_min_ll
          write(6,*)' in grid, rlon_max_ll=',rlon_max_ll
          write(6,*)' in grid, nlon,nlat=',nlon,nlat
       end if
    
    
!note different allocation for glat_an,glat_an
!need to flip cmaq input glat,glon
  

    
       do k=1,nlon
         do i=1,nlat
            glat_an(k,i)=glat(k,i)*deg2rad
            glon_an(k,i)=glon(k,i)*deg2rad            
         end do
       end do


       call init_general_transform(glat_an,glon_an)

end subroutine gridmod_extract



subroutine tintrp3(f,gout,dxin,dyin,dzin)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    intrp3      linear interpolation in 4 dims
!   prgmmr: parrish          org: np22                date: 1990-10-11
!
! abstract: linear interpolate in 4 dimensions (x,y,z)
!
! program history log:
!   1990-10-11  parrish
!   1998-04-05  weiyu yang
!   1999-08-24  derber, j., treadon, r., yang, w., first frozen mpp version
!   2004-05-18  kleist, documentation
!   2005-02-02  treadon - use ione from constants
!   2008-04-03  safford - rm unused vars         
!   2009-01-23  todling - dim on gridtime is nflds
!
!   input argument list:
!     f        - input interpolator
!     dx,dy,dz - input x,y,z-coords of interpolation points (grid units)
!     n        - number of interpolatees
!
!   output argument list:
!     g        - output interpolatees
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$ 
  use kinds, only: r_kind,i_kind
  use constants, only: zero,one
  implicit none

! Declare passed variables
!  integer(i_kind)                             ,intent(in   ) :: n !number of interpolatees (always 1)
  real(r_kind),dimension(nlon,nlat,nsig),intent(in   ) :: f
  real(r_kind)                          ,intent(in   ) :: dxin,dyin,dzin !dx=x location, dy=y location, dz=closest GRID UNIT height location
  real(r_kind)                          ,intent(  out) :: gout

! Declare local variables

  integer(i_kind) m1,i,ix1,iy1,ix,ixp,iyp,n
  integer(i_kind) iy,iz,izp,itime,itimep,j
  real(r_kind) delx,delyp,delxp
  real(r_kind) dely,delz,delzp
  
  real(r_kind),dimension(1) :: dx,dy,dz,g  
  n=1
  dx(1)=dxin
  dy(1)=dyin
  dz(1)=dzin
  g=zero
  
 
  do i=1,n
     ix1=int(dx(i))
     iy1=int(dy(i))
     
     iz=int(dz(i))
     
     ix1=max(1,min(ix1,nlat)); iz=max(1,min(iz,nsig))  
     delx=dx(i)-float(ix1)
     dely=dy(i)-float(iy1)
     delz=dz(i)-float(iz)
     delx=max(zero,min(delx,one)); delz=max(zero,min(delz,one))
     ix=ix1-istart+2_i_kind
     iy=iy1-jstart+2_i_kind
     
     if(iy<1) then
        iy1=iy1+nlon
        iy=iy1-jstart+2_i_kind
     end if
     if(iy>lon1+1) then
        iy1=iy1-nlon
        iy=iy1-jstart+2_i_kind
     end if
     ixp=ix+1; iyp=iy+1
     
     izp=min(iz+1,nsig)
     if(ix1==nlat) then
        ixp=ix
     end if


     delxp=one-delx; delyp=one-dely
     delzp=one-delz
     g(i) =(f(ix ,iy ,iz )*delxp*delyp*delzp &
          + f(ixp,iy ,iz )*delx*delyp*delzp &
          + f(ix ,iyp,iz )*delxp*dely *delzp &
          + f(ixp,iyp,iz )*delx*dely *delzp &
          + f(ix ,iy ,izp)*delxp*delyp*delz  &
          + f(ixp,iy ,izp)*delx*delyp*delz &
          + f(ix ,iyp,izp)*delxp*dely *delz &
          + f(ixp,iyp,izp)*delx*dely *delz) !No time interpolation done so the second half of this has been removed from the original routine
  end do
  gout=g(1)  !  carley add since we are only interested in a single point


  return
end subroutine tintrp3

subroutine tintrp2a(f,gout,dxin,dyin,nlevs)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    tintrp2a
!   prgmmr: parrish          org: np22                date: 1990-10-11
!
! abstract: linear time interpolate in 3 dimensions (x,y,time) over 
!           n levs
!
! program history log:
!   1990-10-11  parrish
!   1998-04-05  weiyu yang
!   1999-08-24  derber, j., treadon, r., yang, w., first frozen mpp version
!   2003-12-22  kleist, modified to perform 2-d interpolation over a
!                      specified number of vertical levels
!   2004-05-18  kleist, documentation
!   2005-02-02  treadon - use ione from constants
!   2006-04-03  derber  - optimize
!   2008-04-03  safford - rm unused vars
!   2009-01-23  todling - dim on gridtime is nflds
! 
!   input argument list:
!     f        - input interpolator
!     dx,dy    - input x,y,z-coords of interpolation points (grid units)
!     obstime  - time to interpolate to
!     gridtime - grid guess times to interpolate from
!     n        - number of interpolatees
!     nlevs    - number of vertical levels over which to perform the 
!                2-d intrpolation 
!     mype     - mpi task id
!     nflds    - number of guess times available to interpolate from
!
!   output argument list:
!     g        - output interpolatees
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_kind,i_kind
 ! use gridmod, only: istart,jstart,nlon,nlat,lon1,lon2,lat2
  use constants, only: zero,one
  implicit none

! Declare passed variables
  integer(i_kind)                              ,intent(in   ) :: nlevs
  real(r_kind),dimension(nlon,nlat,nsig),intent(in   ) :: f
  real(r_kind),                                 intent(in   ) :: dxin,dyin
  real(r_kind),dimension(nlevs)         ,intent(  out) :: gout

! Declare local variables
  integer(i_kind) i,ix1,iy1,ix,ixp,iyp,n
  integer(i_kind) iy,itime,itimep,j,k
  real(r_kind) delx,delyp,delxp
  real(r_kind) dely
  
  real(r_kind),dimension(1) :: dx,dy
  real(r_kind),dimension(nlevs,1) :: g  
  n=1
  dx(1)=dxin
  dy(1)=dyin
  g=zero  
  
  

  do i=1,n
     ix1=int(dx(i))
     iy1=int(dy(i))
     ix1=max(1,min(ix1,nlat))  
     delx=dx(i)-float(ix1)
     dely=dy(i)-float(iy1)
     delx=max(zero,min(delx,one))
     ix=ix1-istart+2_i_kind
     iy=iy1-jstart+2_i_kind
     if(iy<1) then
        iy1=iy1+nlon
        iy=iy1-jstart+2_i_kind
     end if
     if(iy>lon1+1) then
        iy1=iy1-nlon
        iy=iy1-jstart+2_i_kind
     end if
     ixp=ix+1; iyp=iy+1
     if(ix1==nlat) then
        ixp=ix
     end if
     
     delxp=one-delx; delyp=one-dely
     do k=1,nlevs
        g(k,i)=(f(ix,iy,k)*delxp*delyp+f(ixp,iy,k)*delx*delyp &
              +  f(ix,iyp,k)*delxp*dely+f(ixp,iyp,k)*delx*dely) !No time interpolation done so the second half of this has been removed from the original routine

 
     end do ! end loop over vertical levs
  end do ! end loop over number of locations
  
  gout=g(:,1)  !carley add since we are only interested in 1 interpolatee (n)
  

  return
end subroutine tintrp2a








subroutine tintrp2a_single_level(fin,gout,dxin,dyin)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    tintrp2a
!   prgmmr: parrish          org: np22                date: 1990-10-11
!
! abstract: linear time interpolate in 2 dimensions (x,y) over 
!           1 level with 1 interpolatee !!!CARLEY MODIFIED
!
! program history log:
!   1990-10-11  parrish
!   1998-04-05  weiyu yang
!   1999-08-24  derber, j., treadon, r., yang, w., first frozen mpp version
!   2003-12-22  kleist, modified to perform 2-d interpolation over a
!                      specified number of vertical levels
!   2004-05-18  kleist, documentation
!   2005-02-02  treadon - use ione from constants
!   2006-04-03  derber  - optimize
!   2008-04-03  safford - rm unused vars
!   2009-01-23  todling - dim on gridtime is nflds
! 
!   input argument list:
!     f        - input interpolator
!     dx,dy    - input x,y,z-coords of interpolation points (grid units)
!     obstime  - time to interpolate to
!     gridtime - grid guess times to interpolate from
!     n        - number of interpolatees
!     nlevs    - number of vertical levels over which to perform the 
!                2-d intrpolation 
!     mype     - mpi task id
!     nflds    - number of guess times available to interpolate from
!
!   output argument list:
!     g        - output interpolatees
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_kind,i_kind
 ! use gridmod, only: istart,jstart,nlon,nlat,lon1,lon2,lat2
  use constants, only: zero,one
  implicit none

! Declare passed variables
  real(r_kind),dimension(nlon,nlat),intent(in   ) :: fin
  real(r_kind),                                 intent(in   ) :: dxin,dyin
  real(r_kind),                                 intent(  out) :: gout

! Declare local variables
  integer(i_kind) i,ix1,iy1,ix,ixp,iyp,n,nlevs
  integer(i_kind) iy,itime,itimep,j,k
  real(r_kind) delx,delyp,delxp
  real(r_kind) dely
  
  real(r_kind),dimension(1) :: dx,dy
  real(r_kind),dimension(1,1) :: g
  real(r_kind),dimension(nlon,nlat,1) :: f
  
  f=zero
  f(:,:,1)=fin    
  n=1
  nlevs=1
  dx(1)=dxin
  dy(1)=dyin
  g=zero  
  
  

  do i=1,n
     ix1=int(dx(i))
     iy1=int(dy(i))
     ix1=max(1,min(ix1,nlat))  
     delx=dx(i)-float(ix1)
     dely=dy(i)-float(iy1)
     delx=max(zero,min(delx,one))
     ix=ix1-istart+2_i_kind
     iy=iy1-jstart+2_i_kind
     if(iy<1) then
        iy1=iy1+nlon
        iy=iy1-jstart+2_i_kind
     end if
     if(iy>lon1+1) then
        iy1=iy1-nlon
        iy=iy1-jstart+2_i_kind
     end if
     ixp=ix+1; iyp=iy+1
     if(ix1==nlat) then
        ixp=ix
     end if
     
     delxp=one-delx; delyp=one-dely
     do k=1,nlevs
        g(k,i)=(f(ix,iy,k)*delxp*delyp+f(ixp,iy,k)*delx*delyp &
              +  f(ix,iyp,k)*delxp*dely+f(ixp,iyp,k)*delx*dely)!No time interpolation done so the second half of this has been removed from the original routine
 
     end do ! end loop over vertical levs
  end do ! end loop over number of locations
  
  gout=g(1,1)  !carley add since we are only interested in 1 interpolatee (n)
  

  return
end subroutine tintrp2a_single_level



 subroutine init_general_transform(glats,glons)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    init_general_transform
!   prgmmr:  parrish
!
! abstract:  set up constants to allow conversion between earth lat lon and analysis grid units.
!     There is no need to specify details of the analysis grid projection.  All that is required
!     is the earth latitude and longitude in radians of each analysis grid point.
!
! program history log:
!   2009-08-04  lueken - added subprogram doc block
!   2010-09-08  parrish - replace computation of wind rotation reference angle cos_beta_ref,sin_beta_ref
!                          with new, more accurate and robust version which works for any orientation
!                          of the analysis grid on the sphere (only restriction for now is that
!                          x-y coordinate of analysis grid is right handed).
!
!   input argument list:
!    glons,glats - lons,lats of input grid points of dimesion nlon,nlat
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

  use constants, only: zero,one,half,pi
  implicit none

  real(r_kind)   ,intent(in   ) :: glats(nlon,nlat),glons(nlon,nlat)

  real(r_kind),parameter:: rbig =1.0e30_r_kind
  real(r_kind) xbar_min,xbar_max,ybar_min,ybar_max
  real(r_kind) clon,slon,r_of_lat,xbar,ybar
  integer(i_kind) i,j,istart0,iend,iinc,itemp,ilast,jlast
  real(r_kind),allocatable:: clata(:,:),slata(:,:),clona(:,:),slona(:,:)
  real(r_kind) clat0,slat0,clon0,slon0
  real(r_kind) clat_m1,slat_m1,clon_m1,slon_m1
  real(r_kind) clat_p1,slat_p1,clon_p1,slon_p1
  real(r_kind) x,y,z,xt,yt,zt,xb,yb,zb
  real(r_kind) rlonb_m1,clonb_m1,slonb_m1
  real(r_kind) rlonb_p1,clonb_p1,slonb_p1
  real(r_kind) crot,srot

  pihalf=half*pi

!  define xtilde, ytilde grid, transform

!      glons,glats are lons, lats of input grid points of dimension nlon,nlat
  call get_xytilde_domain(nlon,nlat,glons,glats,nxtilde,nytilde, &
                   xbar_min,xbar_max,ybar_min,ybar_max)
  allocate(i0_tilde(nxtilde,nytilde),j0_tilde(nxtilde,nytilde))
  allocate(ip_tilde(nxtilde,nytilde),jp_tilde(nxtilde,nytilde))
  allocate(xtilde0(nlon,nlat),ytilde0(nlon,nlat))

! define atilde_x, btilde_x, atilde_y, btilde_y

  btilde_x   =(nxtilde -one     )/(xbar_max-xbar_min)
  btilde_xinv=(xbar_max-xbar_min)/(nxtilde -one     )
  atilde_x   =one-btilde_x*xbar_min
  btilde_y   =(nytilde -one     )/(ybar_max-ybar_min)
  btilde_yinv=(ybar_max-ybar_min)/(nytilde -one     )
  atilde_y   =one-btilde_y*ybar_min

! define xtilde0,ytilde0
  do j=1,nlat
     do i=1,nlon
        r_of_lat=pihalf+sign_pole*glats(i,j)
        clon=cos(glons(i,j)+rlambda0)
        slon=sin(glons(i,j)+rlambda0)
        xbar=r_of_lat*clon
        ybar=r_of_lat*slon
        xtilde0(i,j)=atilde_x+btilde_x*xbar
        ytilde0(i,j)=atilde_y+btilde_y*ybar
     end do
  end do

!  now get i0_tilde, j0_tilde, ip_tilde,jp_tilde
  ilast=1 ; jlast=1
  istart0=nxtilde
  iend=1
  iinc=-1
  do j=1,nytilde
     itemp=istart0
     istart0=iend
     iend=itemp
     iinc=-iinc
     ybar=j
     do i=istart0,iend,iinc
        xbar=i
        call nearest_3(ilast,jlast,i0_tilde(i,j),j0_tilde(i,j), &
                       ip_tilde(i,j),jp_tilde(i,j),xbar,ybar,nlon,nlat,xtilde0,ytilde0)
     end do
  end do

!   new, more accurate and robust computation of cos_beta_ref and sin_beta_ref which is independent
!     of sign_pole and works for any orientation of grid on sphere (only restriction for now is that
!     x-y coordinate of analysis grid is right handed).
  allocate(clata(nlon,nlat),slata(nlon,nlat),clona(nlon,nlat),slona(nlon,nlat))
  allocate(cos_beta_ref(nlon,nlat),sin_beta_ref(nlon,nlat))
  do j=1,nlat
     do i=1,nlon
        clata(i,j)=cos(glats(i,j))
        slata(i,j)=sin(glats(i,j))
        clona(i,j)=cos(glons(i,j))
        slona(i,j)=sin(glons(i,j))
     end do
  end do
  do j=1,nlat
     do i=2,nlon-1

!     do all interior lon points to 2nd order accuracy

!   transform so pole is at rlat0,rlon0 and 0 meridian is tangent to earth latitude at rlat0,rlon0.

        clat0=clata(i,j) ; slat0=slata(i,j) ; clon0=clona(i,j) ; slon0=slona(i,j)

!    now obtain new coordinates for m1 and p1 points.

        clat_m1=clata(i-1,j) ; slat_m1=slata(i-1,j) ; clon_m1=clona(i-1,j) ; slon_m1=slona(i-1,j)
        clat_p1=clata(i+1,j) ; slat_p1=slata(i+1,j) ; clon_p1=clona(i+1,j) ; slon_p1=slona(i+1,j)

        x=clat_m1*clon_m1 ; y=clat_m1*slon_m1 ; z=slat_m1
        xt=x*clon0+y*slon0 ; yt=-x*slon0+y*clon0 ; zt=z
        yb=zt*clat0-xt*slat0
        xb=yt
        zb=xt*clat0+zt*slat0

        rlonb_m1=atan2(-yb,-xb)   !  the minus signs here are so line for m1 is directed same
        clonb_m1=cos(rlonb_m1)
        slonb_m1=sin(rlonb_m1)

        x=clat_p1*clon_p1 ; y=clat_p1*slon_p1 ; z=slat_p1
        xt=x*clon0+y*slon0 ; yt=-x*slon0+y*clon0 ; zt=z
        yb=zt*clat0-xt*slat0
        xb=yt
        zb=xt*clat0+zt*slat0
        rlonb_p1=atan2(yb,xb)
        clonb_p1=cos(rlonb_p1)
        slonb_p1=sin(rlonb_p1)
        crot=half*(clonb_m1+clonb_p1)
        srot=half*(slonb_m1+slonb_p1)
        cos_beta_ref(i,j)=crot*clon0-srot*slon0
        sin_beta_ref(i,j)=srot*clon0+crot*slon0
     end do
!               now do i=1 and i=nlon at 1st order accuracy
     i=1

!   transform so pole is at rlat0,rlon0 and 0 meridian is tangent to earth latitude at rlat0,rlon0.

        clat0=clata(i,j) ; slat0=slata(i,j) ; clon0=clona(i,j) ; slon0=slona(i,j)
!    now obtain new coordinates for m1 and p1 points.

        clat_p1=clata(i+1,j) ; slat_p1=slata(i+1,j) ; clon_p1=clona(i+1,j) ; slon_p1=slona(i+1,j)

        x=clat_p1*clon_p1 ; y=clat_p1*slon_p1 ; z=slat_p1
        xt=x*clon0+y*slon0 ; yt=-x*slon0+y*clon0 ; zt=z
        yb=zt*clat0-xt*slat0
        xb=yt
        zb=xt*clat0+zt*slat0
        rlonb_p1=atan2(yb,xb)
        clonb_p1=cos(rlonb_p1)
        slonb_p1=sin(rlonb_p1)
        crot=clonb_p1
        srot=slonb_p1
        cos_beta_ref(i,j)=crot*clon0-srot*slon0
        sin_beta_ref(i,j)=srot*clon0+crot*slon0

     i=nlon

!   transform so pole is at rlat0,rlon0 and 0 meridian is tangent to earth latitude at rlat0,rlon0.

        clat0=clata(i,j) ; slat0=slata(i,j) ; clon0=clona(i,j) ; slon0=slona(i,j)

!    now obtain new coordinates for m1 and p1 points.

        clat_m1=clata(i-1,j) ; slat_m1=slata(i-1,j) ; clon_m1=clona(i-1,j) ; slon_m1=slona(i-1,j)

        x=clat_m1*clon_m1 ; y=clat_m1*slon_m1 ; z=slat_m1
        xt=x*clon0+y*slon0 ; yt=-x*slon0+y*clon0 ; zt=z
        yb=zt*clat0-xt*slat0
        xb=yt
        zb=xt*clat0+zt*slat0

        rlonb_m1=atan2(-yb,-xb)   !  the minus signs here are so line for m1 is directed same
        clonb_m1=cos(rlonb_m1)
        slonb_m1=sin(rlonb_m1)

        crot=clonb_m1
        srot=slonb_m1
        cos_beta_ref(i,j)=crot*clon0-srot*slon0
        sin_beta_ref(i,j)=srot*clon0+crot*slon0
  end do

end subroutine init_general_transform

 subroutine get_xytilde_domain(nx0,ny0,rlons0,rlats0, &
                                  nx,ny,xminout,xmaxout,yminout,ymaxout)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    get_xytilde_domain
!   prgmmr:
!
! abstract:
!
! program history log:
!   2009-08-04  lueken - added subprogram doc block
!
!   input argument list:
!    nx0,ny0
!    rlons0,rlats0
!
!   output argument list:
!    nx,ny
!    xminout,xmaxout,yminout,ymaxout
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

   use constants, only: one,deg2rad,half,zero,r10
!  define parameters for xy domain which optimally overlays input grid

  implicit none
  integer(i_kind),intent(in   ) :: nx0,ny0
  real(r_kind)   ,intent(in   ) :: rlons0(nx0,ny0),rlats0(nx0,ny0)

  integer(i_kind),intent(  out) :: nx,ny
  real(r_kind)   ,intent(  out) :: xminout,xmaxout,yminout,ymaxout

  real(r_kind),parameter:: r37=37.0_r_kind

  real(r_kind) area,areamax,areamin,extra,rlats0max,rlats0min,testlambda
  real(r_kind) xthis,ythis
  integer(i_kind) i,ip1,j,jp1,m

  real(r_kind) coslon0(nx0,ny0),sinlon0(nx0,ny0)
  real(r_kind) coslat0(nx0,ny0),sinlat0(nx0,ny0)
  real(r_kind) count,delbar
  real(r_kind) dx,dy,disti,distj,distmin,distmax
  real(r_kind) xmin,xmax,ymin,ymax

!  get range of lats for input grid

  rlats0max=maxval(rlats0) ; rlats0min=minval(rlats0)

!   assign hemisphere ( parameter sign_pole )

  sign_pole = zero
  if(rlats0min>-r37*deg2rad) sign_pole=-one   !  northern hemisphere xy domain
  if(rlats0max< r37*deg2rad) sign_pole= one   !  southern hemisphere xy domain
  ! if neither condition satisfied (rlat0max > 37N, rlat0min < 37S), try 
  ! this... 
  if (sign_pole == zero) then 
     if (abs(rlats0max) > abs(rlats0min)) then 
        sign_pole=-one  ! NH domain 
     else 
        sign_pole=one   ! SH 
     endif 
  endif


!   get optimum rotation angle rlambda0

  areamin= huge(areamin)
  areamax=-huge(areamax)
  do m=0,359
     testlambda=m*deg2rad
     xmax=-huge(xmax)
     xmin= huge(xmin)
     ymax=-huge(ymax)
     ymin= huge(ymin)
     do j=1,ny0,ny0-1
        do i=1,nx0
           xthis=(pihalf+sign_pole*rlats0(i,j))*cos(rlons0(i,j)+testlambda)
           ythis=(pihalf+sign_pole*rlats0(i,j))*sin(rlons0(i,j)+testlambda)
           xmax=max(xmax,xthis)
           ymax=max(ymax,ythis)
           xmin=min(xmin,xthis)
           ymin=min(ymin,ythis)
        end do
     end do
     do j=1,ny0
        do i=1,nx0,nx0-1
           xthis=(pihalf+sign_pole*rlats0(i,j))*cos(rlons0(i,j)+testlambda)
           ythis=(pihalf+sign_pole*rlats0(i,j))*sin(rlons0(i,j)+testlambda)
           xmax=max(xmax,xthis)
           ymax=max(ymax,ythis)
           xmin=min(xmin,xthis)
           ymin=min(ymin,ythis)
        end do
     end do
     area=(xmax-xmin)*(ymax-ymin)
     areamax=max(area,areamax)
     if(area<areamin) then
        areamin =area
        rlambda0=testlambda
        xmaxout =xmax
        xminout =xmin
        ymaxout =ymax
        yminout =ymin
     end if
  end do


!   now determine resolution of input grid and choose nx,ny of xy grid accordingly
!                 (currently hard-wired at 1/2 the average input grid increment)

  do j=1,ny0
     do i=1,nx0
        coslon0(i,j)=cos(one*rlons0(i,j)) ; sinlon0(i,j)=sin(one*rlons0(i,j))
        coslat0(i,j)=cos(one*rlats0(i,j)) ; sinlat0(i,j)=sin(one*rlats0(i,j))
     end do
  end do

  delbar=zero
  count =zero
  do j=1,ny0-1
     jp1=j+1
     do i=1,nx0-1
        ip1=i+1
        disti=acos(sinlat0(i,j)*sinlat0(ip1,j)+coslat0(i,j)*coslat0(ip1,j)* &
                  (sinlon0(i,j)*sinlon0(ip1,j)+coslon0(i,j)*coslon0(ip1,j)))
        distj=acos(sinlat0(i,j)*sinlat0(i,jp1)+coslat0(i,j)*coslat0(i,jp1)* &
                  (sinlon0(i,j)*sinlon0(i,jp1)+coslon0(i,j)*coslon0(i,jp1)))
        distmax=max(disti,distj)
        distmin=min(disti,distj)
        delbar=delbar+distmax
        count=count+one
     end do
  end do
  delbar=delbar/count
  dx=half*delbar
  dy=dx

!   add extra space to computational grid to push any boundary problems away from
!     area of interest

  extra=r10*dx
  xmaxout=xmaxout+extra
  xminout=xminout-extra
  ymaxout=ymaxout+extra
  yminout=yminout-extra
  nx=1+(xmaxout-xminout)/dx
  ny=1+(ymaxout-yminout)/dy
 
 end subroutine get_xytilde_domain

 subroutine nearest_3(ilast,jlast,i0,j0,ip,jp,x,y,nx0,ny0,x0,y0)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    nearest_3
!   prgmmr:
!
! abstract: find closest 3 points to (x,y) on grid defined by x0,y0
!
! program history log:
!   2009-08-04  lueken - added subprogram doc block
!
!   input argument list:
!    ilast,jlast
!    nx0,ny0
!    x,y
!    x0,y0
!
!   output argument list:
!    ilast,jlast
!    i0,j0
!    ip,jp
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

  implicit none

  integer(i_kind),intent(inout) :: ilast,jlast
  integer(i_kind),intent(  out) :: i0,j0
  integer(i_byte),intent(  out) :: ip,jp
  integer(i_kind),intent(in   ) :: nx0,ny0
  real(r_kind)   ,intent(in   ) :: x,y
  real(r_kind)   ,intent(in   ) :: x0(nx0,ny0),y0(nx0,ny0)
 
  real(r_kind) dista,distb,dist2,dist2min
  integer(i_kind) i,inext,j,jnext

  do
     i0=ilast
     j0=jlast
     dist2min=huge(dist2min)
     inext=0
     jnext=0
     do j=max(j0-1,1),min(j0+1,ny0)
        do i=max(i0-1,1),min(i0+1,nx0)
           dist2=(x-x0(i,j))**2+(y-y0(i,j))**2
           if(dist2<dist2min) then
              dist2min=dist2
              inext=i
              jnext=j
           end if
        end do
     end do
     if(inext==i0.and.jnext==j0) exit
     ilast=inext
     jlast=jnext
  end do

!  now find which way to go in x for second point

  ip=0
  if(i0==nx0)  ip=-1
  if(i0==1) ip=1
  if(ip==0) then
     dista=(x-x0(i0-1,j0))**2+(y-y0(i0-1,j0))**2
     distb=(x-x0(i0+1,j0))**2+(y-y0(i0+1,j0))**2
     if(distb<dista) then
        ip=1
     else
        ip=-1
     end if
  end if

!  repeat for y for 3rd point

  jp=0
  if(j0==ny0  ) jp=-1
  if(j0==1 ) jp=1
  if(jp==0) then
     dista=(x-x0(i0,j0-1))**2+(y-y0(i0,j0-1))**2
     distb=(x-x0(i0,j0+1))**2+(y-y0(i0,j0+1))**2
     if(distb<dista) then
        jp=1
     else
        jp=-1
     end if
  end if

  ilast=i0
  jlast=j0
    
 end subroutine nearest_3
 

!-------------------------------------------------------------------------
!    NOAA/NCEP, National Centers for Environmental Prediction GSI        !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  rotate_wind_ll2xy ---  Rotate earth vector wind
!
! !INTERFACE:
!
  subroutine rotate_wind_ll2xy(u0,v0,u,v,rlon0,x,y)

! !USES:

    use constants, only: one,two,pi,rad2deg,one_tenth
    implicit none

! !INPUT PARAMETERS:

    real(r_kind),intent(in   ) :: u0,v0        ! earth wind component
    real(r_kind),intent(in   ) :: rlon0        ! earth   lon (radians)
    real(r_kind),intent(in   ) :: x,y          ! local x,y coordinate (grid units)

! !OUTPUT PARAMETERS:

    real(r_kind),intent(  out) :: u,v          ! rotated coordinate of winds

! !DESCRIPTION: to convert earth vector wind components to corresponding
!           local x,y coordinate
!
! !REVISION HISTORY:
!   2003-09-30  parrish
!   2004-05-13  kleist, documentation
!   2004-07-15  todling, protex-compliant prologue
!   2010-09-08  parrish, remove sign_pole variable--no longer needed, due to more accurate and
!                 robust computation of reference wind rotation angle defined by
!                 cos_beta_ref, sin_beta_ref.
!
! !REMARKS:
!   language: f90
!   machine:  ibm rs/6000 sp; SGI Origin 2000; Compaq/HP
!
! !AUTHOR:
!   parrish          org: np22                date: 2003-09-30
!
!EOP
!-------------------------------------------------------------------------

  real(r_kind) beta,delx,delxp,dely,delyp
  real(r_kind) sin_beta,cos_beta
  integer(i_kind) ix,iy

!  interpolate departure from longitude part of angle between earth positive east and local positive x

  ix=x
  iy=y
  ix=max(1,min(ix,nlon-1))
  iy=max(1,min(iy,nlat-1))
  delx=x-ix
  dely=y-iy
  delxp=one-delx
  delyp=one-dely
  cos_beta=cos_beta_ref(ix  ,iy  )*delxp*delyp+cos_beta_ref(ix+1,iy  )*delx *delyp+ &
           cos_beta_ref(ix  ,iy+1)*delxp*dely +cos_beta_ref(ix+1,iy+1)*delx *dely
  sin_beta=sin_beta_ref(ix  ,iy  )*delxp*delyp+sin_beta_ref(ix+1,iy  )*delx *delyp+ &
           sin_beta_ref(ix  ,iy+1)*delxp*dely +sin_beta_ref(ix+1,iy+1)*delx *dely
  beta=atan2(sin_beta,cos_beta)

!  now rotate;

  u= u0*cos(beta-rlon0)+v0*sin(beta-rlon0)
  v=-u0*sin(beta-rlon0)+v0*cos(beta-rlon0)

 end subroutine rotate_wind_ll2xy 
 
 
 subroutine destroy_interp
 implicit none  
  deallocate(i0_tilde,j0_tilde)
  deallocate(ip_tilde,jp_tilde)
  deallocate(xtilde0,ytilde0,cos_beta_ref,sin_beta_ref,glon,glat,glon_an,glat_an) 
 end subroutine destroy_interp
 
 

end module interp_util
