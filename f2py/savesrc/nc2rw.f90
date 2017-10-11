!-------------------------------------------------------------------------
!    NOAA/NCEP, National Centers for Environmental Prediction GSI        !
!-------------------------------------------------------------------------
!
!
!   COPIED FROM GRIDMOD.F90 ON 01/29/2012 FROM BRACNH EXP-ENSVAR-UPDATES
!        THIS VERSION HAS SOME MODIFICATIONS TO FIX INTERPOLATION 
!        ISSUES.
!
!
!BOP
!
! !IROUTINE:  tll2xy --- convert earth lon-lat to x-y grid coordinates
!
! !INTERFACE:
!
  subroutine tll2xy(rlon,rlat,x,y,outside)

! !USES:

    use constants, only: one
    use interp_util
    implicit none

    real*8,intent(in   ) :: rlon  ! earth longitude (radians)
    real*8,intent(in   ) :: rlat  ! earth latitude  (radians)

! !OUTPUT PARAMETERS:

    real*8,intent(  out) :: x  ! x-grid coordinate (grid units)
    real*8,intent(  out) :: y  ! y-grid coordinate (grid units)
    logical     ,intent(  out) :: outside     ! .false., then point is inside x-y domain
                                              ! .true.,  then point is outside x-y domain

! !DESCRIPTION: to convert earth lon-lat to x-y grid units of a 
!           general regional rectangular domain.  Also, decide if
!           point is inside this domain.  As a result, there is
!           no restriction on type of horizontal coordinate for
!           a regional run, other than that it not have periodicity
!           or polar singularities.
!           This is done by first converting rlon, rlat to an
!           intermediate coordinate xtilde,ytilde, which has
!           precomputed pointers and constants for final conversion
!           to the desired x,y via 3 point inverse interpolation.
!           All of the information needed is derived from arrays
!           specifying earth latitude and longitude of every point
!           on the input grid.  Currently, the input x-y grid that
!           this is based on must be non-staggered.  This restriction
!           will eventually be lifted so we can run directly from
!           model grids that are staggered without first resorting
!           to interpolation of the guess to a non-staggered grid.
!
! !REVISION HISTORY:
!   2003-08-28  parrish
!   2004-05-13  kleist, documentation
!   2004-07-15  todling, protex-compliant prologue
!   2004-07-23  parrish - new routine

!
! !REMARKS:
!   language: f90
!   machine:  ibm rs/6000 sp; SGI Origin 2000; Compaq/HP
!
! !AUTHOR:
!   parrish          org: np22                date: 2003-08-28
!
!EOP
!-------------------------------------------------------------------------

    real*8 clon,slon,r_of_lat,xtilde,ytilde
    real*8 dtilde,etilde
    real*8 d1tilde,d2tilde,e1tilde,e2tilde,detinv
    integer*4 itilde,jtilde
    integer*4 i0,j0,ip,jp

!   first compute xtilde, ytilde

    clon=cos(rlon+rlambda0)
    slon=sin(rlon+rlambda0)
    r_of_lat=pihalf+sign_pole*rlat

    xtilde=atilde_x+btilde_x*r_of_lat*clon
    ytilde=atilde_y+btilde_y*r_of_lat*slon

!  next get interpolation information

    itilde=max(1,min(nint(xtilde),nxtilde))
    jtilde=max(1,min(nint(ytilde),nytilde))

    i0     =   i0_tilde(itilde,jtilde)
    j0     =   j0_tilde(itilde,jtilde)
    ip     =i0+ip_tilde(itilde,jtilde)
    jp     =j0+jp_tilde(itilde,jtilde)
    dtilde =xtilde-xtilde0(i0,j0)
    etilde =ytilde-ytilde0(i0,j0)
    d1tilde=(xtilde0(ip,j0)-xtilde0(i0,j0))*(ip-i0)
    d2tilde=(xtilde0(i0,jp)-xtilde0(i0,j0))*(jp-j0)
    e1tilde=(ytilde0(ip,j0)-ytilde0(i0,j0))*(ip-i0)
    e2tilde=(ytilde0(i0,jp)-ytilde0(i0,j0))*(jp-j0)
    detinv =one/(d1tilde*e2tilde-d2tilde*e1tilde)
    x = i0+detinv*(e2tilde*dtilde-d2tilde*etilde)
    y = j0+detinv*(d1tilde*etilde-e1tilde*dtilde)
    if (i0 == ip .and. j0 == jp) then ! ob at center of domain. 
       x = i0; y = j0 
    else 
       x = i0+detinv*(e2tilde*dtilde-d2tilde*etilde) 
       y = j0+detinv*(d1tilde*etilde-e1tilde*dtilde) 
    endif

    outside=x < rlon_min_dd .or. x > rlon_max_dd .or. &
            y < rlat_min_dd .or. y > rlat_max_dd

 end subroutine tll2xy



SUBROUTINE invtllv(ALM,APH,TLMO,CTPH0,STPH0,TLM,TPH)
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    invtllv
!
!   prgrmmr:
!
! abstract:  inverse of tllv:  input ALM,APH is rotated lon,lat
!                   output is earth lon,lat, TLM,TPH
!
! program history log:
!   2008-03-25  safford -- add subprogram doc block, rm unused uses
!
!   input argument list:
!     alm   -- input earth longitude
!     aph   -- input earth latitude
!     tlmo  -- input earth longitude of rotated grid origin (radrees)
!     ctph0 -- cos(earth lat of rotated grid origin)
!     stph0 -- sin(earth lat of rotated grid origin)
!
!   output argument list:
!     tlm   -- rotated grid longitude
!     tph   -- rotated grid latitude
!
! attributes:
!   language:  f90
!   machine:
!
!$$$ end documentation block

  use interp_util
  implicit none

  real*8,intent(in   ) :: alm,aph,tlmo,ctph0,stph0
  real*8,intent(  out) :: tlm,tph

  real*8:: relm,srlm,crlm,sph,cph,cc,anum,denom

  RELM=ALM
  SRLM=SIN(RELM)
  CRLM=COS(RELM)
  SPH=SIN(APH)
  CPH=COS(APH)
  CC=CPH*CRLM
  ANUM=CPH*SRLM
  DENOM=CTPH0*CC-STPH0*SPH
  TLM=tlmo+ATAN2(ANUM,DENOM)
  TPH=ASIN(CTPH0*SPH+STPH0*CC)
  
END SUBROUTINE invtllv
