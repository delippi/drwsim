module constants
!$$$   module documentation block
!                .      .    .                                       .
! module:    constants
!   prgmmr: treadon          org: np23                date: 2003-09-25
!
! abstract:  This module contains the definition of various constants
!            used in the gsi code
!
! program history log:
!   2003-09-25  treadon - original code
!   2004-03-02  treadon - allow global and regional constants to differ
!   2004-06-16  treadon - update documentation
!   2004-10-28  treadon - replace parameter tiny=1.e-12 with tiny
!                         and tiny_single
!   2004-11-16  treadon - add huge_single, huge parameters
!   2005-01-27  cucurull - add ione
!   2005-08-24  derber   - move cg_term to constants from qcmod
!   2006-03-07  treadon  - add rd_over_cp_mass
!   2006-05-18  treadon  - add huge
!   2006-06-06       su  - add var-qc wgtlim, change value to 0.25 (ECMWF)
!   2006-07-28  derber   - add r1000
!   2007-03-20  rancic   - add r3600
!   2009-02-05  cucurull - modify refractive indexes for gpsro data
!   2010-08-25  cucurull - add constants to compute compressibility factor
!                        - add option to use Rueger/Bevis refractive index coeffs
!   2010-12-20 pagowski  - add max_varname_length=12
!   2010-04-01 li        - add maximum diurnal thermocline thickness
!
! Subroutines Included:
!   sub init_constants_derived - compute derived constants
!   sub init_constants         - set regional/global constants
!   sub gps_constants          - set Rueger/Bevis refractive index coefficients
!
! Variable Definitions:
!   see below
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$ end documentation block

  implicit none

! set default as private
!  private
! set subroutines as public
  public :: init_constants_derived
  public :: init_constants
  public :: gps_constants
! set passed variables to public
  public :: one,two,ione,half,zero,izero,deg2rad,pi,three,quarter,one_tenth
  public :: rad2deg,zero_quad,r3600,r1000,r60inv,five,four,rd_over_cp,grav
  public :: rd,rozcon,rearth_equator,zero_single,tiny_r_kind,tiny_single,ten
  public :: omega,rcp,rearth,fv,h300,cp,cg_term,tpwcon,xb,ttp,psatk,xa,tmix
  public :: xai,xbi,psat,eps,omeps,wgtlim,one_quad,two_quad,epsq,climit,epsm1,hvap
  public :: hsub,cclimit,el2orc,elocp,h1000,cpr,pcpeff0,pcpeff2,delta,pcpeff1
  public :: factor1,c0,pcpeff3,factor2,dx_inv,dx_min,rhcbot,rhctop,hfus,ke2
  public :: rrow,cmr,cws,r60,huge_i_kind,huge_r_kind,t0c,rd_over_cp_mass
  public :: somigliana,grav_equator,grav_ratio,flattening,semi_major_axis
  public :: n_b,n_a,eccentricity,huge_single,constoz,g_over_rd,amsua_clw_d2
  public :: amsua_clw_d1,n_c,rd_over_g,zero_ilong
  public :: r10,r100,sqrt_tiny_r_kind,r2000,r4000
  public :: r0_01,r0_02,r0_03,r0_04,r0_05,r400,r2400
  public :: cpf_a0, cpf_a1, cpf_a2, cpf_b0, cpf_b1, cpf_c0, cpf_c1, cpf_d, cpf_e
  public :: psv_a, psv_b, psv_c, psv_d
  public :: ef_alpha, ef_beta, ef_gamma
  public :: max_varname_length
  public :: z_w_max

! Declare derived constants
  integer*4:: huge_i_kind
  integer*4, parameter :: max_varname_length=12
  real*4:: tiny_single, huge_single
  real*8:: xai, xa, xbi, xb, dldt, rozcon,ozcon,fv, tpwcon,eps, rd_over_g
  real*8:: el2orc, g_over_rd, rd_over_cp, cpr, omeps, epsm1, factor2
  real*8:: factor1, huge_r_kind, tiny_r_kind, deg2rad, pi, rad2deg, cg_term
  real*8:: eccentricity_linear, cv, rv, rd_over_cp_mass, cliq, rd, cp_mass
  real*8:: eccentricity, grav, rearth, r60inv
  real*8:: sqrt_tiny_r_kind
  real*8:: n_a, n_b, n_c

! Define constants common to global and regional applications
  real*8,parameter::  rearth_equator= 6.37813662e6  ! equatorial earth radius          (m)
  real*8,parameter::  omega  = 7.2921e-5            !  angular velocity of earth       (1/s)
  real*8,parameter::  cp     = 1.0046e+3            !  specific heat of air @pressure  (J/kg/K)
  real*8,parameter::  cvap   = 1.8460e+3            !  specific heat of h2o vapor      (J/kg/K)
  real*8,parameter::  csol   = 2.1060e+3            !  specific heat of solid h2o (ice)(J/kg/K)
  real*8,parameter::  hvap   = 2.5000e+6            !  latent heat of h2o condensation (J/kg)
  real*8,parameter::  hfus   = 3.3358e+5            !  latent heat of h2o fusion       (J/kg)
  real*8,parameter::  psat   = 6.1078e+2            !  pressure at h2o triple point    (Pa)
  real*8,parameter::  t0c    = 2.7315e+2            !  temperature at zero celsius     (K)
  real*8,parameter::  ttp    = 2.7316e+2            !  temperature at h2o triple point (K)
  real*8,parameter::  jcal   = 4.1855e+0            !  joules per calorie              ()
  real*8,parameter::  stndrd_atmos_ps = 1013.25e2   ! 1976 US standard atmosphere ps   (Pa)

! Numeric constants
  integer*4,parameter::  izero  = 0
  integer*4,parameter::  ione   = 1

  integer*4,parameter::  zero_ilong = 0

  real*4,parameter::  zero_single= 0.0

  real*8,parameter::  zero      = 0.0
  real*8,parameter::  r0_01     = 0.01
  real*8,parameter::  r0_02     = 0.02
  real*8,parameter::  r0_03     = 0.03
  real*8,parameter::  r0_04     = 0.04
  real*8,parameter::  r0_05     = 0.05
  real*8,parameter::  one_tenth = 0.10
  real*8,parameter::  quarter   = 0.25
  real*8,parameter::  one       = 1.0
  real*8,parameter::  two       = 2.0
  real*8,parameter::  three     = 3.0
  real*8,parameter::  four      = 4.0
  real*8,parameter::  five      = 5.0
  real*8,parameter::  ten       = 10.0
  real*8,parameter::  r10       = 10.0
  real*8,parameter::  r60       = 60.
  real*8,parameter::  r100      = 100.0
  real*8,parameter::  r400      = 400.0
  real*8,parameter::  r1000     = 1000.0
  real*8,parameter::  r2000     = 2000.0
  real*8,parameter::  r2400     = 2400.0
  real*8,parameter::  r4000     = 4000.0
  real*8,parameter::  r3600     = 3600.0

! maximum diurnal thermocline thickness
  real*8,parameter:: z_w_max   = 30.0

  real*16,parameter::  zero_quad = 0.0
  real*16,parameter::  one_quad  = 1.0
  real*16,parameter::  two_quad  = 2.0

! Constants for compressibility factor (Davis et al 1992)
  real*8,parameter::  cpf_a0 =  1.58123e-6 ! K/Pa
  real*8,parameter::  cpf_a1 = -2.9331e-8  ! 1/Pa
  real*8,parameter::  cpf_a2 =  1.1043e-10 ! 1/K 1/Pa
  real*8,parameter::  cpf_b0 =  5.707e-6   ! K/Pa
  real*8,parameter::  cpf_b1 = -2.051e-8   ! 1/Pa
  real*8,parameter::  cpf_c0 =  1.9898e-4  ! K/Pa
  real*8,parameter::  cpf_c1 = -2.376e-6   ! 1/Pa
  real*8,parameter::  cpf_d  =  1.83e-11   ! K2/Pa2
  real*8,parameter::  cpf_e  = -0.765e-8   ! K2/Pa2

! Constants for vapor pressure at saturation
  real*8,parameter::  psv_a =  1.2378847e-5       !  (1/K2)
  real*8,parameter::  psv_b = -1.9121316e-2       !  (1/K)
  real*8,parameter::  psv_c = 33.93711047         !
  real*8,parameter::  psv_d = -6.3431645e+3       !  (K)

! Constants for enhancement factor to calculating the mole fraction of water vapor
  real*8,parameter::  ef_alpha = 1.00062           !
  real*8,parameter::  ef_beta  = 3.14e-8           !  (1/Pa)
  real*8,parameter::  ef_gamma = 5.6e-7            !  (1/K2)

! Parameters below from WGS-84 model software inside GPS receivers.
  real*8,parameter::  semi_major_axis = 6378.1370e3     !                     (m)
  real*8,parameter::  semi_minor_axis = 6356.7523142e3  !                     (m)
  real*8,parameter::  grav_polar      = 9.8321849378    !                     (m/s2)
  real*8,parameter::  grav_equator    = 9.7803253359    !                     (m/s2) 
  real*8,parameter::  earth_omega     = 7.292115e-5     !                     (rad/s)
  real*8,parameter::  grav_constant   = 3.986004418e14  !                     (m3/s2)

! Derived geophysical constants
  real*8,parameter::  flattening = (semi_major_axis-semi_minor_axis)/semi_major_axis
  real*8,parameter::  somigliana = &
       (semi_minor_axis/semi_major_axis) * (grav_polar/grav_equator) - one
  real*8,parameter::  grav_ratio = (earth_omega*earth_omega * &
       semi_major_axis*semi_major_axis * semi_minor_axis) / grav_constant 

! Derived thermodynamic constants
  real*8,parameter::  dldti = cvap-csol
  real*8,parameter::  hsub = hvap+hfus
  real*8,parameter::  psatk = psat*0.001
  real*8,parameter::  tmix = ttp-20.
  real*8,parameter::  elocp = hvap/cp
  real*8,parameter::  rcp  = one/cp

! Constants used in GFS moist physics
  real*8,parameter::  h300 = 300.
  real*8,parameter::  half = 0.5
  real*8,parameter::  cclimit = 0.001
  real*8,parameter::  climit = 1.e-20
  real*8,parameter::  epsq = 2.e-12
  real*8,parameter::  h1000 = r1000
  real*8,parameter::  rhcbot=0.85
  real*8,parameter::  rhctop=0.85
  real*8,parameter::  dx_max=-8.8818363
  real*8,parameter::  dx_min=-5.2574954
  real*8,parameter::  dx_inv=one/(dx_max-dx_min)
  real*8,parameter::  c0=0.002
  real*8,parameter::  delta=0.6077338
  real*8,parameter::  pcpeff0=1.591
  real*8,parameter::  pcpeff1=-0.639
  real*8,parameter::  pcpeff2=0.0953
  real*8,parameter::  pcpeff3=-0.00496
  real*8,parameter::  cmr = one/0.0003
  real*8,parameter::  cws = 0.025
  real*8,parameter::  ke2 = 0.00002
  real*8,parameter::  row = r1000
  real*8,parameter::  rrow = one/row

! Constant used to process ozone
  real*8,parameter::  constoz = 604229.0

! Constants used in cloud liquid water correction for AMSU-A
! brightness temperatures
  real*8,parameter::  amsua_clw_d1 = 0.754
  real*8,parameter::  amsua_clw_d2 = -2.265

! Constants used for variational qc
  real*8,parameter::  wgtlim = quarter  ! Cutoff weight for concluding that obs has been
                                     ! rejected by nonlinear qc. This limit is arbitrary
                                     ! and DOES NOT affect nonlinear qc. It only affects
                                     ! the printout which "counts" the number of obs that
                                     ! "fail" nonlinear qc.  Observations counted as failing
                                     ! nonlinear qc are still assimilated.  Their weight
                                     ! relative to other observations is reduced. Changing
                                     ! wgtlim does not alter the analysis, only
                                     ! the nonlinear qc data "count"

contains

  subroutine init_constants_derived
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    init_constants_derived          set derived constants
!     prgmmr:    treadon          org: np23           date: 2004-12-02
!
! abstract:  This routine sets derived constants
!
! program history log:
!   2004-12-02  treadon
!   2005-03-03  treadon - add implicit none
!   2008-06-04  safford - rm unused vars
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp
!
!$$$
    implicit none

!   Trigonometric constants
    pi      = acos(-one)
    deg2rad = pi/180.0
    rad2deg = one/deg2rad
    cg_term = (sqrt(two*pi))/two                  ! constant for variational qc
    tiny_r_kind = tiny(zero)
    sqrt_tiny_r_kind = r10*sqrt(tiny_r_kind)
    huge_r_kind = huge(zero)
    tiny_single = tiny(zero_single)
    huge_single = huge(zero_single)
    huge_i_kind = huge(izero)
    r60inv=one/r60

!   Geophysical parameters used in conversion of geopotential to
!   geometric height
    eccentricity_linear = sqrt(semi_major_axis**2 - semi_minor_axis**2)
    eccentricity = eccentricity_linear / semi_major_axis

    return
  end subroutine init_constants_derived

  subroutine init_constants(regional)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    init_constants       set regional or global constants
!     prgmmr:    treadon          org: np23           date: 2004-03-02
!
! abstract:  This routine sets constants specific to regional or global
!            applications of the gsi
!
! program history log:
!   2004-03-02  treadon
!   2004-06-16  treadon, documentation
!   2004-10-28  treadon - use intrinsic TINY function to set value
!                         for smallest machine representable positive
!                         number
!   2004-12-03  treadon - move derived constants to init_constants_derived
!   2005-03-03  treadon - add implicit none
!
!   input argument list:
!     regional - if .true., set regional gsi constants;
!                otherwise (.false.), use global constants
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp
!
!$$$
    implicit none

    logical,intent(in   ) :: regional

    real*8 reradius,g,r_d,r_v,cliq_wrf

!   Define regional constants here
    if (regional) then

!      Name given to WRF constants
       reradius = one/6370.e03
       g        = 9.81
       r_d      = 287.04
       r_v      = 461.6
       cliq_wrf = 4190.0
       cp_mass  = 1004.67

!      Transfer WRF constants into unified GSI constants
       rearth = one/reradius
       grav   = g
       rd     = r_d
       rv     = r_v
       cv     = cp-r_d
       cliq   = cliq_wrf
       rd_over_cp_mass = rd / cp_mass

!   Define global constants here
    else
       rearth = 6.3712e+6
       grav   = 9.80665e+0
       rd     = 2.8705e+2
       rv     = 4.6150e+2
       cv     = 7.1760e+2
       cliq   = 4.1855e+3
       cp_mass= zero
       rd_over_cp_mass = zero
    endif


!   Now define derived constants which depend on constants
!   which differ between global and regional applications.

!   Constants related to ozone assimilation
    ozcon = grav*21.4e-9
    rozcon= one/ozcon

!   Constant used in vertical integral for precipitable water
    tpwcon = 100.0/grav

!   Derived atmospheric constants
    fv         = rv/rd-one    ! used in virtual temperature equation 
    dldt       = cvap-cliq
    xa         = -(dldt/rv)
    xai        = -(dldti/rv)
    xb         = xa+hvap/(rv*ttp)
    xbi        = xai+hsub/(rv*ttp)
    eps        = rd/rv
    epsm1      = rd/rv-one
    omeps      = one-eps
    factor1    = (cvap-cliq)/rv
    factor2    = hvap/rv-factor1*t0c
    cpr        = cp*rd
    el2orc     = hvap*hvap/(rv*cp)
    rd_over_g  = rd/grav
    rd_over_cp = rd/cp
    g_over_rd  = grav/rd

    return
  end subroutine init_constants

  subroutine gps_constants(use_compress)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    gps_constants     set Bevis or Rueger refractive index coeff
!     prgmmr:    cucurull          org: np23           date: 2010-08-25
!
! abstract:  This routine sets constants for the refractivity equation. GSI uses Bevis 
!            coefficients when the compressibility factors option is turned off 
!            and uses Rueger coefficients otherwise.
!
! program history log:
!   2010-08-25  cucurull
!   2010-08-25  cucurull, documentation
!
!   input argument list:
!     compress - if .true., set Rueger coefficients;
!                otherwise (.false.), use Bevis coefficients
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp
!
!$$$
    implicit none

    logical,intent(in   ) :: use_compress

!   Define refractive index coefficients here
    if (use_compress) then

       ! Constants for gpsro data (Rueger 2002)
       n_a = 77.6890   ! K/mb
       n_b = 3.75463e+5  ! K^2/mb
       n_c = 71.2952   ! K/mb
    else
       ! Constants for gpsro data (Bevis et al 1994)
       n_a = 77.60     ! K/mb
       n_b = 3.739e+5  ! K^2/mb
       n_c = 70.4      ! K/mb
    endif

    return
  end subroutine gps_constants

end module constants
