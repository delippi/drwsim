program drwsim

  use kinds,     only: r_kind,r_double,i_kind,r_single
  use constants, only: zero,izero,half,one,ione,two,deg2rad,rearth,rad2deg,&
                       one_tenth,r1000,r60,r60inv,r100,r400,init_constants_derived,&
                       init_constants,tiny_r_kind
  use constants, only: flattening,semi_major_axis,grav_ratio,grav,wgtlim,&
                       grav_equator,eccentricity,somigliana
  use interp_util
  use netcdf

  implicit none

  integer(i_kind)     :: nread,ndata,nodata,mindat

  ! Declare local parameters
  real(r_kind),parameter :: four_thirds = 4.0_r_kind / 3.0_r_kind
  real(r_kind),parameter :: r8          = 8.0_r_kind
  real(r_kind),parameter :: r360        = 360.0_r_kind

  !--Derived data type declaration

  type :: radar
     character(4) :: radid
     integer(i_kind) :: vcpnum
     integer(i_kind) :: year
     integer(i_kind) :: month
     integer(i_kind) :: day
     integer(i_kind) :: hour
     integer(i_kind) :: minute
     integer(i_kind) :: second
     real(r_kind) :: radlat
     real(r_kind) :: radlon
     real(r_kind) :: radhgt
     real(r_kind) :: fstgatdis
     real(r_kind) :: gateWidth
     real(r_kind) :: elev_angle
     integer(i_kind) :: num_beam
     integer(i_kind) :: num_gate
     real(r_kind) :: nyq_vel
     real(r_kind) :: azim(360)          !Dims are fixed to facilitate column max calculations
     real(r_kind) :: field(10000,360)   !Dims are fixed to facilitate column max calculations
  end type radar

!--Counters for diagnostics
 integer(i_kind) :: num_missing=izero, &      !counts 
                    numbadtime=izero,num_badtilt=izero, &
                    num_badrange=izero,ibadazm=izero                                                                                                                                                            

!--General declarations
  integer(i_kind) :: ierror,lunrad,i,j,k,v,na,nb,nelv,nvol, &
                     ikx,mins_an,mins_ob,lu
  integer(i_kind) :: maxobs,nchanl,ilat,ilon,maxgate,rad_nelv,cm

  integer(i_kind),dimension(5) :: obdate,iadate
  real(r_kind),dimension(nsig) :: zges,hges
  real(r_kind) :: azm,cosazm,sinazm,cosazm_rw,sinazm_rw,costilt_rw,cosazm_earth,sinazm_earth
  real(r_kind) :: ugesin,vgesin
  real(r_kind) :: zsges,ddiff,observation,rms,bias,sumrms,sumbias
  real(r_kind) :: sin2,termg,termr,termrg,zob,dpres
  real(r_kind) :: a,b,c,ha,epsh,h,aactual,a43,thistilt,x
  real(r_kind) :: thistiltr,selev0,celev0,thisrange,this_stahgt,thishgt
  real(r_kind) :: celev,selev,gamma,thisazimuthr,rlon0, &
                  clat0,slat0,dlat,dlon,thiserr,thislon,thislat, &
                  rlonloc,rlatloc,rlonglob,rlatglob,timeb,rad_per_meter, &
                  delta_gate,delta_az,radar_x,radar_y,radar_lon,radar_lat
  real(r_kind) :: fcstgesin
  real(r_kind) :: radar_twindow                                          


  character(8) cstaid
  character(4) this_staid
  logical   :: outside


  !---------SETTINGS---------!
  integer(i_kind) :: maxobrange=250000_i_kind  ! Range (m) *within* which to use observations 
  integer(i_kind) :: minobrange=20000_i_kind   ! Range (m) *outside* of which
  real(r_kind)    :: mintilt=0.0_r_kind        ! Only use tilt(elevation) angles (deg) >= this number 
  real(r_kind)    :: maxtilt=5.5_r_kind        ! Do no use tilt(elevation) angles (deg) <= this number
  integer(i_kind) :: ithin=4_i_kind            ! Gates to skip for ob thinning (must be >=1)
  real(r_kind)    :: radartwindow=2.5_r_kind   ! Time window within which to grid observations   (minutes)
  character(4)    :: staid='KOUN'
  real(r_kind)    :: mindbz=0_r_kind
  real(r_kind)    :: tilts=0_r_kind
  integer(i_kind) :: azimuths=360_i_kind
  integer(i_kind) :: gatespc=250_i_kind
  integer(i_kind) :: gates=400_i_kind
  !----------------------------------------------!

  !---------NETCDF VARS---------!
  integer(i_kind) :: ncid3d,ncid2d,ncidg,ier
  integer(i_kind) :: pVarId,uVarId,vVarId,wVarId,ghVarId,dbzVarId
  integer(i_kind) :: numsig,numlat,numlon,numtime
  real(r_kind),allocatable  ::    pcoord(    :  ) ! (nsig)
  real(r_kind),allocatable  ::     ges_u(:,:,:,:) ! (nlon,nlat,nsig,ntime)
  real(r_kind),allocatable  ::     ges_v(:,:,:,:) ! (nlon,nlat,nsig,ntime)
  real(r_kind),allocatable  ::     ges_w(:,:,:,:) ! (nlon,nlat,nsig,ntime)
  real(r_kind),allocatable  :: geop_hgtl(:,:,:,:) ! (nlon,nlat,nsig,ntime)
  real(r_kind),allocatable  ::     ges_z(:,:,:,:) ! (nlon,nlat,nsig,ntime)
  real(r_kind),allocatable  ::   ges_dbz(:,:,  :) ! (nlon,nlat,    ,ntime)
  !----------------------------------------------!

  character(len=120)  :: datapath
  character(len=120)  :: nesteddata3d,nesteddata2d,nestedgrid
  character(len=120)  :: bufroutfile

  namelist/drw/iadate,radartwindow,nesteddata3d,nesteddata2d,nestedgrid,bufroutfile,&
               mintilt,maxtilt,staid,mindbz,tilts,maxobrange,minobrange,azimuths,ithin,&
               gatespc,datapath


!--------------------------------------------------------------------------------------!
!                            END OF ALL DECLARATIONS
!                            !
!--------------------------------------------------------------------------------------!

  !--Set up the constants module used here
  call init_constants_derived
  call init_constants(.true.)    !initialize regional constants
  mindat=ione

  !----READ NAMELIST----
  read(*,drw)
  !----READ NAMELIST----

  gatespc=gatespc*ithin
  gates=int(maxobrange/gatespc) !calculate num gates based on namelist settings.
  nesteddata3d=trim(datapath) // trim(nesteddata3d) !concat strings
  nesteddata2d=trim(datapath) // trim(nesteddata2d)
  nestedgrid  =trim(datapath) // trim(nestedgrid)

  print *, "----NC FILE SPECS--------"
  print *, "nlat",nlat
  print *, "nlon",nlon
  print *, "nsig",nsig
  print *, "grid spacing = 3km"

  print *, "----Namelist settings----"
  print *, "staid       ",staid
  print *, "tilts       ",tilts
  print *, "azimuths    ",azimuths
  print *, "gates       ",gates
  print *, "ithin       ",ithin
  print *, "gatespc     ",gatespc
  print *, "maxobrange  ",maxobrange
  print *, "minobrange  ",minobrange
  print *, "mintilt     ",mintilt
  print *, "maxtilt     ",maxtilt
  print *, "mindbz      ",mindbz
  print *, "radartwindow",radartwindow
  print *, "datapath:",datapath
  print *, "nesteddata3d:",nesteddata3d
  print *, "nesteddata2d:",nesteddata2d
  print *, "nestedgrid:",nestedgrid

  print *, "deg2rad",deg2rad
  print *, "rearth",rearth
  if(deg2rad>0) print *, "(derived constants succesfully initialized)"
  if(rearth==6370.e03_r_kind) print *, "(regional consts sucess init)"

  if (minobrange >= maxobrange) then
     STOP 'MINIMUM OB RANGE >= MAXIMUM OB RANGE. PROGRAM STOPPING NOW drwsim.f90'
  end if
  if (ithin < 1) then
     STOP 'ithin MUST BE >=1! CHECK NAMELIST AND RESET.'
  end if


!------ OPEN 3D NETCDF FILE ------------------------------
  ncid3d=30
  ! Open the 3D netcdf file
  ier=nf90_open(nesteddata3d,nf90_NoWrite,ncid3d)
  if(ier /= nf90_NoErr) STOP "ERROR READ NETCDF STOP!"

  ! u
  allocate(ges_u(nlon,nlat,nsig,ntime))
  print *, "Reading ges_u"
  ges_u=0.0_r_kind !initialize
  ier=nf90_inq_varid(ncid3d,"ucomp",uVarId)
  if(ier /= nf90_NoErr) STOP "ERROR INQ ges_u NETCDF STOP!"
  ier=nf90_get_var(ncid3d,uVarId,ges_u)
  if(ier /= nf90_NoErr) STOP "ERROR GET ges_u NETCDF STOP!"
  ! check the order of the array (nlon,nlat,nsig,ntime)
  !ier=nf90_inquire_dimension(ncid3d, 1, len = numlon)
  !ier=nf90_inquire_dimension(ncid3d, 2, len = numlat)
  !ier=nf90_inquire_dimension(ncid3d, 3, len = numsig)
  !ier=nf90_inquire_dimension(ncid3d, 5, len = numtime)
  !print *, numtime,numsig,numlat,numlon
  !print *, ntime,nsig,nlat,nlon



   print *, maxval(ges_u(:,:,:,:)) 
   

   print *, "end of program" 
end program drwsim
