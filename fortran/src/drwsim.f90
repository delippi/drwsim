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

  integer(i_kind)     :: ndata,mindat

  ! Declare local parameters
  real(r_kind),parameter :: four_thirds = 4.0_r_kind / 3.0_r_kind
  real(r_kind),parameter :: r8          = 8.0_r_kind
  real(r_kind),parameter :: r360        = 360.0_r_kind

!--General declarations
  integer(i_kind) :: n,iret
  integer(i_kind),dimension(4) :: iadate,intdate
  character(4) :: yyyy
  character(2) :: mm,dd,hh,mn
  character(10)             :: idate
  real(r_kind),dimension(nsig) :: zges,hges
  real(r_kind) :: cosazm,sinazm,costilt,sintilt
  real(r_kind) :: ugesin,vgesin,wgesin,dbzgesin
  real(r_kind) :: zsges,psfcsges
  real(r_kind) :: sin2,termg,termr,termrg,dpres
  real(r_kind) :: b,c,ha,epsh,h,aactual,a43,thistilt
  real(r_kind) :: thistiltr,selev0,celev0,thisrange,this_stahgt,thishgt
  real(r_kind) :: celev,selev,gamma,thisazimuthr,thisazimuth,rlon0,rlat0,stahgt, &
                  clat0,slat0,dlat,dlon,thislon,thislat, &
                  rlonloc,rlatloc,rlonglob,rlatglob,rad_per_meter, &
                  radar_x,radar_y,radar_lon,radar_lat
  real(r_kind),allocatable :: delz(:,:),height(:,:)
  real(r_kind),allocatable :: drwpol(:,:,:) !tilt,azm,gate
  integer(i_kind) :: irid,itilt,iazm,igate,itime,iazm90,isig


  character(4) this_staid
  character(5) str_gatespc
  logical   :: diagprint,inside,bufrisopen,radar_location
  integer   :: diagverbose

!  type(radar),allocatable :: strct_in_vel(:,:),rad(:)


  !---------DEFAULT SETTINGS---------!
  integer(i_kind) :: maxobrange=250000_i_kind  ! Range (m) *within* which to use observations 
  integer(i_kind) :: minobrange=20000_i_kind   ! Range (m) *outside* of which
  real(r_kind)    :: mintilt=0.0_r_kind        ! Only use tilt(elevation) angles (deg) >= this number 
  real(r_kind)    :: maxtilt=5.5_r_kind        ! Do no use tilt(elevation) angles (deg) <= this number
  integer(i_kind) :: ithin=4_i_kind            ! Gates to skip for ob thinning (must be >=1)
  character(4)    :: staid='KOUN'              ! default station ID to use
  real(r_kind)    :: mindbz=0_r_kind           ! minimum dbz value needed at a location to create a drw
  real(r_kind),dimension(25)    :: tilts=0_r_kind ! initialize the tilts to zero-degrees
  integer(i_kind) :: azimuths=360_i_kind       ! number of azimuths
  integer(i_kind) :: gatespc=250_i_kind        ! gate spacing (meters)
  integer(i_kind) :: numgates=400_i_kind       ! number of gates
  integer(i_kind) :: ntime=1_i_kind            ! number of times from nc file
  integer(i_kind) :: nelv=1_i_kind             ! number of elvation tilts
  !----------------------------------------------!

  !---------NETCDF VARS---------!
  integer(i_kind) :: ncid3d,ncid2d,ncidgs,ier
  integer(i_kind) :: tVarId,pVarId,uVarId,vVarId,wVarId,ghVarId,dbzVarId !3d
  integer(i_kind) :: hgtVarID,psfcVarID                                  !2d
  integer(i_kind) :: glonVarId,glatVarId                                 !gs
  real(r_kind),allocatable  ::      time(      :) ! (    ,    ,    ,ntime)
  real(r_kind),allocatable  ::    pcoord(    :  ) ! (    ,    ,nsig,     )
  real(r_kind),allocatable  ::     ges_u(:,:,:,:) ! (nlon,nlat,nsig,ntime)
  real(r_kind),allocatable  ::     ges_v(:,:,:,:) ! (nlon,nlat,nsig,ntime)
  real(r_kind),allocatable  ::     ges_w(:,:,:,:) ! (nlon,nlat,nsig,ntime)
  real(r_kind),allocatable  :: geop_hgtl(:,:,:,:) ! (nlon,nlat,nsig,ntime)
  real(r_kind),allocatable  ::   ges_dbz(:,:,:,:) ! (nlon,nlat,nsig,ntime)
  real(r_kind),allocatable  ::     ges_z(:,:    ) ! (nlon,nlat,    ,     )
  real(r_kind),allocatable  ::  ges_psfc(:,:    ) ! (nlon,nlat,    ,     )
!  real(r_kind),allocatable  ::      lons(:,:    ) ! (nlon,nlat,    ,     )
!  real(r_kind),allocatable  ::      lats(:,:    ) ! (nlon,nlat,    ,     )
  !----------------------------------------------!

  !---------GLOBAL RADAR CSV FILE VARS---------!
  integer(i_kind)                         :: ii,numradars
  real(r_kind),dimension(:),allocatable   :: dflat,dflon,dfheight
  character(12),dimension(:),allocatable  :: dfid
  !----------------------------------------------!


  !---------BUFR VARS--------------------------!
  integer(i_kind) :: itiltbufr,iazmbufr,igatebufr,iazmbufr90
  character(80) :: bufrfilename,hdstr,obstr
  real(r_double) :: hdr(12)
  real(r_double),allocatable :: obs(:,:)
  character(8) :: chdr,subset
  equivalence (hdr(1),chdr)

  !---------L2RWBUFR CSV TABLE VARS------------!
  character(10) :: message_type
  character(8) :: cdummy

  !-------------TIMER VARS---------------------!
  integer(i_kind) :: time_array_0(8), time_array_1(8),hrs,mins,secs
  real(r_kind)    :: start_time,finish_time,total_time


  character(len=120)  :: datapath
  character(len=120)  :: nesteddata3d,nesteddata2d,nestedgrid,radarcsv
  character(len=120)  :: bufroutfile

  namelist/drw/iadate,ntime,nelv,nesteddata3d,nesteddata2d,nestedgrid,bufroutfile,&
               mintilt,maxtilt,staid,mindbz,tilts,maxobrange,minobrange,azimuths,ithin,&
               gatespc,datapath,diagprint,diagverbose,radarcsv


!--------------------------------------------------------------------------------------!
!                            END OF ALL DECLARATIONS
!                            !
!--------------------------------------------------------------------------------------!

  !--Call Timer
  call date_and_time(values=time_array_0)
  start_time = time_array_0(5)*3600 + time_array_0(6)*60 + time_array_0(7) + time_array_0(8)*0.001

  !--Set up the constants module used here
  call init_constants_derived
  call init_constants(.true.)    !initialize regional constants
  mindat=ione
  diagprint=.false.
  diagverbose=0

  !----READ NAMELIST----
  read(*,drw)
  !----READ NAMELIST----



  gatespc=gatespc*ithin
  numgates=int(maxobrange/gatespc) !calculate num gates based on namelist settings.
  nesteddata3d=trim(datapath) // trim(nesteddata3d) !concat strings
  nesteddata2d=trim(datapath) // trim(nesteddata2d)
  nestedgrid  =trim(datapath) // trim(nestedgrid)

  if(diagprint .and. diagverbose >= 1) then
     write(*,*) "----NC FILE SPECS--------"
     write(*,*) "nlat",nlat
     write(*,*) "nlon",nlon
     write(*,*) "nsig",nsig
     write(*,*) "grid spacing = 3km"

     write(*,*) "----Namelist settings----"
     write(*,*) "staid       ",staid
     write(*,*) "tilts       ",tilts
     write(*,*) "azimuths    ",azimuths
     write(*,*) "numgates    ",numgates
     write(*,*) "ithin       ",ithin
     write(*,*) "gatespc     ",gatespc
     write(*,*) "maxobrange  ",maxobrange
     write(*,*) "minobrange  ",minobrange
     write(*,*) "mintilt     ",mintilt
     write(*,*) "maxtilt     ",maxtilt
     write(*,*) "mindbz      ",mindbz
     write(*,*) "datapath:",datapath
     write(*,*) "nesteddata3d:",nesteddata3d
     write(*,*) "nesteddata2d:",nesteddata2d
     write(*,*) "radarcsv:",radarcsv
     write(*,*) "nestedgrid:",nestedgrid
     write(*,*) "diagprint: ",diagprint
     write(*,*) "diagverbose: ",diagverbose
  end if

  if(diagprint .and. diagverbose >= 3) then
     write(*,*) "deg2rad",deg2rad
     write(*,*) "rearth",rearth
     if(deg2rad>0) write(*,*) "(derived constants succesfully initialized)"
     if(rearth==6370.e03_r_kind) write(*,*) "(regional consts sucess init)"
  end if

  if (minobrange >= maxobrange) then
     STOP 'MINIMUM OB RANGE >= MAXIMUM OB RANGE. PROGRAM STOPPING NOW drwsim.f90'
  end if
  if (ithin < 1) then
     STOP 'ithin MUST BE >=1! CHECK NAMELIST AND RESET.'
  end if


!------ OPEN GLOBAL RADAR LIST ---------------------------
  write(*,*) "Reading Global Radar List"
  numradars=154
  allocate(dfid(numradars),dflat(numradars),dflon(numradars),dfheight(numradars))
  open(40,file=trim(radarcsv))
  read(40,'(a2,1x,a3,1x,a3,1x,a6)') cdummy !read 1st line which is just a header.
  do ii=1,numradars
     read(40,'(a12,1x,2f12.4,1x,f6.2)') dfid(ii),dflat(ii),dflon(ii),dfheight(ii)
     dfid(ii)=trim(dfid(ii))
     if(diagprint .and. diagverbose >= 3) write(*,*) dfid(ii),dflat(ii),dflon(ii),dfheight(ii)
  end do
  close(40)

!------ OPEN 3D NETCDF FILE ------------------------------
  ncid3d=30
  write(*,*) "Opening and reading nesteddata3d"
  ier=nf90_open(nesteddata3d,nf90_NoWrite,ncid3d)
  if(ier /= nf90_NoErr) STOP "ERROR READ NETCDF STOP!"

  ! time
  allocate(time(ntime))
  write(*,*) "Reading time"
  pcoord=0.0_r_kind !initialize
  ier=nf90_inq_varid(ncid3d,"time",tVarId)
  if(ier /= nf90_NoErr) STOP "ERROR INQ time NETCDF STOP!"
  ier=nf90_get_var(ncid3d,tVarId,time)
  if(ier /= nf90_NoErr) STOP "ERROR GET time NETCDF STOP!"

  ! pcoord
  allocate(pcoord(nsig))
  write(*,*) "Reading pcoord"
  pcoord=0.0_r_kind !initialize
  ier=nf90_inq_varid(ncid3d,"pfull",pVarId)
  if(ier /= nf90_NoErr) STOP "ERROR INQ pcoord NETCDF STOP!"
  ier=nf90_get_var(ncid3d,pVarId,pcoord)
  if(ier /= nf90_NoErr) STOP "ERROR GET pcoord NETCDF STOP!"
  !--convert from mb to Pa
  pcoord=pcoord*100.0_r_kind

  ! u
  allocate(ges_u(nlon,nlat,nsig,ntime))
  write(*,*) "Reading ges_u"
  ges_u=0.0_r_kind !initialize
  ier=nf90_inq_varid(ncid3d,"ucomp",uVarId)
  if(ier /= nf90_NoErr) STOP "ERROR INQ ges_u NETCDF STOP!"
  ier=nf90_get_var(ncid3d,uVarId,ges_u)
  if(ier /= nf90_NoErr) STOP "ERROR GET ges_u NETCDF STOP!"

  ! v
  allocate(ges_v(nlon,nlat,nsig,ntime))
  write(*,*) "Reading ges_v"
  ges_v=0.0_r_kind !initialize
  ier=nf90_inq_varid(ncid3d,"vcomp",vVarId)
  if(ier /= nf90_NoErr) STOP "ERROR INQ ges_v NETCDF STOP!"
  ier=nf90_get_var(ncid3d,vVarId,ges_v)
  if(ier /= nf90_NoErr) STOP "ERROR GET ges_v NETCDF STOP!"

  ! w
  allocate(ges_w(nlon,nlat,nsig,ntime))
  write(*,*) "Reading ges_w"
  ges_w=0.0_r_kind !initialize
  ier=nf90_inq_varid(ncid3d,"w",    wVarId)
  if(ier /= nf90_NoErr) STOP "ERROR INQ ges_w NETCDF STOP!"
  ier=nf90_get_var(ncid3d,wVarId,ges_w)
  if(ier /= nf90_NoErr) STOP "ERROR GET ges_w NETCDF STOP!"

  ! geop_hgtl - note this field may not be the correct one to use...
  allocate(geop_hgtl(nlon,nlat,nsig,ntime))
  write(*,*) "Reading geop_hgtl"
  geop_hgtl=0.0_r_kind !initialize
  ier=nf90_inq_varid(ncid3d,"delz", ghVarId)
  if(ier /= nf90_NoErr) STOP "ERROR INQ geop_hgtl NETCDF STOP!"
  ier=nf90_get_var(ncid3d,ghVarId,geop_hgtl)
  if(ier /= nf90_NoErr) STOP "ERROR GET geop_hgtl NETCDF STOP!"
!  geop_hgtl=geop_hgtl*grav

!  ! dbz - which doesn't exist in the file yet...
!  allocate(ges_dbz(nlon,nlat,nsig,ntime))
!  write(*,*) "Reading ges_dbz"
!  ges_dbz=0.0_r_kind !initialize
!  ier=nf90_inq_varid(ncid3d,"reflectivity", dbzVarId)
!  if(ier /= nf90_NoErr) STOP "ERROR INQ ges_dbz NETCDF STOP!"
!  ier=nf90_get_var(ncid3d,dbzVarId,ges_dbz)
!  if(ier /= nf90_NoErr) STOP "ERROR GET ges_dbz NETCDF STOP!"

  ! close file
  write(*,*) "Closing nesteddata3d"
  ier=nf90_close(ncid3d)
  if(ier /= nf90_NoErr) STOP "ERROR CLOSE 3D NETCDF STOP!"

!------ OPEN 2D NETCDF FILE ------------------------------
  ncid2d=20
  write(*,*) "Opening and reading nesteddata2d"
  ier=nf90_open(nesteddata2d,nf90_NoWrite,ncid2d)
  if(ier /= nf90_NoErr) STOP "ERROR READ NETCDF STOP!"

  ! z - note this field may not be the correct one to use...
  allocate(ges_z(nlon,nlat))
  write(*,*) "Reading ges_z"
  ges_z=0.0_r_kind !initialize
  ier=nf90_inq_varid(ncid2d,"HGTsfc", hgtVarID)
  if(ier /= nf90_NoErr) STOP "ERROR INQ hgtsfc NETCDF STOP!"
  ier=nf90_get_var(ncid2d,hgtVarID,ges_z)
  if(ier /= nf90_NoErr) STOP "ERROR GET hgtsfc NETCDF STOP!"
  !ges_z=ges_z/grav

  ! psfc
  allocate(ges_psfc(nlon,nlat))
  write(*,*) "Reading ges_psfc"
  ges_psfc=0.0_r_kind !initialize
  ier=nf90_inq_varid(ncid2d,"PRESsfc", psfcVarID)
  if(ier /= nf90_NoErr) STOP "ERROR INQ psfc NETCDF STOP!"
  ier=nf90_get_var(ncid2d,psfcVarID,ges_psfc)
  if(ier /= nf90_NoErr) STOP "ERROR GET psfc NETCDF STOP!"

  ! close file
  write(*,*) "Closing nesteddata2d"
  ier=nf90_close(ncid2d)
  if(ier /= nf90_NoErr) STOP "ERROR CLOSE 2D NETCDF STOP!"

!------ OPEN GRID SPEC NETCDF FILE -----------------------
  ncidgs=21
  write(*,*) "Opening and reading nestedgrid"
  ier=nf90_open(nestedgrid,nf90_NoWrite,ncidgs)
  if(ier /= nf90_NoErr) STOP "ERROR READ NETCDF STOP!"
  
  ! glon
  write(*,*) "Reading glon"
  allocate(glon(nlon,nlat))
  ier=nf90_inq_varid(ncidgs,"grid_lont", glonVarId)
  if(ier /= nf90_NoErr) STOP "ERROR INQ glon NETCDF STOP!"
  ier=nf90_get_var(ncidgs,glonVarId,glon)
  if(ier /= nf90_NoErr) STOP "ERROR GET glon NETCDF STOP!"

  ! glat
  write(*,*) "Reading glat"
  allocate(glat(nlon,nlat))
  ier=nf90_inq_varid(ncidgs,"grid_latt", glatVarId)
  if(ier /= nf90_NoErr) STOP "ERROR INQ glat NETCDF STOP!"
  ier=nf90_get_var(ncidgs,glatVarId,glat)
  if(ier /= nf90_NoErr) STOP "ERROR GET glat NETCDF STOP!"

  ! close file
  write(*,*) "Closing nestedgrid"
  ier=nf90_close(ncidgs)
  if(ier /= nf90_NoErr) STOP "ERROR CLOSE GRID SPEC NETCDF STOP!"
 

  !--I think you need to add each "delz" up to the get height at the associated
  !pressure level. Also, start from the bottom (level=63) and work your way up.
  allocate(height(nlon,nlat),delz(nlon,nlat)) 
  do itime=1,ntime
     height = zero
     do isig=nsig,1,-1
        delz                      = geop_hgtl(:,:,isig,itime)
        geop_hgtl(:,:,isig,itime) = height + delz
        height                    = geop_hgtl(:,:,isig,itime)
     end do
  end do

  !--subtract terrain height
  do itime=1,ntime 
     do isig=1,nsig
        geop_hgtl(:,:,isig,itime) = geop_hgtl(:,:,isig,itime) - ges_z(:,:)
     end do
  end do

   if(diagprint .and. diagverbose >= 2) then
      write(*,*)"time     ", maxval(     time(      :)), minval(     time(      :)) 
      write(*,*)"pcoord   ", maxval(   pcoord(    :  )), minval(   pcoord(    :  )) 
      write(*,*)"ges_u    ", maxval(    ges_u(:,:,:,:)), minval(    ges_u(:,:,:,:)) 
      write(*,*)"ges_v    ", maxval(    ges_v(:,:,:,:)), minval(    ges_v(:,:,:,:))
      write(*,*)"ges_w    ", maxval(    ges_w(:,:,:,:)), minval(    ges_w(:,:,:,:)) 
      write(*,*)"geop_hgtl", maxval(geop_hgtl(:,:,:,:)), minval(geop_hgtl(:,:,:,:)) 
!      write(*,*)"ges_dbz  ", maxval(  ges_dbz(:,:,:,:)), minval(  ges_dbz(:,:,:,:)) 
      write(*,*)"glon     ", maxval(     glon(:,:    )), minval(     glon(:,:    ))
      write(*,*)"glat     ", maxval(     glat(:,:    )), minval(     glat(:,:    ))
      write(*,*)"ges_z    ", maxval(    ges_z(:,:    )), minval(    ges_z(:,:    )) 
      write(*,*)"ges_psfc ", maxval( ges_psfc(:,:    )), minval( ges_psfc(:,:    )) 
   end if

   !--Make sure we're not missing fields.
   if(maxval(     time(      :)) == 0 .and. minval(     time(      :)) == 0) STOP 'MISSING TIME' 
   if(maxval(   pcoord(    :  )) == 0 .and. minval(   pcoord(    :  )) == 0) STOP 'MISSING PCOORD' 
   if(maxval(    ges_u(:,:,:,:)) == 0 .and. minval(    ges_u(:,:,:,:)) == 0) STOP 'MISSING U'  
   if(maxval(    ges_v(:,:,:,:)) == 0 .and. minval(    ges_v(:,:,:,:)) == 0) STOP 'MISSING V' 
   if(maxval(    ges_w(:,:,:,:)) == 0 .and. minval(    ges_w(:,:,:,:)) == 0) STOP 'MISSING W' 
   if(maxval(geop_hgtl(:,:,:,:)) == 0 .and. minval(geop_hgtl(:,:,:,:)) == 0) STOP 'MISSING GEOPOTHGT' 
!   if(maxval(  ges_dbz(:,:,:,:)) == 0 .and. minval(  ges_dbz(:,:,:,:)) == 0) STOP 'MISSING DBZ' 
   if(maxval(     glon(:,:    )) == 0 .and. minval(     glon(:,:    )) == 0) STOP 'MISSING GLON' 
   if(maxval(     glat(:,:    )) == 0 .and. minval(     glat(:,:    )) == 0) STOP 'MISSING GLAT' 
   if(maxval(    ges_z(:,:    )) == 0 .and. minval(    ges_z(:,:    )) == 0) STOP 'MISSING HGTSFC' 
   if(maxval( ges_psfc(:,:    )) == 0 .and. minval( ges_psfc(:,:    )) == 0) STOP 'MISSING PRESsfc'
   
  !set up information needed to interpolate model to observation
  !*********************************
  call gridmod_extract
  !*********************************

   if (diagprint .and. diagverbose >= 3) then
     write(6,*)'----Model grid info diagnostic print----'
     write(6,*)'nlon,nlat,rlambda0,pihalf,sign_pole:',nlon,nlat,rlambda0,pihalf,sign_pole
     write(6,*)'atilde_x,atilde_y,btilde_x,btilde_y:',atilde_x,atilde_y,btilde_x,btilde_y
     write(6,*)'rlon_min_dd,rlon_max_dd,rlat_min_dd,rlat_max_dd,nxtilde,nytilde:',&
                rlon_min_dd,rlon_max_dd,rlat_min_dd,rlat_max_dd,nxtilde,nytilde
     write(6,*)'Max/min i0_tilde,j0_tilde:',maxval(i0_tilde),minval(i0_tilde),maxval(j0_tilde),minval(j0_tilde)
     write(6,*)'Max/min ip_tilde,jp_tilde:',maxval(ip_tilde),minval(ip_tilde),maxval(jp_tilde),minval(jp_tilde)
     write(6,*)'Max/min xtilde0,ytilde0:',maxval(xtilde0),minval(xtilde0),maxval(ytilde0),minval(ytilde0)
     write(6,*)'----End model grid info diagnostic print----'
  end if


  ndata=izero

  !-Obtain analysis time in minutes since reference date
  !call w3fs21(iadate,mins_an)  !mins_an -integer number of mins snce 01/01/1978
  !rmins_an=mins_an             !convert to real number

  loopOVERtime: do itime=1,ntime
     if(itime > 1) call newdate(iadate,1,iadate) ! don't increment the first time.
     loopOVERradars: do irid=1,numradars 
        allocate(drwpol(nelv,360,numgates))
        drwpol=-999.0_r_kind !Initialize/Reset the drw polar field
        radar_location=.true. ! logical to only compute radar x,y once later in loop.
        this_staid=adjustl(trim(dfid(irid)))
        ifKGRK: if(this_staid==trim(adjustl(staid))) then
           stahgt=dfheight(irid)
           rlon0=dflon(irid)*deg2rad
           rlat0=dflat(irid)*deg2rad
           clat0=cos(rlat0)
           slat0=sin(rlat0)
           loopOVERtilts:    do itilt=1,nelv
              bufrisopen=.false.   !Initialize bufr file as closed.
              thistilt=tilts(itilt)
              thistiltr=thistilt*deg2rad
              celev0=cos(thistiltr)
              selev0=sin(thistiltr)
              loopOVERazimuths: do iazm=1,360
                 1000 format(a5,1x,i4,i2.2,i2.2,i2.2,4x,a6,1x,a4,4x,a5,1x,f4.1,a1,f4.1,4x,a4,1x,i3,a8)
                 write(*,1000),"Date:",iadate(1),iadate(2),iadate(3),iadate(4),&
                               "Radar:",adjustl(trim(dfid(irid))),&
                               "Tilt:",tilts(itilt),"/",tilts(nelv),&
                               "Azm:",iazm,"/360-deg"
                 thisazimuth=90.0_r_kind-float(iazm) ! 90-azm to be consistent with l2rwbufr
                 if(thisazimuth>=r360) thisazimuth=thisazimuth-r360
                 if(thisazimuth<zero) thisazimuth=thisazimuth+r360
                 thisazimuthr=thisazimuth*deg2rad
                 loopOVERgates: do igate=1,numgates     
                    inside=.false.
                    if(igate*gatespc >= minobrange .and. igate*gatespc <= maxobrange) inside=.true.
                    ifinside: if(inside) then
                       !--Find observation height using method from read_l2bufr_mod.f90 
                       thisrange=igate*gatespc
                       aactual=(rearth+stahgt)
                       a43=aactual*four_thirds
                       b   = thisrange*(thisrange+two*aactual*selev0)
                       c   = sqrt(aactual*aactual+b)
                       ha  = b/(aactual+c)
                       epsh=(thisrange*thisrange-ha*ha)/(r8*aactual)
                       h=ha-epsh
                       thishgt=stahgt+h
                       dpres=thishgt !store the absolute ob height (m) in dpres.                   
                       !--Find observation location using method fromread_l2bufr_mod.f90
                       !-Get corrected tilt angle
                       celev=celev0
                       selev=selev0
                       celev=a43*celev0/(a43+h)
                       selev=(thisrange*thisrange+h*h+two*a43*h)/(two*thisrange*(a43+h))
                       gamma=half*thisrange*(celev0+celev)
                       !-Get earth lat lon of ob
                       rad_per_meter=one/rearth
                       rlonloc=rad_per_meter*gamma*cos(thisazimuthr)
                       rlatloc=rad_per_meter*gamma*sin(thisazimuthr)
                       call invtllv(rlonloc,rlatloc,rlon0,clat0,slat0,rlonglob,rlatglob)
                       !-Find grid relative location of the radar.
                       if(radar_location) then
                          radar_lat=dflat(irid) !lat/lons stored as deg.
                          radar_lon=dflon(irid)
                          if(radar_lon>=r360) radar_lon=radar_lon-r360 !fix if needed.
                          if(radar_lon<zero) radar_lon=radar_lon+r360
                          radar_lon=radar_lon*deg2rad !convert to radians.
                          radar_lat=radar_lat*deg2rad
                          call tll2xy(radar_lon,radar_lat,radar_x,radar_y)
                          if(diagprint .and. diagverbose >= 3) write(*,*) "Radar x,y location is:",radar_x,radar_y 
                          if(diagprint .and. diagverbose >= 3) write(*,*) "Radar lon,lat is     :",radar_lon*rad2deg,radar_lat*rad2deg
                          radar_location=.false. ! turn off get radar x/y until next radar is processed.
                       end if
                       !-Find grid relative location of the ob.
                       thislat=rlatglob*rad2deg
                       thislon=rlonglob*rad2deg
                       if(thislon>=r360) thislon=thislon-r360
                       if(thislon<zero) thislon=thislon+r360
                       thislat=thislat*deg2rad
                       thislon=thislon*deg2rad
                       call tll2xy(thislon,thislat,dlon,dlat)
                       if(diagprint .and. diagverbose >= 5) write(*,*) "ob x,y location is   :",dlon,dlat
                       if(diagprint .and. diagverbose >= 5) write(*,*) "ob lon,lat is        :",thislon*rad2deg,thislat*rad2deg

!**********************************DO WIND ROTATION IF RADIAL WIND****!
!   SHOULD NOT BE NEEDED FOR FV3???
!**********************************DO WIND ROTATION IF RADIAL WIND****!
                    
                       call tintrp2a_single_level(ges_z,zsges,dlon,dlat)
                       call tintrp2a_single_level(ges_psfc,psfcsges,dlon,dlat)
                       dpres=dpres-zsges      !  remove terrain height from ob absolute height
                       call tintrp2a(geop_hgtl(:,:,:,itime),hges,dlon,dlat,nsig)  
!                  !    Convert geopotential height at layer midpoints to
!                  !    geometric height using
!                  !    equations (17, 20, 23) in MJ Mahoney's note "A
!                  !    discussion of various
!                  !    measures of altitude" (2001).  Available on the web at
!                  !    http://mtp.jpl.nasa.gov/notes/altitude/altitude.html
!                  !    http://archive.is/reKhz
!                  !
!                  !    termg  = equation 17
!                  !    termr  = equation 21
!                  !    termrg = first term in the denominator of equation 23
!                  !    zges   = equation 23
                       sin2  = sin(thislat)*sin(thislat)
                       termg = grav_equator * ((one+somigliana*sin2)/sqrt(one-eccentricity*eccentricity*sin2))
                       termr = semi_major_axis /(one + flattening + grav_ratio-  two*flattening*sin2)
                       termrg = (termg/grav)*termr
                       do n=1,nsig
                          zges(n) = (termr*hges(n)) / (termrg-hges(n))  ! eq (23)
                       end do
                       if(diagprint .and. diagverbose >=1 .and. iazm==360 .and. igate==20) then
                           write(*,*) "dpres   = ",dpres,igate
                           write(*,*) "thishgt = ",thishgt,igate
                           write(*,*) "zsges   = ",zsges,igate
                       end if
                  !    Convert observation height (in dpres) from meters to grid relative units.
                       !call grdcrd(dpres,1,zges,nsig,-1)
                       call grdcrd1(dpres,zges,nsig,-1)
                    
                       if(diagprint .and. diagverbose >=1 .and. iazm==360 .and. igate==20) then
                           write(*,*) "After grdcrd"
                           write(*,*) "dpres   = ",dpres,igate
                           write(*,*) "thishgt = ",thishgt,igate
                           write(*,*) "zsges   = ",zsges,igate
                           write(*,*) "zges       = ",zges
                           write(*,*) "pcoord     = ",pcoord
                           write(*,*) "ges_psfc   = ",psfcsges
                           !write(*,*) "geop_hgtl  = ",geop_hgtl(thislon,thislat,:,1)
                       end if
                     ! Interpolate guess dbz to observation location - cycle if below threshold.
!                       call tintrp3(ges_dbz(:,:,:,itime),dbzgesin,dlon,dlat,dpres)
!                       dbzCheck: if(dbzgesin >= mindbz .or. .true.) then
                     !    Interpolate guess wind to observation location                                  
                          call tintrp3(ges_u(:,:,:,itime),ugesin,dlon,dlat,dpres)
                          call tintrp3(ges_v(:,:,:,itime),vgesin,dlon,dlat,dpres) !should these be dlat,dlon?
                          call tintrp3(ges_w(:,:,:,itime),wgesin,dlon,dlat,dpres)
                     !    Convert guess u,v,w wind components to radial value
                          cosazm  = cos(thisazimuthr)! cos(azimuth angle)                       
                          sinazm  = sin(thisazimuthr)! sin(azimuth angle)                       
                          costilt = cos(thistiltr)   ! cos(tilt angle)
                          sintilt = sin(thistiltr)   ! sin(tilt angle)
                     !-------------WIND FORWARD MODEL-----------------------------------------!   
                          iazm90=90-iazm
                          if(iazm90>=r360) iazm90=iazm90-r360
                          if(iazm90< zero) iazm90=iazm90+r360
                          drwpol(itilt,iazm,igate) = ugesin*cosazm*costilt  +vgesin*sinazm*costilt  +wgesin*sintilt 
                          ndata=ndata+1
!                       end if dbzCheck
                    end if ifinside
                 end do loopOVERgates
              end do loopOVERazimuths
           end do loopOVERtilts

           if(diagprint .and. diagverbose >= 2) print*,minval(drwpol),maxval(drwpol)


           !-------------BUFFERIZE--------------------------------------------------!
           !    At this point we have observations from all radars at every scan
           ! angle at a single time. We will put this information in its own bufr
           ! file hence this is contained within loopOVERtime.
           !
           write(*,*)"Writing bufr file for ",trim(dfid(irid))
           hdstr='SSTN CLON CLAT SELV ANEL YEAR MNTH DAYS HOUR MINU QCRW ANAZ'
           obstr='DIST125M DMVR DVSW'                     !NL2RW--level 2 radial wind.
           open(41,file='l2rwbufr.table.csv')        
           read(41,'(a10)') cdummy !read 1st line which is just a header.
           do ii=0,23 !00z to 23z -- this starts on the second line of the file.
              if(ii<iadate(4) .or. ii>iadate(4)) then
                  read(41,'(a10)') cdummy
              else if(ii==iadate(4)) then
                  read(41,'(a10)') message_type
                  message_type=trim(message_type) 
              end if
           end do
           close(41)
           if(diagprint .and. diagverbose >= 2) write(*,*) message_type
           call w3ai15(iadate(1),yyyy,1,4,'')
           call w3ai15(iadate(2),  mm,1,2,'')
           call w3ai15(iadate(3),  dd,1,2,'')
           call w3ai15(iadate(4),  hh,1,2,'')
           call w3ai15(iadate(5),  mn,1,2,'')
           idate=trim(yyyy)//trim(mm)//trim(dd)//trim(hh)
           subset=trim(adjustl(message_type))
           chdr   = dfid(irid)       !SSTN - RADAR STATION IDENTIFIER -- uses same memory location as hdr(1)
           hdr(2) = dflon(irid)      !CLON - LONGITUDE (COARSE ACCURACY)
           hdr(3) = dflat(irid)      !CLAT - LATITUDE (COARSE ACCURACY)
           hdr(4) = dfheight(irid)   !SELV - HEIGHT OF STATION
          !hdr(5) - tilt loop below.
           hdr(6) = iadate(1)        !YEAR - YEAR
           hdr(7) = iadate(2)        !MNTH - MONTH
           hdr(8) = iadate(3)        !DAYS - DAY
           hdr(9) = iadate(4)        !HOUR - HOUR 
           hdr(10)= 00               !MINU - MINUTE
           hdr(11)= 1                !QCRW - QUALITY MARK FOR WINDS ALONG RADIAL LINE
          !hdr(12)- azm loop below.
           bufrfilename=trim(idate)//'_fv3.t'//trim(hh)//'z.drw.bufr'
           bufrtilt: do itiltbufr=1,nelv
              intdate=iadate(1)*1000000 + iadate(2)*10000 + iadate(3)*100 + iadate(4) ! int(yyyymmddhh)
              hdr(5) = tilts(itiltbufr) 
              if(.not.bufrisopen) then !open a new message for each tilt 
                 open(unit=10,file=trim(bufrfilename),status='unknown',action='write',form='unformatted')
                 open(unit=11,file='l2rwbufr.table',status='old',action='read',form='formatted')
                 call openbf(10,'OUT',11)
                 bufrisopen=.true.
              end if
              call openmb(10,trim(subset),intdate)
              bufrazm: do iazmbufr=1,360
                 iazmbufr90=90-iazmbufr
                 if(iazmbufr90>=r360) iazmbufr90=iazmbufr90-r360
                 if(iazmbufr90< zero) iazmbufr90=iazmbufr90+r360
                 allocate(obs(3,numgates))
                 obs=-999.0_r_kind ! Initialize as missing values
                 hdr(12)=float(iazmbufr90)
                 bufrgate: do igatebufr=1,numgates
                    obs(1,igatebufr) = igatebufr !DISTANCE (FROM ANTENNA TO GATE CENTER) IN UNITS OF 250M
                    obs(2,igatebufr) = drwpol(itiltbufr,iazmbufr90,igatebufr) !DOPPLER MEAN RADIAL VELOC 
                    obs(3,igatebufr) = 1.0_r_kind                       !DOPPLER VELOCITY SPECTRAL WIDTH
                 end do bufrgate
                 ! encode radial velocity
                 call ufbint(10,hdr,12,1,iret,trim(hdstr))
                 call ufbint(10,obs, 3,numgates,iret,trim(obstr))
                 call writsb(10)
                 deallocate(obs)
              end do bufrazm
              call closmg(10) ! close bufr message
           end do bufrtilt
           call closbf(10) ! close bufr file
           close(10)       ! close bufr file
           close(11)       ! close l2rwbufr.table
           bufrisopen=.false.
        
        end if ifKGRK
        deallocate(drwpol)
     end do loopOVERradars 
  end do loopOVERtime

  !-Call Timer
  call date_and_time(values=time_array_1)
  finish_time = time_array_1(5)*3600 + time_array_1(6)*60 + time_array_1(7) + time_array_1(8)*0.001
  total_time=finish_time-start_time
  hrs =int(        total_time/3600.0       )
  mins=int(    mod(total_time,3600.0)/60.0 )
  secs=int(mod(mod(total_time,3600.0),60.0))
  write(*,'(a14,i2.2,a1,i2.2,a1,i2.2)')"Elapsed time:   ",hrs,":",mins,":",secs
  write(*,*) "Elapsed time (s) =", total_time
  write(*,*) "end of program"
  write(*,*) 'numgates =',numgates
  write(*,*) 'numelevs =',nelv
  write(*,*) 'numazims =',azimuths
  write(*,*) 'numtimes =',ntime
  write(*,*) 'ndata    =',ndata
end program drwsim


  subroutine tll2xy(rlon,rlat,x,y)

! !USES:

    use constants, only: one
    use interp_util
    implicit none

    real(r_kind),intent(in   ) :: rlon  ! earth longitude (radians)
    real(r_kind),intent(in   ) :: rlat  ! earth latitude  (radians)

! !OUTPUT PARAMETERS:

    real(r_kind),intent(  out) :: x  ! x-grid coordinate (grid units)
    real(r_kind),intent(  out) :: y  ! y-grid coordinate (grid units)
!    logical     ,intent(  out) :: outside     ! .false., then point is inside x-y domain
                                              ! .true.,  then point is outside
                                              ! x-y domain

    real(r_kind) clon,slon,r_of_lat,xtilde,ytilde
    real(r_kind) dtilde,etilde
    real(r_kind) d1tilde,d2tilde,e1tilde,e2tilde,detinv
    integer(i_kind) itilde,jtilde
    integer(i_kind) i0,j0,ip,jp

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

!    outside=x < rlon_min_dd .or. x > rlon_max_dd .or. &
!            y < rlat_min_dd .or. y > rlat_max_dd

 end subroutine tll2xy



SUBROUTINE invtllv(ALM,APH,TLMO,CTPH0,STPH0,TLM,TPH)
  use kinds, only:  r_kind
  use interp_util
  implicit none

  real(r_kind),intent(in   ) :: alm,aph,tlmo,ctph0,stph0
  real(r_kind),intent(  out) :: tlm,tph

  real(r_kind):: relm,srlm,crlm,sph,cph,cc,anum,denom

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

subroutine newdate(indate,nhr,outdate)
! This routine takes a date and an amount of time in hours to increment and
! outputs the new date. indate must be dimension 4 with years starting in 1st
! dimension.
! input:   indate (yyyymmddhh) 
!             nhr (        hh)
! output: outdate (yyyymmddhh)

  use kinds

  implicit none

  integer(i_kind),dimension(4),intent(in   ) :: indate
  integer(i_kind)             ,intent(in   ) :: nhr
  integer(i_kind),dimension(4),intent(  out) :: outdate

  !--local declarations
  integer(i_kind) :: yyyy,mm,dd,hh,maxdd

  yyyy=indate(1)
  mm=indate(2)
  dd=indate(3)
  hh=indate(4)

  if(mm==01 .or. mm==03 .or. mm==05 .or. mm==07 .or. mm==08 .or. mm==10 .or. mm==12) then
     maxdd=31
  else if(mm==02 .and. mod(yyyy,4)==0) then
     maxdd=28
  else if(mm==02 .and. mod(yyyy,4)> 0) then
     maxdd=29
  else if(mm==04 .or. mm==06 .or. mm==09 .or. mm==11) then
     maxdd=30
  end if
  
  hh=hh+nhr
  do while (hh >= 24)
     if(hh >= 24) then
        hh=hh-24
        dd=dd+1
        if(dd>maxdd) then
          dd=01
          mm=mm+1
          if(mm>12) then
             mm=01
             yyyy=yyyy+1
          end if
        end if
     end if
  end do
  
  outdate(1)=yyyy
  outdate(2)=mm
  outdate(3)=dd
  outdate(4)=hh

end subroutine
