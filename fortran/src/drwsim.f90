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
                     ikx,mins_an,mins_ob,lu,n,iret
  integer(i_kind) :: maxobs,nchanl,ilat,ilon,maxgate,rad_nelv,cm

  integer(i_kind),dimension(5) :: obdate,iadate,intdate
  character(4) :: yyyy
  character(2) :: mm,dd,hh,mn,fhr
  character(10)             :: idate
  real(r_kind),dimension(nsig) :: zges,hges
  real(r_kind) :: azm,cosazm,sinazm,costilt,sintilt,cosazm_earth,sinazm_earth
  real(r_kind) :: ugesin,vgesin,wgesin
  real(r_kind) :: zsges,ddiff,observation,rms,bias,sumrms,sumbias
  real(r_kind) :: sin2,termg,termr,termrg,zob,dpres
  real(r_kind) :: a,b,c,ha,epsh,h,aactual,a43,thistilt,x
  real(r_kind) :: thistiltr,selev0,celev0,thisrange,this_stahgt,thishgt
  real(r_kind) :: celev,selev,gamma,thisazimuthr,thisazimuth,rlon0,rlat0,stahgt, &
                  clat0,slat0,dlat,dlon,thiserr,thislon,thislat, &
                  rlonloc,rlatloc,rlonglob,rlatglob,timeb,rad_per_meter, &
                  delta_gate,delta_az,radar_x,radar_y,radar_lon,radar_lat
  real(r_kind) :: fcstgesin
  real(r_kind),allocatable :: drwpol(:,:,:) !tilt,azm,gate
  real(r_kind) :: radar_twindow                                          
  real(r_kind) :: rmins_an
  real(r_kind) :: rdummy
  integer(i_kind):: idummy
  integer(i_kind) :: irid,itilt,iazm,igate,iazm90


  character(8) cstaid
  character(4) this_staid
  logical   :: outside,diagprint,inside,bufrisopen,endbufr
  integer   :: diagverbose

  type(radar),allocatable :: strct_in_vel(:,:),rad(:)


  !---------SETTINGS---------!
  integer(i_kind) :: maxobrange=250000_i_kind  ! Range (m) *within* which to use observations 
  integer(i_kind) :: minobrange=20000_i_kind   ! Range (m) *outside* of which
  real(r_kind)    :: mintilt=0.0_r_kind        ! Only use tilt(elevation) angles (deg) >= this number 
  real(r_kind)    :: maxtilt=5.5_r_kind        ! Do no use tilt(elevation) angles (deg) <= this number
  integer(i_kind) :: ithin=4_i_kind            ! Gates to skip for ob thinning (must be >=1)
  real(r_kind)    :: radartwindow=2.5_r_kind   ! Time window within which to grid observations   (minutes)
  character(4)    :: staid='KOUN'
  real(r_kind)    :: mindbz=0_r_kind
  real(r_kind),dimension(25)    :: tilts=0_r_kind
  integer(i_kind) :: azimuths=360_i_kind
  integer(i_kind) :: gatespc=250_i_kind
  integer(i_kind) :: numgates=400_i_kind
  !----------------------------------------------!

  !---------NETCDF VARS---------!
  integer(i_kind) :: ncid3d,ncid2d,ncidgs,ier
  integer(i_kind) :: pVarId,uVarId,vVarId,wVarId,ghVarId !3d
  integer(i_kind) :: dbzVarId                            !2d
  integer(i_kind) :: glonVarId,glatVarId                 !gs
  integer(i_kind) :: numsig,numlat,numlon,numtimes
  real(r_kind),allocatable  ::    pcoord(    :  ) ! (nsig)
  real(r_kind),allocatable  ::     ges_u(:,:,:,:) ! (nlon,nlat,nsig,ntime)
  real(r_kind),allocatable  ::     ges_v(:,:,:,:) ! (nlon,nlat,nsig,ntime)
  real(r_kind),allocatable  ::     ges_w(:,:,:,:) ! (nlon,nlat,nsig,ntime)
  real(r_kind),allocatable  :: geop_hgtl(:,:,:,:) ! (nlon,nlat,nsig,ntime)
  real(r_kind),allocatable  ::     ges_z(:,:   ) ! (nlon,nlat,          )
  real(r_kind),allocatable  ::   ges_dbz(:,:,  :) ! (nlon,nlat,    ,ntime)
  real(r_kind),allocatable  ::      lons(:,:    ) ! (nlon,nlat           )
  real(r_kind),allocatable  ::      lats(:,:    ) ! (nlon,nlat           )
  !----------------------------------------------!

  !---------GLOBAL RADAR CSV FILE VARS---------!
  integer(i_kind)                         :: ii,numradars
  real(r_kind),dimension(:),allocatable   :: dflat,dflon,dfheight
  character(12),dimension(:),allocatable  :: dfid
  !----------------------------------------------!


  !---------BUFR VARS--------------------------!
  integer(i_kind) :: itiltbufr,iazmbufr,igatebufr
  character(80) :: bufrfilename,hdstr,obstr
  real(r_double) :: hdr(12)
  real(r_double),allocatable :: obs(:,:)
  character(8) :: chdr,subset
  equivalence (hdr(1),chdr)

  !---------L2RWBUFR CSV TABLE VARS------------!
  character(10) :: message_type
  character(8) :: cdummy



  character(len=120)  :: datapath
  character(len=120)  :: nesteddata3d,nesteddata2d,nestedgrid,radarcsv
  character(len=120)  :: bufroutfile

  namelist/drw/iadate,radartwindow,nesteddata3d,nesteddata2d,nestedgrid,bufroutfile,&
               mintilt,maxtilt,staid,mindbz,tilts,maxobrange,minobrange,azimuths,ithin,&
               gatespc,datapath,diagprint,diagverbose,radarcsv


!--------------------------------------------------------------------------------------!
!                            END OF ALL DECLARATIONS
!                            !
!--------------------------------------------------------------------------------------!

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
     print *, "----NC FILE SPECS--------"
     print *, "nlat",nlat
     print *, "nlon",nlon
     print *, "nsig",nsig
     print *, "grid spacing = 3km"

     print *, "----Namelist settings----"
     print *, "staid       ",staid
     print *, "tilts       ",tilts
     print *, "azimuths    ",azimuths
     print *, "numgates    ",numgates
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
     print *, "radarcsv:",radarcsv
     print *, "nestedgrid:",nestedgrid
     print *, "diagprint: ",diagprint
     print *, "diagverbose: ",diagverbose
  end if

  if(diagprint .and. diagverbose >= 3) then
     print *, "deg2rad",deg2rad
     print *, "rearth",rearth
     if(deg2rad>0) print *, "(derived constants succesfully initialized)"
     if(rearth==6370.e03_r_kind) print *, "(regional consts sucess init)"
  end if

  if (minobrange >= maxobrange) then
     STOP 'MINIMUM OB RANGE >= MAXIMUM OB RANGE. PROGRAM STOPPING NOW drwsim.f90'
  end if
  if (ithin < 1) then
     STOP 'ithin MUST BE >=1! CHECK NAMELIST AND RESET.'
  end if


!------ OPEN GLOBAL RADAR LIST ---------------------------
  if(diagprint .and. diagverbose >= 1) print *, "Reading Global Radar List"
  numradars=154
  allocate(dfid(numradars),dflat(numradars),dflon(numradars),dfheight(numradars))
  open(40,file=trim(radarcsv))!,status='old',action='read',iostat=ierror,form='formatted')
  do ii=1,numradars
     read(40,'(a12,1x,2f12.4,1x,f6.2)') dfid(ii),dflat(ii),dflon(ii),dfheight(ii)
     dfid(ii)=trim(dfid(ii))
     if(diagprint .and. diagverbose >= 2) print *, dfid(ii),dflat(ii),dflon(ii),dfheight(ii)
  end do
  close(40)

!------ OPEN 3D NETCDF FILE ------------------------------
  ncid3d=30
  if(diagprint .and. diagverbose >= 1) print *, "Opening and reading nesteddata3d"
  ier=nf90_open(nesteddata3d,nf90_NoWrite,ncid3d)
  if(ier /= nf90_NoErr) STOP "ERROR READ NETCDF STOP!"

  ! pcoord
  allocate(pcoord(nsig))
  if(diagprint .and. diagverbose >= 2) print *, "Reading pcoord"
  pcoord=0.0_r_kind !initialize
  ier=nf90_inq_varid(ncid3d,"pfull",pVarId)
  if(ier /= nf90_NoErr) STOP "ERROR INQ pcoord NETCDF STOP!"
  ier=nf90_get_var(ncid3d,pVarId,pcoord)
  if(ier /= nf90_NoErr) STOP "ERROR GET pcoord NETCDF STOP!"

  ! u
  allocate(ges_u(nlon,nlat,nsig,ntime))
  if(diagprint .and. diagverbose >= 2) print *, "Reading ges_u"
  ges_u=0.0_r_kind !initialize
  ier=nf90_inq_varid(ncid3d,"ucomp",uVarId)
  if(ier /= nf90_NoErr) STOP "ERROR INQ ges_u NETCDF STOP!"
  ier=nf90_get_var(ncid3d,uVarId,ges_u)!,start=(/1,1,1,idummy/))
  if(ier /= nf90_NoErr) STOP "ERROR GET ges_u NETCDF STOP!"

  ! v
  allocate(ges_v(nlon,nlat,nsig,ntime))
  if(diagprint .and. diagverbose >= 2) print *, "Reading ges_v"
  ges_v=0.0_r_kind !initialize
  ier=nf90_inq_varid(ncid3d,"vcomp",vVarId)
  if(ier /= nf90_NoErr) STOP "ERROR INQ ges_v NETCDF STOP!"
  ier=nf90_get_var(ncid3d,vVarId,ges_v)
  if(ier /= nf90_NoErr) STOP "ERROR GET ges_v NETCDF STOP!"

  ! w
  allocate(ges_w(nlon,nlat,nsig,ntime))
  if(diagprint .and. diagverbose >= 2) print *, "Reading ges_w"
  ges_w=0.0_r_kind !initialize
  ier=nf90_inq_varid(ncid3d,"w",    wVarId)
  if(ier /= nf90_NoErr) STOP "ERROR INQ ges_w NETCDF STOP!"
  ier=nf90_get_var(ncid3d,wVarId,ges_w)
  if(ier /= nf90_NoErr) STOP "ERROR GET ges_w NETCDF STOP!"

  ! z
  allocate(ges_z(nlon,nlat))
  if(diagprint .and. diagverbose >= 2) print *, "Reading ges_z"
  ges_z=0.0_r_kind !initialize

  ! geop_hgtl
  allocate(geop_hgtl(nlon,nlat,nsig,ntime))
  if(diagprint .and. diagverbose >= 2) print *, "Reading geop_hgtl"
  geop_hgtl=0.0_r_kind !initialize
  ier=nf90_inq_varid(ncid3d,"delz", ghVarId)
  if(ier /= nf90_NoErr) STOP "ERROR INQ geop_hgtl NETCDF STOP!"
  ier=nf90_get_var(ncid3d,ghVarId,geop_hgtl)
  if(ier /= nf90_NoErr) STOP "ERROR GET geop_hgtl NETCDF STOP!"

  ! close file
  if(diagprint .and. diagverbose >= 1) print *, "Closing nesteddata3d"
  ier=nf90_close(ncid3d)
  if(ier /= nf90_NoErr) STOP "ERROR CLOSE 3D NETCDF STOP!"

!------ OPEN 2D NETCDF FILE ------------------------------
  ncid2d=20
  if(diagprint .and. diagverbose >= 1) print *, "Opening and reading nesteddata2d"
  ier=nf90_open(nesteddata2d,nf90_NoWrite,ncid2d)
  if(ier /= nf90_NoErr) STOP "ERROR READ NETCDF STOP!"

  ! dbz - which doesn't exist in the file yet...
  allocate(ges_dbz(nlon,nlat,ntime))
  if(diagprint .and. diagverbose >= 2) print *, "Reading ges_dbz"
  ges_dbz=0.0_r_kind !initialize
  ier=nf90_inq_varid(ncid2d,"PRATEsfc", dbzVarId)
  if(ier /= nf90_NoErr) STOP "ERROR INQ ges_dbz NETCDF STOP!"
  ier=nf90_get_var(ncid2d,dbzVarId,ges_dbz)
  if(ier /= nf90_NoErr) STOP "ERROR GET ges_dbz NETCDF STOP!"

  ! close file
  if(diagprint .and. diagverbose >= 1) print *, "Closing nesteddata2d"
  ier=nf90_close(ncid2d)
  if(ier /= nf90_NoErr) STOP "ERROR CLOSE 2D NETCDF STOP!"

!------ OPEN GRID SPEC NETCDF FILE -----------------------
  ncidgs=21
  if(diagprint .and. diagverbose >= 2) print *, "Opening and reading nestedgrid"
  ier=nf90_open(nestedgrid,nf90_NoWrite,ncidgs)
  if(ier /= nf90_NoErr) STOP "ERROR READ NETCDF STOP!"
  
  ! glon
  if(diagprint .and. diagverbose >= 2) print *, "Reading glon"
  allocate(glon(nlon,nlat))
  ier=nf90_inq_varid(ncidgs,"grid_lont", glonVarId)
  if(ier /= nf90_NoErr) STOP "ERROR INQ glon NETCDF STOP!"
  ier=nf90_get_var(ncidgs,glonVarId,glon)
  if(ier /= nf90_NoErr) STOP "ERROR GET glon NETCDF STOP!"

  ! glat
  if(diagprint .and. diagverbose >= 2) print *, "Reading glat"
  allocate(glat(nlon,nlat))
  ier=nf90_inq_varid(ncidgs,"grid_latt", glatVarId)
  if(ier /= nf90_NoErr) STOP "ERROR INQ glat NETCDF STOP!"
  ier=nf90_get_var(ncidgs,glatVarId,glat)
  if(ier /= nf90_NoErr) STOP "ERROR GET glat NETCDF STOP!"

  ! close file
  if(diagprint .and. diagverbose >= 1) print *, "Closing nestedgrid"
  ier=nf90_close(ncidgs)
  if(ier /= nf90_NoErr) STOP "ERROR CLOSE GRID SPEC NETCDF STOP!"
  

   if(diagprint .and. diagverbose >= 3) then
      print *, maxval(    ges_u(:,:,:,:)) 
      print *, maxval(    ges_v(:,:,:,:)) 
      print *, maxval(    ges_w(:,:,:,:)) 
      print *, maxval(geop_hgtl(:,:,:,:)) 
      print *, maxval(  ges_dbz(:,:,  :)) 
      print *, maxval(     glon(:,:    ))
      print *, maxval(     glat(:,:    ))
      print *, maxval(    ges_z(:,:    )) 
   end if
   
  !set up information needed to interpolate model to observation
  !*********************************
  ! I DON'T THINK WE NEED TO DO THIS SINCE THEY"RE ALREADY LON/LAT NOT LAT/LON
  call gridmod_extract
  if(diagprint .and. diagverbose>=3)  print *, "back from gridmod_extract"
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
  nread=izero

  !-Obtain analysis time in minutes since reference date

  call w3fs21(iadate,mins_an)  !mins_an -integer number of mins snce 01/01/1978
  rmins_an=mins_an             !convert to real number

 
  nelv=25
  rad_nelv=25 
  allocate(drwpol(1,360,numgates))

  loopOVERradars: do irid=1,numradars 
     drwpol=-999.0_r_kind !Initialize/Reset the drw polar field
     endbufr=.false.
     ifKGRK: if(trim(dfid(irid))==" KGRK") then ! The space before KGRK is needed...
        stahgt=dfheight(irid)
        rlon0=dflon(irid)*deg2rad
        rlat0=dflat(irid)*deg2rad
        clat0=cos(rlat0)
        slat0=sin(rlat0)
        loopOVERtilts:    do itilt=1,1!nelv
           bufrisopen=.false.   !Initialize bufr file as closed.
           thistilt=tilts(itilt)
           thistiltr=thistilt*deg2rad
           celev0=cos(thistiltr)
           selev0=sin(thistiltr)
           loopOVERazimuths: do iazm=1,360
              !iazm90=iazm-90.0_r_kind
              thisazimuthr=iazm*deg2rad
              loopOVERgates: do igate=1,numgates     
                 if(diagprint .and. diagverbose >= 2) print *, igate,iazm,itilt,irid
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
                    radar_lat=dflat(irid)
                    radar_lon=dflon(irid)
                    if(radar_lon>=360) radar_lon=radar_lon-360
                    if(radar_lon<zero) radar_lon=radar_lon+360
                    radar_lon=radar_lon*deg2rad
                    radar_lat=radar_lat*deg2rad
                    call tll2xy(radar_lon,radar_lat,radar_x,radar_y)
                    if(diagprint .and. diagverbose >= 2) print *, "Radar x,y location on the model grid is:",&
                                                                   radar_x,radar_y 
                    !-Find grid relative location of the ob.
                    !if(diagprint .and. diagverbose >= 2) print*,"one:",rlatglob,rlonglob
                    thislat=rlatglob*rad2deg
                    thislon=rlonglob*rad2deg
                    if(thislon>=360) thislon=thislon-360
                    if(thislon<zero) thislon=thislon+360
                    thislat=thislat*deg2rad
                    thislon=thislon*deg2rad
                    call tll2xy(thislon,thislat,dlon,dlat)
                    call tintrp2a_single_level(ges_z,zsges,dlat,dlon)
                    if(diagprint .and. diagverbose >= 3) print *,"Back from tintrp2a_single_level"
                    !Subtract off the terrain height
                    dpres=dpres-zsges      !  remove terrain height from ob absolute height
                    call tintrp2a(geop_hgtl,hges,dlat,dlon,nsig)  
                    if(diagprint .and. diagverbose >= 3) print *,"Back from tintrp2a"
               !    Convert geopotential height at layer midpoints to
               !    geometric height using
               !    equations (17, 20, 23) in MJ Mahoney's note "A
               !    discussion of various
               !    measures of altitude" (2001).  Available on the web at
               !    http://mtp.jpl.nasa.gov/notes/altitude/altitude.html
               !
               !    termg  = equation 17
               !    termr  = equation 21
               !    termrg = first term in the denominator of equation 23
               !    zges   = equation 23
                    sin2  = sin(thislat)*sin(thislat)
                    termg = grav_equator * &
                         ((one+somigliana*sin2)/sqrt(one-eccentricity*eccentricity*sin2))
                    termr = semi_major_axis /(one + flattening + grav_ratio-  &
                         two*flattening*sin2)
                    termrg = (termg/grav)*termr
                    do n=1,nsig
                       zges(n) = (termr*hges(n)) / (termrg-hges(n))  ! eq (23)
                    end do
               !    Convert observation height (in dpres) from meters to
               !    grid relative
               !    units.  Save the observation height in zob for later
               !    use.
                    zob = dpres
                    call grdcrd(dpres,1,zges,nsig,1)
                    if(diagprint .and. diagverbose >= 3) print *,"Back from tintrp2a"
               !    Interpolate guess u and v to observation location                                  
                    call tintrp3(ges_u,ugesin,dlat,dlon,dpres)
                    call tintrp3(ges_v,vgesin,dlat,dlon,dpres)
                    call tintrp3(ges_w,wgesin,dlat,dlon,dpres)
               !    Convert guess u,v wind components to radial value
               !    consident with obs    
                    cosazm  = cos(thisazimuthr)! cos(azimuth angle)                       
                    sinazm  = sin(thisazimuthr)! sin(azimuth angle)                       
                    costilt = cos(thistiltr)   ! cos(tilt angle)
                    sintilt = sin(thistiltr)   ! sin(tilt angle)
               !-------------WIND FORWARD MODEL-----------------------------------------!   
                    drwpol(itilt,iazm,igate) = ugesin*cosazm*costilt  +vgesin*sinazm*costilt  +wgesin*sintilt 
   
                 end if ifinside
              end do loopOVERgates
           end do loopOVERazimuths
        end do loopOVERtilts

        if(diagprint .and. diagverbose >= 2) print*,minval(drwpol),maxval(drwpol)


        !-------------BUFFERIZE--------------------------------------------------!
        if(diagprint .and. diagverbose >= 1) print *,"Writing bufr file for ",trim(dfid(irid))
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
        if(diagprint .and. diagverbose >= 2) print *, message_type
        call w3ai15(iadate(1),yyyy,1,4,'')
        call w3ai15(iadate(2),  mm,1,2,'')
        call w3ai15(iadate(3),  dd,1,2,'')
        call w3ai15(iadate(4),  hh,1,2,'')
        call w3ai15(iadate(5),  mn,1,2,'')
        fhr='00'
        idate=trim(yyyy)//trim(mm)//trim(dd)//trim(hh)
        subset=trim(adjustl(message_type))
        chdr   = dfid(irid)       !SSTN - RADAR STATION IDENTIFIER -- uses same memory location as hdr(1)
        hdr(2) = dflon(irid)      !CLON - LONGITUDE (COARSE ACCURACY)
        hdr(3) = dflat(irid)      !CLAT - LATITUDE (COARSE ACCURACY)
        hdr(4) = dfheight(irid)   !SELV - HEIGHT OF STATION
        !hdr(5) goes in tilt loop below.
        hdr(6) = iadate(1)        !YEAR - YEAR
        hdr(7) = iadate(2)        !MNTH - MONTH
        hdr(8) = iadate(3)        !DAYS - DAY
        hdr(9) = iadate(4)        !HOUR - HOUR 
        hdr(10)= 00               !MINU - MINUTE
        hdr(11)= 1                !QCRW - QUALITY MARK FOR WINDS ALONG RADIAL LINE
        !hdr(12) goes in azm loop below.
        bufrfilename=trim(idate)//'_fv3.t'//trim(hh)//'z.drw.bufr'
        bufrtilt: do itiltbufr=1,1!nelv
           hdr(5) = tilts(itiltbufr) 
           if(.not.bufrisopen) then !open a new message for each tilt 
              open(unit=10,file=trim(bufrfilename),status='unknown',action='write',form='unformatted')
              open(unit=11,file='l2rwbufr.table',status='old',action='read',form='formatted')
              call openbf(10,'OUT',11)
              bufrisopen=.true.
           end if
           intdate=iadate(1)*1000000 + iadate(2)*10000 + iadate(3)*100 + iadate(4) ! int(yyyymmddhh)
           call openmb(10,trim(subset),intdate)
           bufrazm: do iazmbufr=1,360
              allocate(obs(3,numgates))
              obs=-999.0_r_kind ! Initialize as missing values
              hdr(12)=float(iazmbufr)
              bufrgate: do igatebufr=1,numgates
                 obs(1,igatebufr) = igatebufr     !DISTANCE (FROM ANTENNA TO GATE CENTER) IN UNITS OF 250M
                 obs(2,igatebufr) = drwpol(itiltbufr,iazmbufr,igatebufr)  !DOPPLER MEAN RADIAL VELOC            
                 obs(3,igatebufr) = 1.0_r_kind                            !DOPPLER VELOCITY SPECTRAL WIDTH
              end do bufrgate
              ! encode radial velocity
              call ufbint(10,hdr,12,1,iret,trim(hdstr))
              call ufbint(10,obs, 3,numgates,iret,trim(obstr))
              call writsb(10)
              deallocate(obs)
           end do bufrazm
           call closmg(10) ! close bufr message
        end do bufrtilt
        call closbf(10) !close bufr file
        close(10)       ! close bufr file
        close(11)       ! close l2rwbufr.table
        bufrisopen=.false.


     end if ifKGRK
  end do loopOVERradars 


  print *, "end of program" 
end program drwsim


  subroutine tll2xy(rlon,rlat,x,y)!,outside)

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
