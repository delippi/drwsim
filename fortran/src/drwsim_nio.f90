program drwsim

  use kinds,     only: r_kind,r_double,i_kind,r_single
  use constants, only: zero,izero,half,one,ione,two,deg2rad,rearth,rad2deg,&
                       one_tenth,r1000,r60,r60inv,r100,r400,init_constants_derived,&
                       init_constants,tiny_r_kind
  use constants, only: flattening,semi_major_axis,grav_ratio,grav,wgtlim,&
                       grav_equator,eccentricity,somigliana
  use interp_util
  use netcdf
  use nemsio_module, only: nemsio_gfile,nemsio_open,nemsio_close,nemsio_getfilehead,&
                         nemsio_getheadvar,nemsio_realkind,nemsio_intkind,&
                         nemsio_readrecv,nemsio_init,nemsio_setheadvar,nemsio_writerecv,&
                         nemsio_readrec

  implicit none

  integer(i_kind)     :: ndata
  integer(i_kind)     :: nradars
  integer(i_kind)     :: mindat

  ! Declare local parameters
  real(r_kind),parameter :: four_thirds = 4.0_r_kind / 3.0_r_kind
  real(r_kind),parameter :: r8          = 8.0_r_kind
  real(r_kind),parameter :: r360        = 360.0_r_kind

!--General declarations
  integer(i_kind) :: iret
  integer(i_kind),dimension(4) :: iadate
  integer(i_kind),dimension(4) :: intdate
  character(4) :: yyyy
  character(2) :: mm
  character(2) :: dd
  character(2) :: hh
  character(2) :: mn
  character(10):: cdate
  
  real(r_kind) :: cosazm
  real(r_kind) :: sinazm
  real(r_kind) :: costilt
  real(r_kind) :: sintilt
  real(r_kind) :: ugesin
  real(r_kind) :: vgesin
  real(r_kind) :: wgesin
  real(r_kind) :: dbzgesin
  real(r_kind) :: zsges
  real(r_kind) :: dpres,zob
  real(r_kind) :: b,c,ha,epsh,h,aactual,a43,thistilt
  real(r_kind) :: thistiltr,selev0,celev0,thisrange,thishgt
  real(r_kind) :: celev,selev,gamma,thisazimuthr,thisazimuth,rlon0,rlat0,stahgt, &
                  clat0,slat0,dlat,dlon,thislon,thislat, &
                  rlonloc,rlatloc,rlonglob,rlatglob,rad_per_meter, &
                  radar_x,radar_y,radar_lon,radar_lat
  !real(r_kind),allocatable :: delz(:,:),height(:,:),delp(:,:)
  real(r_kind),allocatable :: drwpol(:,:,:) !tilt,azm,gate
  integer(i_kind) :: irid,itilt,iazm,igate,itime,iazm90,isig

  character(4) this_staid
  logical   :: diagprint
  logical   :: inside
  logical   :: bufrisopen
  logical   :: radar_location
  integer   :: diagverbose

  !---------DEFAULT SETTINGS---------!
  character(6)    :: datatype = 'NEMSIO'       ! Input data format for reading (NEMSIO or *NETCDF) *doesn't work
  integer(i_kind) :: maxobrange=250000_i_kind  ! Range (m) *within* which to use observations 
  integer(i_kind) :: minobrange=20000_i_kind   ! Range (m) *outside* of which
  real(r_kind)    :: mintilt=0.0_r_kind        ! Only use tilt(elevation) angles (deg) >= this number 
  real(r_kind)    :: maxtilt=5.5_r_kind        ! Do no use tilt(elevation) angles (deg) <= this number
  integer(i_kind) :: ithin=4_i_kind            ! Gates to skip for ob thinning (must be >=1)
  character(4)    :: staid='KOUN'              ! default station ID to use
  real(r_kind)    :: mindbz=-999_r_kind        ! minimum dbz value needed at a location to create a drw
  real(r_kind),allocatable :: tilts(:)
  integer(i_kind) :: vcpid=212_i_kind          ! default volume coverage pattern (VCP). 
  integer(i_kind) :: azimuths=360_i_kind       ! number of azimuths
  integer(i_kind) :: gatespc=250_i_kind        ! gate spacing (meters)
  integer(i_kind) :: numgates=400_i_kind       ! number of gates
  integer(i_kind) :: ntime=1_i_kind            ! number of times from nc file
  integer(i_kind) :: nelv=1_i_kind             ! number of elvation tilts
  logical         :: use_dbz=.true.            ! check for dbz at obs location?
  !----------------------------------------------!

  !---------NETCDF/NEMSIO VARS---------!
  real(r_kind),allocatable  ::     ges_u(:,:,:,:) ! (nlon,nlat,nsig,ntime)
  real(r_kind),allocatable  ::     ges_v(:,:,:,:) ! (nlon,nlat,nsig,ntime)
  real(r_kind),allocatable  ::     ges_w(:,:,:,:) ! (nlon,nlat,nsig,ntime)
  real(r_kind),allocatable  ::      delz(:,:,:,:) ! (nlon,nlat,nsig,ntime)
  real(r_kind),allocatable  ::      delp(:,:,:,:) ! (nlon,nlat,nsig,ntime)
  real(r_kind),allocatable  ::   ges_dbz(:,:,:,:) ! (nlon,nlat,nsig,ntime)
  !----------------------------------------------!

  ! ALL NEMSIO ADDITIONS HERE BEFORE SORTING 
  type(nemsio_gfile) :: gfile
  integer(i_kind) :: nlonsin,nlatsin,nlevsin,idvc
  integer(i_kind) :: idate(4)
  integer(i_kind) :: nfhour
  integer(i_kind) :: ilat
  integer(i_kind) :: ilon
  integer(i_kind) :: k
  integer(i_kind) :: azmspc
  real(r_single),allocatable,dimension(:,:,:  ) :: nems_vcoord
  real(r_single),allocatable,dimension(:      ) :: lats
  real(r_single),allocatable,dimension(:      ) :: lons
  real(r_single),allocatable,dimension(:      ) :: rlats
  real(r_single),allocatable,dimension(:      ) :: rlons
  real(r_kind),  allocatable,dimension(:,:    ) :: ges_ps
  real(r_kind),  allocatable,dimension(:      ) :: ak
  real(r_kind),  allocatable,dimension(:      ) :: bk
  real(r_kind),  allocatable,dimension(:,:,:  ) :: presi
  real(r_kind),  allocatable,dimension(:,:,:  ) :: hgti
  real(r_kind),  allocatable,dimension(:,:    ) :: ges_z
  real(r_kind),  allocatable,dimension(:,:    ) :: work
  real(r_kind),  allocatable,dimension(:,:,:,:) :: zges
  real(r_kind),  allocatable,dimension(:,:,:,:) :: pges
  real(r_kind),  allocatable,dimension(:      ) :: zges1d_sliced
  real(r_kind),              dimension(nlon*nlat) :: u1d
  real(r_kind),              dimension(nlon*nlat) :: v1d
  real(r_kind),              dimension(nlon*nlat) :: w1d
  real(r_kind),              dimension(nlon*nlat) :: sz1d
  real(r_kind),              dimension(nlon*nlat) :: sp1d
  real(r_kind),              dimension(nlon*nlat) :: dz1d
  real(r_kind),              dimension(nlon*nlat) :: dp1d
  real(r_kind) :: dlonm1,dlonp1,dlatm1,dlatp1,radar_xa,radar_ya,radar_xb,radar_yb
  logical :: geometric_height_switch
  logical :: test_delzdelp_vs_akbk 


  !---------GLOBAL RADAR CSV FILE VARS---------!
  integer(i_kind)                         :: ii
  integer(i_kind)                         :: numradars
  real(r_kind),dimension(:),allocatable   :: dflat
  real(r_kind),dimension(:),allocatable   :: dflon
  real(r_kind),dimension(:),allocatable   :: dfheight
  character(12),dimension(:),allocatable  :: dfid
  !----------------------------------------------!


  !---------BUFR VARS--------------------------!
  integer(i_kind) :: itiltbufr
  integer(i_kind) :: iazmbufr
  integer(i_kind) :: igatebufr
  integer(i_kind) :: iazmbufr90
  character(80)   :: bufrfilename
  character(80)   :: hdstr
  character(80)   :: obstr
  real(r_double)  :: hdr(12)
  real(r_double),allocatable :: obs(:,:)
  character(8) :: chdr
  character(8) :: subset
  equivalence (hdr(1),chdr)

  !---------L2RWBUFR CSV TABLE VARS------------!
  character(10) :: message_type
  character(8) :: cdummy

  !-------------TIMER VARS---------------------!
  integer(i_kind) :: time_array_0(8)
  integer(i_kind) :: time_array_1(8)
  integer(i_kind) :: hrs
  integer(i_kind) :: mins
  integer(i_kind) :: secs
  real(r_kind)    :: start_time
  real(r_kind)    :: finish_time
  real(r_kind)    :: total_time


  character(len=180)  :: datapath
  character(len=180)  :: nesteddata3d
  character(len=180)  :: nesteddata2d
  character(len=180)  :: nestedgrid
  character(len=180)  :: ak_bk
  character(len=180)  :: radarcsv
  character(len=180)  :: nesteddatadbz
  character(len=180)  :: filename

  namelist/drw/datatype,iadate,ntime,staid,ithin,mintilt,maxtilt,maxobrange,minobrange,&
               azimuths,use_dbz,mindbz,gatespc,diagprint,diagverbose,radarcsv,vcpid
  namelist/nc/datapath,nesteddata3d,nesteddata2d,nestedgrid,ak_bk,nesteddatadbz
  namelist/nio/datapath,filename
     

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
  open(11,file='namelist')
  read(11,drw)
  if(datatype == 'NEMSIO') read(11,nio)
  !----READ NAMELIST----

  !----CREATE VOLUME COVERAGE PATTERN (VCP)----start 
  !----SIMPLIFICATION NOTES:
  !    Most precipitation vcps have multiple scans at each layer.
  !    Further, there are SAILS1,2,&3 which jump back to 0.5 every so often.
  !    All simulated obs are assumed to be valid at the background timestamp.
  !    I cannot include this, because I do not have minutely output to do so.
  if(vcpid == 215) then
     nelv=15
     allocate(tilts(nelv))
     tilts(1:nelv)=(/ real(r_kind) :: 0.5,0.9,1.3,1.8,2.4,3.1,4.0,5.1,6.4,8.0,10.0,12.0,14.0,16.7/)
   else if(vcpid == 11 .or. vcpid == 211) then
     !SIMPLIFIED VCP 215 - the actual VCP has 2 scans at each level (28 total tilts).
     nelv=14_i_kind
     allocate(tilts(nelv))
     tilts(1:nelv)=(/ real(r_kind) :: 0.5,1.5,2.4,3.4,4.3,5.3,6.2,7.5,8.7,10.0,12.0,14.0,16.7,19.5/)
   else if(vcpid == 12 .or. vcpid == 212) then
     !SIMPLIFIED VCP 212 - the actual VCP has 2 scans at each level (28 total tilts).
     nelv=14_i_kind
     allocate(tilts(nelv))
     tilts(1:nelv)=(/ real(r_kind) :: 0.5,0.9,1.3,1.8,2.4,3.1,4.0,5.1,6.4,8.0,10.0,12.5,15.6,19.5/)
   else if(vcpid == 21 .or. vcpid == 121 .or. vcpid == 221) then
     !SIMPLIFIED VCP 221 - the actual VCP has 2 scans at each level (18 total tilts).
     nelv=9_i_kind
     allocate(tilts(nelv))
     tilts(1:nelv)=(/ real(r_kind) :: 0.5,1.5,2.4,3.4,4.3,6.0,9.9,14.6,19.5/)
   else if(vcpid == 999 ) then ! my vcp
     nelv=25_i_kind
     allocate(tilts(nelv))
     tilts(1:nelv)=(/ real(r_kind) :: 0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,6.0,7.0,8.0,9.0,10.0,&
                                          11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0/)
   else if(vcpid == 998 ) then ! my vcp for testing
     nelv=1_i_kind
     allocate(tilts(nelv))
     tilts(1:nelv)=(/ real(r_kind) :: 0.5/)
   end if
  !----CREATE VOLUME COVERAGE PATTERN (VCP)----end


  gatespc=gatespc*ithin
  azmspc=360/azimuths
  numgates=int(maxobrange/gatespc) !calculate num gates based on namelist settings.
  filename=trim(datapath) // trim(filename)

  if(diagprint .and. diagverbose >= 1) then
     write(*,*) "----NC FILE SPECS--------"
     write(*,*) "nlat",nlat
     write(*,*) "nlon",nlon
     write(*,*) "nsig",nsig
     1001 format (a16,f5.2,a3)
     write(*,1001) "grid spacing = ~",360.0*111.0/nlon,"-km"

     write(*,*) "----Namelist settings----"
     write(*,*) "radarcsv:   ",radarcsv
     write(*,*) "diagprint:  ",diagprint
     write(*,*) "diagverbose:",diagverbose
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
     write(*,*) "datapath:   ",datapath
     write(*,*) "nemsio file:  ",filename
  end if
  if(modulo(360,azimuths) /= 0) stop "ERROR: !!! AZIMUTHS IN NAMELIST IS NOT A FACTOR OF 360-DEG. !!!"

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
  write(*,*) "Done: Reading Global Radar List"


!------ OPEN NEMSIO FILE FOR READING -----------------------
  if(datatype=='NEMSIO') then
     write(*,*) 'Open and Read NEMSIO files'
     allocate(work(nlon,nlat))
     work    = zero

     call nemsio_init(iret=iret)
     call nemsio_open(gfile,filename,'READ',iret=iret)
     call nemsio_getheadvar(gfile,'idate',idate,iret)
     call nemsio_getheadvar(gfile,'nfhour',nfhour,iret)
     call nemsio_getfilehead(gfile,iret=iret,dimx=nlonsin,dimy=nlatsin,dimz=nlevsin,idvc=idvc)

     ! read sfc height
     write(*,*) 'Reading surface height'
     allocate(ges_z(nlon,nlat))
     ges_z=0.0_r_kind ! initialize
     call nemsio_readrecv(gfile,'hgt','sfc',1,sz1d,iret=iret)
     ges_z(:,:)=reshape(sz1d(:),(/size(work,1),size(work,2)/))
     ges_z=ges_z ! already in correct units?

     ! read sfc pressure
     write(*,*) 'Reading surface pressure'
     allocate(ges_ps(nlonsin,nlatsin))
     ges_ps=0.0_r_kind ! initialize
     call nemsio_readrecv(gfile,'pres','sfc',1,sp1d,iret=iret)
     ges_ps(:,:)=reshape(sp1d(:),(/size(work,1),size(work,2)/))
     ges_ps=ges_ps*0.01_r_kind ! convert ps to millibars.

     ! ak bk --> interface pressure/height calculation
     write(*,*) 'Reading ak bk'
     allocate(nems_vcoord(nsig+1,3,2),presi(nlon,nlat,nsig+1),hgti(nlon,nlat,nsig),ak(nsig+1),bk(nsig+1))
     allocate(pges(nlon,nlat,nsig,1))
     allocate(zges(nlon,nlat,nsig,1))
     pges = 0.0_r_kind ! initialize
     zges = 0.0_r_kind ! initialize
     call nemsio_getfilehead(gfile,iret=iret,vcoord=nems_vcoord)
     if ( idvc == 2 ) then      ! hybrid coordinate
        ak = 0.01_r_kind*nems_vcoord(1:nsig+1,1,1) ! convert to mb
        bk = nems_vcoord(1:nsig+1,2,1)
     endif

     ! This is how you compute the pressure and geometric height at interfaces (not with delz/delp)
     do isig=1,nsig+1
        presi(:,:,isig)=ak(isig)+bk(isig)*ges_ps ! pressure at interfaces
        hgti(:,:,isig)=ak(isig)+bk(isig)*ges_z   ! geometric height at interfaces
     enddo
     test_delzdelp_vs_akbk=.false.
     if(.not. test_delzdelp_vs_akbk) zges(:,:,:,1)=hgti(:,:,:) ! save as correct variable name 
     !###########   TEST DELZ/DELP vs AKBK #################################
     if(test_delzdelp_vs_akbk) then
        allocate( delp(nlon,nlat,nsig,ntime))
        allocate( delz(nlon,nlat,nsig,ntime))
        delp(:,:,:,:)=0.0_r_kind ! initialize
        delz(:,:,:,:)=0.0_r_kind ! initialize
        itime=1
        do isig=1,nsig
           ! dpres
           call nemsio_readrecv(gfile,'dpres','mid layer',isig,dp1d,iret=iret)
           delp(:,:,isig,itime)=reshape(dp1d(:),(/size(work,1),size(work,2)/))
           ! delz
           call nemsio_readrecv(gfile,'delz','mid layer',isig,dz1d,iret=iret)
           delz(:,:,isig,itime)=reshape(dz1d(:),(/size(work,1),size(work,2)/))
        enddo
        do isig=1,nsig
           if(isig==1) then
              zges(:,:,isig,1)=                       delz(:,:,isig,itime)
           else
              zges(:,:,isig,1)=zges(:,:,isig-1,itime)+delz(:,:,isig,itime)
           end if
        enddo
        do isig=nsig,1,-1
           if(isig==nsig) then
              pges(:,:,isig,1)=                       delp(:,:,isig,itime)
           else
              pges(:,:,isig,1)=pges(:,:,isig+1,itime)+delp(:,:,isig,itime)
           endif
        end do
        print *, "pges ",pges(2242,506,1,1) ! print at the location of KGRK
        print *, "zges ",zges(2242,506,1,1)
        print *, "presi",presi(2242,506,1)
        print *, "hgti",hgti(2242,506,1)
        stop 'End of checking delzdelp vs akbk'
     end if ! end test delzdelp vs akbk
     !###################################################################

     ! rlats rlons
     write(*,*) 'Reading lat lon'
     allocate(lats(nlon*nlat),lons(nlon*nlat),rlats(nlat),rlons(nlon))
     call nemsio_getfilehead(gfile,iret=iret,lon=lons,lat=lats)
     do ilon=1,nlon
        rlons(ilon)=lons(ilon)*deg2rad
     end do
     do ilat=1,nlat
        rlats(ilat)=lats(ilat*nlon-(nlon-1))*deg2rad 
     end do
     !lats indicies 1-3072 are all the same and are 89.9
     !lons indicies 1-3072 are increase from 0 to 359.88

     deallocate(ak,bk,lats,lons)


     ! read 3d fields (u,v,w,dbz)
     itime=1
     allocate(ges_u(nlon,nlat,nsig,ntime))
     allocate(ges_v(nlon,nlat,nsig,ntime))
     allocate(ges_w(nlon,nlat,nsig,ntime))
     allocate(ges_dbz(nlon,nlat,nsig,ntime))
        ges_u(:,:,:,:)=0.0_r_kind ! initialize
        ges_v(:,:,:,:)=0.0_r_kind ! initialize
        ges_w(:,:,:,:)=0.0_r_kind ! initialize
        ges_dbz(:,:,:,:)=0.0_r_kind ! initialize

     write(*,*) 'Reading u,v,w level by level'
     do isig=1,nsig
        write(*,*) 'level = ',isig

        ! u
        call nemsio_readrecv(gfile,'ugrd','mid layer',isig,u1d,iret=iret)
        ges_u(:,:,isig,itime)=reshape(u1d(:),(/size(work,1),size(work,2)/))

        ! v
        call nemsio_readrecv(gfile,'ugrd','mid layer',isig,v1d,iret=iret)
        ges_v(:,:,isig,itime)=reshape(v1d(:),(/size(work,1),size(work,2)/))

        ! w
        call nemsio_readrecv(gfile,'ugrd','mid layer',isig,w1d,iret=iret)
        ges_w(:,:,isig,itime)=reshape(w1d(:),(/size(work,1),size(work,2)/))

        ! dbz


     end do

     print *, idate,nfhour

  write(*,*) 'Done: Reading NEMSIO files'
  end if ! NEMSIO READ

  !set up information needed to interpolate model to observation
  !*********************************
  call gridmod_extract
  
  !*********************************

  if (diagprint .and. diagverbose >= 10) then
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
  nradars=izero

  !--Call Timer
  call date_and_time(values=time_array_0)
  start_time = time_array_0(5)*3600 + time_array_0(6)*60 + time_array_0(7) + time_array_0(8)*0.001
  write(*,*) 'STARTING RADAR WIND SIMULATION:'

  loopOVERtime: do itime=1,ntime
     geometric_height_switch=.true.
     if(itime > 1) call newdate(iadate,1,iadate) ! don't increment the first time.
     loopOVERradars: do irid=1,numradars 
        allocate(drwpol(nelv,360,numgates))
        drwpol=-999.0_r_kind !Initialize/Reset the drw polar field
        radar_location=.true. ! logical to only compute radar x,y once later in loop.
        this_staid=adjustl(trim(dfid(irid)))
        ifKGRK: if(this_staid==trim(adjustl(staid)) .or. trim(adjustl(staid))=='all') then
           stahgt=dfheight(irid)
           rlon0=dflon(irid)*deg2rad
           rlat0=dflat(irid)*deg2rad
           clat0=cos(rlat0)
           slat0=sin(rlat0)
           loopOVERtilts:    do itilt=1,nelv
              if(tilts(itilt) > maxtilt .or. tilts(itilt) < mintilt) then
                 !write(*,*) "oops... itilt:",tilts(itilt)," is out of bounds max/mintilt:",maxtilt,mintilt," cycling"
                 cycle
              end if
              bufrisopen=.false.   !Initialize bufr file as closed.
              thistilt=tilts(itilt)
              thistiltr=thistilt*deg2rad
              celev0=cos(thistiltr)
              selev0=sin(thistiltr)
              loopOVERazimuths: do iazm=azmspc,360,azmspc
                 !1000 format(a5,1x,i4,i2.2,i2.2,i2.2,&
                 !         3x,a6,1x,a4,&
                 !         !3x,a5,1x,f4.1,a1,f4.1,&
                 !         3x,a5,1x,i2,a2,i2,1x,a1,f4.1,a1,f4.1,a1,&
                 !         3x,a4,1x,i3,a4,&
                 !         3x,a4,1x,i3,&
                 !         3x,a6,1x,i2,a1,i5,a9)
                 !write(*,1000),"Date:",iadate(1),iadate(2),iadate(3),iadate(4),&
                 !              "Radar:",adjustl(trim(dfid(irid))),&
                 !              !"Tilt:",tilts(itilt),"/",tilts(nelv),&
                 !              "Tilt:",itilt,"of",nelv,"(",tilts(itilt),"/",tilts(nelv),")",&
                 !              "Azm:",iazm,"/360",&
                 !              "VCP:",vcpid,&
                 !              "ObRes:",azmspc,"x",gatespc,"[deg x m]"
                 thisazimuth=90.0_r_kind-float(iazm) ! 90-azm to be consistent with l2rwbufr
                 if(thisazimuth>=r360) thisazimuth=thisazimuth-r360
                 if(thisazimuth<zero) thisazimuth=thisazimuth+r360
                 thisazimuthr=thisazimuth*deg2rad
                 loopOVERgates: do igate=1,numgates     
                    inside=.false. ! is our ob location inside the bounds? preset to false, then check.
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

                       !--Determine the x,y (grid relative location) of the radar location.
                       !--Only need to do this once per radar.
                       radar_location=.false.
                       if(radar_location) then
                          radar_lat=dflat(irid) !lat/lons stored as deg.
                          radar_lon=dflon(irid)
                          if(radar_lon>=r360) radar_lon=radar_lon-r360 !fix if needed.
                          if(radar_lon<zero) radar_lon=radar_lon+r360
                          radar_x=radar_lon
                          radar_y=radar_lat
                          call grdcrd1(radar_y,rlats*rad2deg,nlat,-1) !lats are in descending order
                          call grdcrd1(radar_x,rlons*rad2deg,nlon, 1)
                          radar_xa=int(radar_x)
                          radar_xb=int(radar_x)+2
                          radar_ya=int(radar_y)
                          radar_yb=int(radar_y)+2
                          !call tintrp2a_single_level_sliced(ges_z(radar_xm1:radar_xp1,radar_ym1:radar_yp1),&
                          !                            zsges,radar_x-int(radar_x)+2,radar_y-int(radar_y)+2)
                          radar_location=.false. ! turn off get radar x/y until next radar is processed.
                       end if

                       thislat=rlatglob*rad2deg
                       thislon=rlonglob*rad2deg
                       if(thislon>=r360) thislon=thislon-r360
                       if(thislon<zero) thislon=thislon+r360
                       dlat=thislat
                       dlon=thislon
                       !--Find grid relative location of the ob. !!call tll2xy(thislon,thislat,dlon,dlat)
                       call grdcrd1(dlat,rlats*rad2deg,nlat,-1) !lats are in descending order
                       call grdcrd1(dlon,rlons*rad2deg,nlon, 1)

                       !--Interpolate surface height to grid relative ob location (dlon,dlat). 
                       dlonm1=int(dlon)
                       dlonp1=int(dlon)+2
                       dlatm1=int(dlat)
                       dlatp1=int(dlat)+2
                       zsges=0.0_r_kind
                       call tintrp2a_single_level_sliced(ges_z(dlonm1:dlonp1,dlatm1:dlatp1),&
                                                   zsges,dlon-int(dlon)+2,dlat-int(dlat)+2)
                       !--Remove terrain height from ob absolute height and reject if below ground.
                       if(zsges>=dpres) then
                         !write(*,*) 'zsges =',zsges,'is greater than dpres',dpres,'. Rejecting ob.'
                         cycle
                       end if
                       dpres=dpres-zsges

                       !--Compute the geometric height from delz. (NOT THE
                       !CORRECT METHOD... USE AKBK FOR THIS!!! CHECK THE
                       !BEGGINING OF THE CODE FOR CORRECT METHOD!
                       !--Only need to do this once per run of this entire code.
                       !if(geometric_height_switch) then
                       !   geometric_height_switch=.false. ! turn off after first use.
                       !   !allocate(zges(nlon,nlat,nsig,1))
                       !   zges = 0.0_r_kind ! initialize
                       !   do isig=1,nsig
                       !      if(isig==1) then
                       !         zges(:,:,isig,1)=                       delz(:,:,isig,itime)
                       !      else
                       !         zges(:,:,isig,1)=zges(:,:,isig-1,itime)+delz(:,:,isig,itime)
                       !      end if
                       !   enddo
                       !end if

                       zob=dpres

                       !--Use a reduced (3x3) grid to compute the zges1d profile at ob location.
                       if(allocated(zges1d_sliced)) deallocate(zges1d_sliced)
                       allocate(zges1d_sliced(nsig))
                       zges1d_sliced(nsig)=0.0_r_kind ! initialize/reset after each iteration
                       dlonm1=int(dlon)
                       dlonp1=int(dlon)+2
                       dlatm1=int(dlat)
                       dlatp1=int(dlat)+2
                       do k=1,nsig
                          call tintrp2a_single_level_sliced(zges(dlonm1:dlonp1,dlatm1:dlatp1,k,itime),&
                                                  zges1d_sliced(k),dlon-int(dlon)+2,dlat-int(dlat)+2)
                       end do

                       !--Get vertical coordinate of observation given height of column at ob location.
                       call grdcrd1(dpres,zges1d_sliced,nsig,1) ! get grid coordinates of dpres from zges

                       !--Interpolate guess dbz to observation location - cycle if below threshold.
                       call tintrp3(ges_dbz(:,:,:,itime),dbzgesin,dlon,dlat,dpres)
                       dbzCheck: if(dbzgesin >= mindbz .and. use_dbz) then
                          !--Interpolate guess wind to observation location                                  
                          call tintrp3(ges_u(:,:,:,itime),ugesin,dlon,dlat,dpres)
                          call tintrp3(ges_v(:,:,:,itime),vgesin,dlon,dlat,dpres)
                          call tintrp3(ges_w(:,:,:,itime),wgesin,dlon,dlat,dpres)
                          !--Convert guess u,v,w wind components to radial value
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
                       end if dbzCheck
                    end if ifinside
                 end do loopOVERgates
              end do loopOVERazimuths
           end do loopOVERtilts

           if(diagprint .and. diagverbose >= 1) write(*,*) "min/max drw: ",minval(drwpol),maxval(drwpol)


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
           cdate=trim(yyyy)//trim(mm)//trim(dd)//trim(hh)
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
           !/gpfs/hps/nco/ops/com/rap/para/rap.20180504 or
           !/NCEPPROD/hpssprod/runhistory/rh2018/201805/20180503
           bufrfilename=trim(cdate)//'_fv3.t'//trim(hh)//'z.drw.bufr'  !rap.t00z.nexrad.tm00.bufr_d
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
              bufrazm: do iazmbufr=360/azimuths,360,360/azimuths
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
        nradars=nradars+1
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
  write(*,*) 'numradars =',nradars
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
