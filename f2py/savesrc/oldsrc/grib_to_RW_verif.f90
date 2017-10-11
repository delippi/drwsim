program grib_to_RW_verif
!$$$   subprogram documentation block
!                .      .    .                                       .
!   subprogram: grid_dbz        read level2 raw QC'd radial wind files and
!                                  interpolate forecast to ob location for verification
!
!   
!   prgmmr: carley          org:               date: 2011-09-01
!
! abstract: Reads and processes level 2 horizontal radial wind (m/s) by 
!                radar site.  Data are on radar scan surafces. Also reads, but does
!                not process unfolded radial velocities and reflectivity.  Processing includes
!                finding the lat/lon and height of each observation.  
!
!
!                It is important to note that the program has adopted a substantial
!                    amount of code from gridmod.F90 for the purposes of interpolating
!                    model forecast quantities to grid points.
!
!
! program history log:
!   2012-02-24   carley - bug fixes to colmax calculation.
!                          Now first put all obs within the time window
!                          into one large radar structure and perform the
!                          calculation.
!   2012-05-18   carley  - adapt the column max radar gridding code for
!                          use in verification
!           
!   input argument list:
!
! Variable Definitions:
!
!  a43 - real - (4/3)*(earth radius)   
!  a,b,c,ha,epsh,h,aactual - real - used in computing radar observation height 
!  celev0,selev0 - real- cos and sin of elevation angle (raw)
!  celev,selev - real - corrected cos and sin of elevation angle
!  clat0 - real - cos of radar station latitude
!  cstaid - char - radar station ide
!  dbzerr - real - observation error (obtained from convinfo - dBZ)
!  dlat - real - grid relative latitude of observation (grid units)
!  dlon - real - grid relative longitude of observation (grid units)
!  gamma - real - used in finding observation latlon
!  lunrad - int - unit number for reading radar data from file
!  maxobs - int - max number of obs converted to no precip observations
!  num_missing - int - number of missing observations
!  numbadtime - int - number of elevations outside time window
!  num_badtilt - int - number of elevations outside specified interval
!  num_badrange - int - number of obs outside specified range distance
!  obdate - int - dim(5) - yyyy,mm,dd,hh,minmin of observation
!  outside - logical - if observations are outside the domain -> true
!  rlatglob - real - earth relative latitude of observation (radians)
!  rlatloc - real - latitude of observation on radar-relative projection
!  rlonglob - real - earth relative longitude of observation (radians)
!  rlonloc - real - longitude of observation on radar-relative projection
!  rlon0 - real - radar station longitude (radians)
!  rmins_an - real - analysis time from reference date (minutes)
!  rmins_ob - real -  observation time from reference date (minutes)
!  rstation_id - real - radar station id
!  slat0 - real - sin of radar station latitude
!  thisazimuthr - real - 90deg minues the actual azimuth and converted to radians
!  thiserr - real - observation error
!  thislat - real - latitude of observation
!  thislon - real - longitude of observation
!  thisrange - real - range of observation from radar
!  thishgt - real - observation height
!  this_stahgt - real - radar station height (meters about sea level)
!  this_staid - char - radar station id
!  thistilt - real - radar tilt angle (degrees)
!  thistiltr - real- radar tilt angle (radians)
!  timeb - real - obs time (analyis relative minutes)
!
!  
!
! Derived data types
!
!  radar - derived data type for containing volume scan information
!     nelv- int - number of elevation angles 
!     radid - char*4 - radar ID (e.g. KAMA)
!     vcpnum - int - volume coverage pattern number
!     year - int - UTC
!     day - int - UTC
!     month - int - UTC
!     hour - in - UTC
!     minute - int - UTC
!     second - int - UTC
!     radhgt - real - elevation of the radar above sea level in meters (I believe
!              this includes the height of the antenna as well)
!     radlat - real - latitude location of the radar
!     radlon - real - longitude location of the radar
!     fstgatdis - real - first gate distance (meters)
!     gatewidth - real - gate width (meters)
!     elev_angle - real - radar elevation angle (degrees)
!     num_beam - int - number of beams
!     num_gate - int - number of gates
!     nyq_vel - real - nyquist velocity 
!     azim - real - azimuth angles
!     field - real - radar data variable (reflectivity or velocity)
!
! Defined radar types:
!  strct_in_vel - radar - contains volume scan information related to 
!                         radial velocity
!  strct_in_dbz - radar - contains volume scan information related to
!                         radar reflectivity
!                            raw radial velocity    
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$ end documentation block

  use kinds, only: r_kind,r_double,i_kind,r_single
  use constants, only: zero,izero,half,one,ione,two,deg2rad,rearth,rad2deg, &
                       one_tenth,r1000,r60,r60inv,r100,r400,init_constants_derived, &
		       init_constants,tiny_r_kind
  use constants, only: flattening,semi_major_axis,grav_ratio,grav,wgtlim,&
                       grav_equator,eccentricity,somigliana                                            
  use interp_util
  
  implicit none
  
  integer(i_kind)     :: nread,ndata,nodata,mindat

! Declare local parameters
  real(r_kind),parameter :: four_thirds = 4.0_r_kind / 3.0_r_kind
  real(r_kind),parameter :: r8     = 8.0_r_kind
  real(r_kind),parameter:: r360=360.0_r_kind
  
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
  real(r_kind),dimension(nsig):: zges,hges
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
  real(r_kind) :: dbzerr,rmins_an,rmins_ob,tz_base,tz_top ,radar_twindow                                                    

  
  character(8) cstaid
  character(4) this_staid
  logical,save :: radar_location=.true.

  logical      :: outside,diagprint,lpasttime,lbadlev
    
  type(radar),allocatable :: strct_in_vel(:,:),strct_in_dbz(:,:),rad(:),strct_in_velBAD(:,:)
  type(ob),allocatable    :: verif(:)
  
  !---------SETTINGS---------!
  integer(i_kind) :: maxobrange=250000_i_kind	 ! Range (m) *within* which to use observations - obs *outside* this range are not used
  integer(i_kind) :: minobrange=20000_i_kind 	 ! Range (m) *outside* of which to use observatons - obs *inside* this range are not used
  real(r_kind)    :: mintilt=0.0_r_kind   	 ! Only use tilt(elevation) angles (deg) >= this number 
  real(r_kind)    :: maxtilt=5.5_r_kind         ! Do no use tilt(elevation) angles (deg) <= this number
  integer(i_kind) :: thin_freq=4_i_kind          ! Gates to skip for ob thinning (must be >=1)
  real(r_kind)    :: radartwindow=2.5            ! Time window within which to grid observations   (minutes)                                           
  !----------------------------------------------!

!----GRIB VARS--------------------!
    
 INTEGER(i_kind)   :: jf=nlon*nlat
 INTEGER(i_kind)   :: kpds5,kpds6,kpds7
 INTEGER(i_kind)   :: kpds(200),kgds(200),jpds(200),jgds(200)
 INTEGER(i_kind)   :: n,kf,units
 INTEGER(i_kind)   :: ier,iret
 REAL(r_single)    :: tmp_grib(nlon,nlat) !ABSOLUTELY MUST BE SINGLE PRECISION
 real(r_kind)      :: ges_u(nlat,nlon,nsig),ges_v(nlat,nlon,nsig),geop_hgtl(nlat,nlon,nsig) !GEOP height from grib is NOT AGL
 real(r_kind)      :: ges_z(nlat,nlon)

 LOGICAL*1        :: lb(nlon,nlat)
 CHARACTER(len=120)  :: gbFile,infile,outfile  !infile is the ob file





namelist/rw_verif/iadate,radartwindow,infile,gbFile,outfile,diagprint,mintilt, &
                   maxtilt,thin_freq,maxobrange,minobrange

!--------------------------------------------------------------------------------------!
!                            END OF ALL DECLARATIONS                                   !
!--------------------------------------------------------------------------------------!
  
  !--Set up the constants module used here
  call init_constants_derived
  call init_constants(.true.)    !initialize regional constants
  mindat=ione
  diagprint=.false.
  !Obtain input file name and output (binary) file name
  
  !---Read the namelist
  !  ./grid_radar.x < namelist
  !
  
  !----READ NAMELIST-----!
   read(*,rw_verif)
  !----READ NAMELIST-----!  
     
   write(6,*)'----Namelist settings----'
   write(6,*)'Date:',iadate
   write(6,*)'radartwindow:',radartwindow
   write(6,*)'infile:',trim(infile)
   write(6,*)'outfile:',trim(outfile)
   write(6,*)'maxobrange:',maxobrange
   write(6,*)'minobrange:',minobrange
   write(6,*)'mintilt:',mintilt
   write(6,*)'maxtilt:',maxtilt
   write(6,*)'thin_freq:',thin_freq
   write(6,*)'diagprint:',diagprint
   write(6,*)'gbFile:',trim(gbFile)
   write(6,*)'outfile:',trim(outfile)
   write(6,*)'----END Namelist settings----'
  
    
  ikx=izero
  
  if (minobrange >= maxobrange) then
     STOP 'MININMUM OB RANGE >= MAXIMUM OB RANGE FOR READING Vr - PROGRAM STOPPING FROM grib_to_RW_verif.f90'
  end if
  
  if (thin_freq < 1) then
     STOP 'thin_freq MUST BE >= 1!  CHECK NAMELIST AND RESET!'
  end if
  
    
  !OPEN AND READ GRIB FILE


  ! U-component of the wind ! 

! Wildcards for searching Grib file
jpds = -1
jgds = -1
! Search for this parameter
!  (211) on model levels (109)
jpds(5) = 33 
jpds(6) = 109

 ! initialize
 ges_u=0.0_r_kind
 ! Open and read grib file
 units=35
 CALL BAOPENR(units,TRIM(gbFile),ier)
 write(6,*) 'Reading file ',trim(gbFile),' with ier:',ier
 if (ier /= 0 ) then
    write(6,*) 'ERROR OPENING GRIB FILE! STOPPING!'
    STOP
 end if
 
 do n=1,nsig               
        jpds(7) = n   !corresponds to model vertical level in 'traditional coords - NOT NEMSIO        
        tmp_grib=0_r_single
        CALL GETGB(units,0,jf,0,jpds,jgds,kf,k,kpds,kgds,lb,tmp_grib,iret)
        if (iret /= 0 ) then
           write(6,*) 'ERROR READING GRIB FILE! STOPPING!'
           STOP
        end if        
!        write(6,*) 'Max/maxloc/min/minloc of U-WIND with close:',maxval(tmp_grib),maxloc(tmp_grib),minval(tmp_grib),minloc(tmp_grib)        
        !store and flip arrays so they jive with interpolation routines
        do j=1,nlat
          do i=1,nlon
            ges_u(j,i,n)=tmp_grib(i,j)
          end do
        end do               
 end do  


! V-component of the wind ! 

! Wildcards for searching Grib file
jpds = -1
jgds = -1
! Search for this parameter
!  (211) on model levels (109)
jpds(5) = 34 
jpds(6) = 109

 ! initialize
 ges_v=0.0_r_kind
 do n=1,nsig               
        jpds(7) = n   !corresponds to model vertical level in 'traditional coords - NOT NEMSIO        
        tmp_grib=0_r_single
        CALL GETGB(units,0,jf,0,jpds,jgds,kf,k,kpds,kgds,lb,tmp_grib,iret)
        if (iret /= 0 ) then
           write(6,*) 'ERROR READING GRIB FILE! STOPPING!'
           STOP
        end if        
!        write(6,*) 'Max/maxloc/min/minloc of V-WIND with close:',maxval(tmp_grib),maxloc(tmp_grib),minval(tmp_grib),minloc(tmp_grib)        
        !store and flip arrays so they jive with interpolation routines
        do j=1,nlat
          do i=1,nlon
            ges_v(j,i,n)=tmp_grib(i,j)
          end do
        end do               
 end do  




 
 
  ! READ IN THE TERRAIN HEIGHT ges_z

 ! Wildcards for searching Grib file
 jpds = -1
 jgds = -1
! Search for this parameter
 jpds(5) = 7
 jpds(6) = 1
 jpds(7) = 0
 ges_z=zero
 tmp_grib=0_r_single
 CALL GETGB(units,0,jf,0,jpds,jgds,kf,k,kpds,kgds,lb,tmp_grib,iret)
 if (iret /= 0 ) then
    write(6,*) 'ERROR READING GRIB FILE! STOPPING!'
    STOP
 end if 
! write(6,*) 'Max/maxloc/min/minloc of surface height field with close:',maxval(tmp_grib),maxloc(tmp_grib),minval(tmp_grib),minloc(tmp_grib)
!--store and flip arrays so they jive with interpolation routines - also subtract of the model terrain here as is done in load_prsges
 do j=1,nlat
   do i=1,nlon
     ges_z(j,i)=tmp_grib(i,j)
   end do
 end do   
 

 ! READ IN THE GEOP HEIGHT   AND SUBTRACT THE TERRAIN OFF THE GEOP HEIGHT TO GET AGL
 
!rec 6:813282:date 2007042118 HGT kpds5=7 kpds6=109 kpds7=1 levels=(0,1) grid=255 hybrid lev 1 7hr fcst:
!  HGT=Geopotential height [gpm]
 
 
! Wildcards for searching Grib file
jpds = -1
jgds = -1
! Search for this parameter
! reflectivity (211) on model levels (109)
jpds(5) = 7 
jpds(6) = 109

 ! initialize
 geop_hgtl=0.0_r_kind
 ! Open and read grib file
 CALL BAOPENR(units,TRIM(gbFile),ier)
 do n=1,nsig               
        jpds(7) = n   !corresponds to model vertical level in 'traditional coords - NOT NEMSIO        
        tmp_grib=0_r_single
        CALL GETGB(units,0,jf,0,jpds,jgds,kf,k,kpds,kgds,lb,tmp_grib,iret)
        if (iret /= 0 ) then
           write(6,*) 'ERROR READING GRIB FILE! STOPPING!'
           STOP
        end if        
!        write(6,*) 'Max/maxloc/min/minloc of geop_hgtl field with close:',maxval(tmp_grib),maxloc(tmp_grib),minval(tmp_grib),minloc(tmp_grib)
                
        !store and flip arrays so they jive with interpolation routines - also subtract of the model terrain here as is done in load_prsges
        do j=1,nlat
          do i=1,nlon
            geop_hgtl(j,i,n)=tmp_grib(i,j)-ges_z(j,i)
          end do
        end do 
!                write(6,*) 'Max/maxloc/min/minloc of geop_hgtl AGL field with close:',maxval(geop_hgtl(:,:,n)),maxloc(geop_hgtl(:,:,n)),minval(geop_hgtl(:,:,n)),minloc(geop_hgtl(:,:,n))              
 end do  
  
    
  !Read in lat and lon arrays

  allocate(glon(nlon,nlat)) ! unlike rest of gsi these arrays are lon,lat
  allocate(glat(nlon,nlat))

!-----LONGITUDES

! Wildcards for searching Grib file
 jpds = -1
 jgds = -1
! Search for this parameter
 jpds(5) = 177 
 jpds(6) = 1
 jpds(7) = 0
 glon=zero
 tmp_grib=0_r_single
 CALL GETGB(units,0,jf,0,jpds,jgds,kf,k,kpds,kgds,lb,tmp_grib,iret)
 if (iret /= 0 ) then
    write(6,*) 'ERROR READING GRIB FILE! STOPPING!'
    STOP
 end if 
 glon=tmp_grib
! write(6,*) 'Max/maxloc/min/minloc of glon field with close:',maxval(glon),maxloc(glon),minval(glon),minloc(glon)
      
   !!!!!!!!! DO NOT NEED TO FLIP THE INDICES OF THE LON AND LAT ARRAYS!!!!!!!!
               
!-----LATITUDES
 
! Wildcards for searching Grib file
jpds = -1
jgds = -1
! Search for this parameter
! reflectivity (211) on model levels (109)
jpds(5) = 176 
jpds(6) = 1
jpds(7) = 0

 glat=zero
 tmp_grib=0_r_single
 CALL GETGB(units,0,jf,0,jpds,jgds,kf,k,kpds,kgds,lb,tmp_grib,iret)
 if (iret /= 0 ) then
    write(6,*) 'ERROR READING GRIB FILE! STOPPING!'
    STOP
 end if
 glat=tmp_grib        
! write(6,*) 'Max/maxloc/min/minloc of glat field with close:',maxval(glat),maxloc(glat),minval(glat),minloc(glat)        

 
 ! - CLOSE THE GRIB FILE
 CALL BACLOSE(units,ier)
 if (ier /= 0 ) then
    write(6,*) 'ERROR CLOSING GRIB FILE! STOPPING!'
    STOP
 end if
 
      
!!!!!!!!!! DO NOT NEED TO FLIP THE INDICES OF THE LON AND LAT ARRAYS!!!!!!!!  
  
     
  !set up information needed to interpolate model to observation
  !*********************************
  call gridmod_extract
  !*********************************

   if (diagprint) then
     write(6,*)'----Model grid info diagnostic print----'
     write(6,*)'nlon,nlat,rlambda0,pihalf,sign_pole:',nlon,nlat,rlambda0,pihalf,sign_pole
     write(6,*)'atilde_x,atilde_y,btilde_x,btilde_y:',atilde_x,atilde_y,btilde_x,btilde_y
     write(6,*)'rlon_min_dd,rlon_max_dd,rlat_min_dd,rlat_max_dd,nxtilde,nytilde:',rlon_min_dd,rlon_max_dd,rlat_min_dd,rlat_max_dd,nxtilde,nytilde
     write(6,*)'Max/min i0_tilde,j0_tilde:',maxval(i0_tilde),minval(i0_tilde),maxval(j0_tilde),minval(j0_tilde)
     write(6,*)'Max/min ip_tilde,jp_tilde:',maxval(ip_tilde),minval(ip_tilde),maxval(jp_tilde),minval(jp_tilde)
     write(6,*)'Max/min xtilde0,ytilde0:',maxval(xtilde0),minval(xtilde0),maxval(ytilde0),minval(ytilde0)  
     write(6,*)'----End model grid info diagnostic print----'
  end if


       
  !-next three values are dummy values for now
  nchanl=izero
  ilon=2_i_kind
  ilat=3_i_kind
  
  maxobs=2000000_i_kind    !value taken from read_radar.f90 

  maxgate=10000_i_kind     !maximum number of gates.  Only needed in this application here
                           !  After looking through obs data it seems there are never more
			   !  than 920 gates. 
  ndata=izero
  nread=izero
   

   !--Allocate verification object

  allocate(verif(maxobs))
  verif%omf=zero
  verif%omf2=zero  
  verif%ob=zero     
  verif%ges=zero 
  
  
  
  
  
   
  !-Obtain analysis time in minutes since reference date

  call w3fs21(iadate,mins_an)  !mins_an -integer number of mins snce 01/01/1978
  rmins_an=mins_an             !convert to real number
  
       
  lunrad=31_i_kind
  open(lunrad,file=trim(infile),status='old',action='read', &
       iostat=ierror,form='formatted')
  
 fileopen: if (ierror == izero) then           !Check to make sure file is open - will also fail if file does not exist. Closing endif at end of subroutine.
  
  read(lunrad,'(2i8)') nelv,nvol                    !read number of elevations and number of volumes
  allocate(strct_in_vel(1,nelv))  !not concerned with the velocities here FOR NOW                          CARLEY CARLEY CARLEY CARLEY CARLEY CARLEY CARLEY CARLEY CARLEY CARLEY CARLEY 
  allocate(strct_in_velBAD(1,nelv))
  allocate(strct_in_dbz(1,nelv))
  rad_nelv=nint(radartwindow)*nelv
  rad_nelv=29  !carley debug
  allocate(rad(rad_nelv))

  
  
  !Initialize rad to missing values
   do cm=1,rad_nelv
    rad(cm)%field=999.0_r_kind
   end do
!-----Code to read data is based heavily upon code provided by Kang Nai of OU----!

  cm=0 !initialize counter for column max obs
  volumes: do v=ione,nvol
     
         
     read(lunrad,'(i8)') nelv
  
  elevs:  do k=ione,nelv
          
     !  Processed radial velocity (m/s)
     
        read(lunrad,'(a4)') strct_in_vel(1,k)%radid
        read(lunrad,'(i8)') strct_in_vel(1,k)%vcpnum
        read(lunrad,'(6i8)') strct_in_vel(1,k)%year              &
                         ,strct_in_vel(1,k)%month                &
                         ,strct_in_vel(1,k)%day                  &
                         ,strct_in_vel(1,k)%hour                 &
                         ,strct_in_vel(1,k)%minute               &
                         ,strct_in_vel(1,k)%second
        read(lunrad,'(2f10.3,f10.1)') strct_in_vel(1,k)%radlat   &
                                  ,strct_in_vel(1,k)%radlon      &
                                  ,strct_in_vel(1,k)%radhgt
        read(lunrad,'(2f8.1)') strct_in_vel(1,k)%fstgatdis       &
                           ,strct_in_vel(1,k)%gateWidth
        read(lunrad,'(f8.3)') strct_in_vel(1,k)%elev_angle
        read(lunrad,'(2i8)') strct_in_vel(1,k)%num_beam          &
                         ,strct_in_vel(1,k)%num_gate
        na=strct_in_vel(1,k)%num_beam
        nb=strct_in_vel(1,k)%num_gate
     
        !******allocate arrays within radar data type**********!
!           allocate(strct_in_vel(1,k)%azim(na))
!	The number of gates changes with elevation too, so let's just
!	  make one big, uniform size the convers all, and initialize
!	  it to the missing data value		   
!	   allocate(strct_in_vel(1,k)%field(maxgate,na))
	   strct_in_vel(1,k)%field=999.0_r_kind
        !******************************************************!
          
        read(lunrad,'(f8.3)') strct_in_vel(1,k)%nyq_vel
        read(lunrad,'(15f6.1)') (strct_in_vel(1,k)%azim(j),j=ione,na)
        read(lunrad,'(20f6.1)') ((strct_in_vel(1,k)%field(i,j),i=ione,nb),j=ione,na)
          
        !--Processed radar reflectivity factor (dBZ)
          
        read(lunrad,'(a4)') strct_in_dbz(1,k)%radid
        read(lunrad,'(i8)') strct_in_dbz(1,k)%vcpnum
        read(lunrad,'(6i8)') strct_in_dbz(1,k)%year              &
                         ,strct_in_dbz(1,k)%month                &
                         ,strct_in_dbz(1,k)%day                  &
                         ,strct_in_dbz(1,k)%hour                 &
                         ,strct_in_dbz(1,k)%minute               &
                         ,strct_in_dbz(1,k)%second
        read(lunrad,'(2f10.3,f10.1)') strct_in_dbz(1,k)%radlat   &
                                  ,strct_in_dbz(1,k)%radlon      &
                                  ,strct_in_dbz(1,k)%radhgt
        read(lunrad,'(2f8.1)') strct_in_dbz(1,k)%fstgatdis       &
                           ,strct_in_dbz(1,k)%gateWidth
        read(lunrad,'(f8.3)') strct_in_dbz(1,k)%elev_angle
        read(lunrad,'(2i8)') strct_in_dbz(1,k)%num_beam          &
                         ,strct_in_dbz(1,k)%num_gate
        na=strct_in_dbz(1,k)%num_beam
        nb=strct_in_dbz(1,k)%num_gate

	
	
	!******allocate arrays within radar data type**********!
!        allocate(strct_in_dbz(1,k)%azim(na))       
!	The number of gates changes with elevation too, so let's just
!	  make one big, uniform size the convers all, and initialize
!	  it to the missing data value	
!        allocate(strct_in_dbz(1,k)%field(maxgate,na))
	strct_in_dbz(1,k)%field=999.0_r_kind    
        !******************************************************!
	
	  
        read(lunrad,'(f8.3)') strct_in_dbz(1,k)%nyq_vel
        read(lunrad,'(15f6.1)') (strct_in_dbz(1,k)%azim(j),j=ione,na)
        read(lunrad,'(20f6.1)') ((strct_in_dbz(1,k)%field(i,j),i=ione,nb),j=ione,na)
               
        !--Unprocessed, raw radial velocity (m/s)
          
        read(lunrad,'(a4)') strct_in_velBAD(1,k)%radid
        read(lunrad,'(i8)') strct_in_velBAD(1,k)%vcpnum
        read(lunrad,'(6i8)') strct_in_velBAD(1,k)%year              &
                         ,strct_in_velBAD(1,k)%month                &
                         ,strct_in_velBAD(1,k)%day                  &
                         ,strct_in_velBAD(1,k)%hour                 &
                         ,strct_in_velBAD(1,k)%minute               &
                         ,strct_in_velBAD(1,k)%second
        read(lunrad,'(2f10.3,f10.1)') strct_in_velBAD(1,k)%radlat   &
                                  ,strct_in_velBAD(1,k)%radlon      &
                                  ,strct_in_velBAD(1,k)%radhgt
        read(lunrad,'(2f8.1)') strct_in_velBAD(1,k)%fstgatdis       &
                           ,strct_in_velBAD(1,k)%gateWidth
        read(lunrad,'(f8.3)') strct_in_velBAD(1,k)%elev_angle
        read(lunrad,'(2i8)') strct_in_velBAD(1,k)%num_beam          &
                         ,strct_in_velBAD(1,k)%num_gate
        na=strct_in_velBAD(1,k)%num_beam
        nb=strct_in_velBAD(1,k)%num_gate
     
        !******allocate arrays within radar data type**********!
!        allocate(strct_in_vel(1,k)%azim(na))
!	The number of gates changes with elevation too, so let's just
!	  make one big, uniform size the convers all, and initialize
!	  it to the missing data value	
!        allocate(strct_in_vel(1,k)%field(maxgate,na))
	strct_in_velBAD(1,k)%field=999.0_r_kind
        !******************************************************!
     
        read(lunrad,'(f8.3)') strct_in_velBAD(1,k)%nyq_vel
        read(lunrad,'(15f6.1)') (strct_in_velBAD(1,k)%azim(j),j=ione,na)
        read(lunrad,'(20f6.1)') ((strct_in_velBAD(1,k)%field(i,j),i=ione,nb),j=ione,na)     

    end do elevs

    tiltcheck:  do k=ione,nelv

      ! Must loop through once first and make necessary points missing that exceed the time or
      !  or specified tilt thresholds.  Additionally, this step also 'bins up' the obs within
      !  the time window as well as makes sure the correct tilts are selected and loaded.

     !--Check if observation fits within specified time window--!
      !-Find reference time of observation
     
        obdate(1)=strct_in_vel(1,k)%year
        obdate(2)=strct_in_vel(1,k)%month  
        obdate(3)=strct_in_vel(1,k)%day	 
        obdate(4)=strct_in_vel(1,k)%hour   
        obdate(5)=strct_in_vel(1,k)%minute 
        call w3fs21(obdate,mins_ob)                             !mins_ob -integer number of mins snce 01/01/1978
	rmins_ob=mins_ob                                        !convert to real number
	rmins_ob=rmins_ob+(strct_in_vel(1,k)%second*r60inv)     !convert seconds to minutes and add to ob time
 
      !-Comparison is done in units of minutes
      
        timeb = rmins_ob-rmins_an
        lpasttime=.false.
        lbadlev=.false. 
        if(abs(timeb) > abs(radartwindow)) then	  
           !And set to missing value so avoid affecting computations and set flag to NOT load                  
           strct_in_vel(1,k)%field = 999.0_r_kind
           lbadlev=.true.
           numbadtime=numbadtime+1
           if (timeb > zero ) lpasttime=.true.  !exit condition to stop reading the obsfile - speeds things up a lot!                                     
	end if
        
        thistilt=strct_in_vel(1,k)%elev_angle                  
        if (thistilt > maxtilt .or. thistilt < mintilt) then
           if (.not. lbadlev) num_badtilt=num_badtilt+1  !If it was not a bad level from the time check, it is now           
           !And set to missing value so avoid affecting computations and set flag to NOT load                   
           strct_in_vel(1,k)%field = 999.0_r_kind
           lbadlev=.true.      
        end if     
        if (lpasttime) exit tiltcheck
       
        if (.not. lbadlev) then   !If the data are not missing we should load it up for interpolation CARLEY CARLEY CARLEY CARLEY CARLEY CARLEY CARLEY CARLEY 
           cm=cm+1
           write(6,*) 'Loading time:',obdate,strct_in_vel(1,k)%second,'at tilt number',k,'and cm =',cm         
           rad(cm)=strct_in_vel(1,k)      
       end if
                        
    end do  tiltcheck
    if (lpasttime) exit volumes
  
  end do volumes !volume 
        
 else  !fileopen 
  write(6,*) 'VERIF_RW: ERROR OPENING RADAR RADIAL WIND FILE: ',trim(infile),' IOSTAT ERROR: ',ierror
 end if fileopen    
 
 close(lunrad)

  ! - deallocate fields used to read in radar data
  deallocate(strct_in_vel)
  deallocate(strct_in_velBAD)
  deallocate(strct_in_dbz)

  
 write(6,*) 'Loaded',cm,'elevations for Vr obs space verification'      
   
       
 !------Reading radar data finished for for 1 volume---------------!       
       
       
       
    tilts:  do k=ione,cm  !cm was counted up in the pervious step and corresponds to the number tilts loaded into rad

     
     !*************************IMPORTANT***************************!
     !                                                             !
     !    All data = 999.0 correspond to missing or bad data       !       
     !                                                             !
     !*************************************************************!
  
       !-Time window check has already been performed-!
        
        thistilt=rad(k)%elev_angle
        if (thistilt <= maxtilt .and. thistilt >= mintilt) then 
          
          gates: do i=ione,rad(k)%num_gate,thin_freq    !if thin_freq is 4 and we are using 250m gates then this lets in 1ob per 1 km - similar
	                                                           !  res to dBZ since the azimuthal res is the same

           	  
   	      thisrange=rad(k)%fstgatdis + float(i-ione)*rad(k)%gateWidth
       
             !-Check to make sure observations are within specified range 

              if (thisrange <= maxobrange .and. thisrange >= minobrange) then	    
	  	   
	      azms: do j=ione,rad(k)%num_beam
	   
	           !-Check to see if this is a missing observation 
		    
		    nread=nread+ione
		 
                    if ( rad(k)%field(i,j) >= 999.0_r_kind ) then
	               num_missing=num_missing+ione
	  	       cycle azms                        !No reason to process the ob if it is missing 
      	            end if
		    

                   			
		   !--Find observation height using method from read_l2bufr_mod.f90										       
	         
		    this_stahgt=rad(k)%radhgt
                    aactual=rearth+this_stahgt                    
                    a43=four_thirds*aactual
                    thistiltr=thistilt*deg2rad
                    selev0=sin(thistiltr)
                    celev0=cos(thistiltr)   		   
		    b=thisrange*(thisrange+two*aactual*selev0)
                    c=sqrt(aactual*aactual+b)
                    ha=b/(aactual+c)
                    epsh=(thisrange*thisrange-ha*ha)/(r8*aactual)
                    h=ha-epsh
	            thishgt=this_stahgt+h 


                    dpres=thishgt   !store absolute ob height (m) in dpres
                   

		   !--Find observation location using method from read_l2bufr_mod.f90
		 
		   !-Get corrected tilt angle
	            celev=celev0
	            selev=selev0
           	    celev=a43*celev0/(a43+h)
		    selev=(thisrange*thisrange+h*h+two*a43*h)/(two*thisrange*(a43+h))
	          
		    gamma=half*thisrange*(celev0+celev)
	         
                   !-Get earth lat lon of observation
	         	 
                    rlon0=deg2rad*rad(k)%radlon
	       	    clat0=cos(deg2rad*rad(k)%radlat)
		    slat0=sin(deg2rad*rad(k)%radlat)		  
		    thisazimuthr=(90.0_r_kind-rad(k)%azim(j))*deg2rad   !Storing as 90-azm to
		                                                                   ! be consistent with 
										   ! read_l2bufr_mod.f90
		    rad_per_meter=one/rearth
		    rlonloc=rad_per_meter*gamma*cos(thisazimuthr)
                    rlatloc=rad_per_meter*gamma*sin(thisazimuthr)
		                  
		    call invtllv(rlonloc,rlatloc,rlon0,clat0,slat0,rlonglob,rlatglob)
                 
		    if (radar_location) then
		        radar_lat=rad(k)%radlat
			radar_lon=rad(k)%radlon		   
		        if(radar_lon>=r360) radar_lon=radar_lon-r360
                        if(radar_lon<zero ) radar_lon=radar_lon+r360
			!-Convert to radians
			radar_lon=radar_lon*deg2rad
			radar_lat=radar_lat*deg2rad
			!-Find grid relative location of the radar
		        call tll2xy(radar_lon,radar_lat,radar_x,radar_y,outside)			
		        write(6,*) trim(rad(k)%radid),' Radar x,y location on the model grid is:', radar_x,radar_y 		    
		        radar_location=.false.
		    end if
		 
		 
		    thislat=rlatglob*rad2deg
                    thislon=rlonglob*rad2deg 
  
                   !-Check format of longitude and correct if necessary
                 
		    if(thislon>=r360) thislon=thislon-r360
                    if(thislon<zero ) thislon=thislon+r360
                 
		   !-Convert back to radians                 
         		       		       
		    thislat = thislat*deg2rad
                    thislon = thislon*deg2rad
                 		 
		    !find grid relative lat lon locations of earth lat lon
                 
		    call tll2xy(thislon,thislat,dlon,dlat,outside)
                    if (outside) cycle azms             !If observation is outside the domain
		                                        ! then cycle, but don't increase range right away.
							! Domain could be rectangular, so ob may be out of
						        ! range at one end, but not the other.
                                                        

!**********************************DO WIND ROTATION IF RADIAL WIND***********************************************!                                                           
                       cosazm_earth=cos(thisazimuthr)
                       sinazm_earth=sin(thisazimuthr)
                       call rotate_wind_ll2xy(cosazm_earth,sinazm_earth,cosazm,sinazm,thislon,dlon,dlat)
                       azm=atan2(sinazm,cosazm)
                       if(abs(azm)>r400) then
                          ibadazm=ibadazm+1
                          cycle azms
                       end if                                     
!**********************************DO WIND ROTATION IF RADIAL WIND***********************************************!                                                         
                      
                      
                    call tintrp2a_single_level(ges_z,zsges,dlat,dlon)
                    !Subtract off the terrain height
                    
                    dpres=dpres-zsges      !  remove terrain height from ob absolute height
                    if (dpres<zero) cycle azms !  temporary fix to prevent out of bounds array reference in zges,prsltmp
                    
                    call tintrp2a(geop_hgtl,hges,dlat,dlon,nsig)  !geop_hgtl we will read from grib file (note, need to subtract
                                                                  ! off the terrain field
                      
                    !    Convert geopotential height at layer midpoints to geometric height using
                    !    equations (17, 20, 23) in MJ Mahoney's note "A discussion of various
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
                         termr = semi_major_axis /(one + flattening + grav_ratio -  &
                              two*flattening*sin2)
                         termrg = (termg/grav)*termr
                         do n=1,nsig
                            zges(n) = (termr*hges(n)) / (termrg-hges(n))  ! eq (23)
                         end do

                    !    Convert observation height (in dpres) from meters to grid relative
                    !    units.  Save the observation height in zob for later use.
                         zob = dpres
                         call grdcrd(dpres,1,zges,nsig,1)
                    
                    
                    !IF WIND THEN DO BELOW   
                    !    Interpolate guess u and v to observation location                                  
                         call tintrp3(ges_u,ugesin,dlat,dlon,dpres)
                         call tintrp3(ges_v,vgesin,dlat,dlon,dpres)
                    !    Convert guess u,v wind components to radial value consident with obs    
                         cosazm_rw  = cos(azm)  ! cos(azimuth angle)                       
                         sinazm_rw  = sin(azm)  ! sin(azimuth angle)                       
                         costilt_rw = cos(thistiltr)  ! cos(tilt angle)
                   !-------------WIND FORWARD MODEL-----------------------------------------! 	
                         fcstgesin = (ugesin*cosazm_rw+vgesin*sinazm_rw)*costilt_rw
                              
                              
                         if(dpres < zero .or. dpres > nsig) then
                         ! Ob is below ground or above the domain, cycle
                            cycle azms                    
                         end if       
                              
  
                    
                    !print arbitray point, compare with some grads output to see if the interpolation is right
                    if (diagprint .and. ndata==2000) then
                      print*,dlon,dlat,dpres
                      print*,thislat*rad2deg,thislon*rad2deg,thishgt
                      print*,rad(k)%field(i,j)
                      print*,thisrange,thistilt,rad(k)%azim(j)
                      print*,zsges
                      print*,fcstgesin
                    end if
                    
                    		    
		    this_staid=rad(k)%radid      !Via equivalence in declaration, value is propagated
		                                       !  to rstation_id used below.

		
		    ndata  = min(ndata+ione,maxobs)     
		    nodata = min(nodata+ione,maxobs)  !number of obs not used (no meaning here)
		
                                            
	            observation=rad(k)%field(i,j)		     
		    ddiff=observation-fcstgesin
                   
		    verif(ndata)%omf  = ddiff                 ! (O-F)
		    verif(ndata)%omf2 = ddiff*ddiff           ! (O-F)^2
		    verif(ndata)%ob   = observation           ! radar observation
		    verif(ndata)%ges  = fcstgesin             ! Forecast
		
		     
                  end do azms  !j
              else
                 num_badrange=num_badrange+ione      !If outside acceptable range, increment
	      end if   !Range check	
		
	   end do gates    !i
     
        else
           num_badtilt=num_badtilt+ione           !If outside acceptable tilts, increment
        end if         !Tilt check

     end do tilts !elevations

      

!---all looping done now print diagnostic output
if (diagprint) then
  write(6,*)'VERIF_RW: Reached eof on radar wind file'
  write(6,*)'VERIF_RW: # volumes in input file                    =',nvol
  write(6,*)'VERIF_RW: # elevations per volume                    =',nelv
  write(6,*)'VERIF_RW: # bad azimuths                             =',ibadazm
  write(6,*)'VERIF_RW: # of missing data                          =',num_missing
  write(6,*)'VERIF_RW: # azimuths outside specif. range           =',num_badrange
  write(6,*)'VERIF_RW: # outside specif. tilts                    =',num_badtilt
  write(6,*)'VERIF_RW: # outside time range                       =',numbadtime
  write(6,*)'VERIF_RW: # obs used                                 =',ndata
  write(6,*)'VERIF_RW: # obs read                                 =',nread 
end if
!---Load the obs into new ob data type and call interpolation function--!
  deallocate(rad)

 
 !Finally, caluclate the verification stats
 sumbias=zero
 sumrms=zero
 rms=zero
 bias=zero
 
  ! verification step here.  Write out text file RMS,BIAS,ndata
  do i=1,ndata
      sumrms=sumrms+verif(i)%omf2
      sumbias=sumbias+verif(i)%omf     
  end do
  
  rms=sqrt( (one/ndata)*sumrms)
  bias=(one/ndata)*sumbias
  
  write(6,*)'RMS,BIAS,NDATA: ',rms,bias,ndata 
  
  
  deallocate(verif) 

  OPEN (UNIT=10, FILE=trim(outfile), STATUS='unknown', ACTION='WRITE', IOSTAT=ierror)
  write(10,*) rms,bias,ndata
  close(10)



!deallocate interpolation fields
 call destroy_interp
 
end program grib_to_RW_verif


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

    real(r_kind),intent(in   ) :: rlon  ! earth longitude (radians)
    real(r_kind),intent(in   ) :: rlat  ! earth latitude  (radians)

! !OUTPUT PARAMETERS:

    real(r_kind),intent(  out) :: x  ! x-grid coordinate (grid units)
    real(r_kind),intent(  out) :: y  ! y-grid coordinate (grid units)
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

    outside=x < rlon_min_dd .or. x > rlon_max_dd .or. &
            y < rlat_min_dd .or. y > rlat_max_dd

 end subroutine tll2xy


! ----- old tll2xy routine: may have some interpolation issues

!-------------------------------------------------------------------------
!    NOAA/NCEP, National Centers for Environmental Prediction GSI        !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  tll2xy --- convert earth lon-lat to x-y grid coordinates
!
! !INTERFACE:
!
  subroutine tll2xy_old(rlon,rlat,x,y,outside)

! !USES:

    use constants, only: one
    use interp_util
    implicit none

    real(r_kind),intent(in   ) :: rlon  ! earth longitude (radians)
    real(r_kind),intent(in   ) :: rlat  ! earth latitude  (radians)

! !OUTPUT PARAMETERS:

    real(r_kind),intent(  out) :: x  ! x-grid coordinate (grid units)
    real(r_kind),intent(  out) :: y  ! y-grid coordinate (grid units)
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

    outside=x < rlon_min_dd .or. x > rlon_max_dd .or. &
            y < rlat_min_dd .or. y > rlat_max_dd

 end subroutine tll2xy_old


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
