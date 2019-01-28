program drwsim

  use kinds,     only: r_kind,r_double,i_kind,r_single
  use constants, only: zero,izero,half,one,ione,two,deg2rad,rearth,rad2deg,&
                       one_tenth,r1000,r60,r60inv,r100,r400,init_constants_derived,&
                       init_constants,tiny_r_kind
  use constants, only: flattening,semi_major_axis,grav_ratio,grav,wgtlim,&
                       grav_equator,eccentricity,somigliana
  use constants, only: one,grav,half,huge_single
  use constants, only: eps,rd,t0c,fv
  use constants, only: cpf_a0,cpf_a1,cpf_a2,cpf_b0,cpf_b1,cpf_c0,cpf_c1,cpf_d,cpf_e
  use constants, only: psv_a,psv_b,psv_c,psv_d
  use constants, only: ef_alpha,ef_beta,ef_gamma
  use constants, only: rd_over_cp

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
  real(r_kind),parameter :: PI          = 3.141592653589793238462
  real(r_kind),parameter :: four_thirds = 4.0_r_kind / 3.0_r_kind
  real(r_kind),parameter :: r8          = 8.0_r_kind
  real(r_kind),parameter :: r360        = 360.0_r_kind
  real(r_kind),parameter :: r720        = 720.0_r_kind

!--General declarations
  integer(i_kind) :: iret
  integer(i_kind) :: io
  integer(i_kind),dimension(4) :: iadate
  integer(i_kind),dimension(4) :: intdate
  character(4) :: yyyy
  character(2) :: mm
  character(2) :: dd
  character(2) :: hh
  character(2) :: mn
  character(10):: cdate
  
  real(r_kind) :: oberr,temp,temp1,temp2,clock,mu,mu2,it,itt,var_err,std_err,sum_err,var
  integer(i_kind) :: n
  real(r_kind),allocatable :: ob_err(:)
  integer(i_kind),dimension(1) :: seed
  real(r_kind) :: cosazm
  real(r_kind) :: sinazm
  real(r_kind) :: costilt
  real(r_kind) :: sintilt
  real(r_kind) :: cosazm_costilt,sinazm_costilt
  real(r_kind) :: ugesin
  real(r_kind) :: vgesin
  real(r_kind) :: wgesin
  real(r_kind) :: dbzgesin
  real(r_kind) :: zsges
  real(r_kind) :: sin2,termg,termr,termrg,dpres,zob
  real(r_kind) :: b,c,ha,epsh,h,aactual,a43,thistilt
  real(r_kind) :: thistiltr,selev0,celev0,thisrange,thishgt
  real(r_kind) :: celev,selev,gamma,thisazimuthr,thisazimuth,rlon0,rlat0,stahgt, &
                  clat0,slat0,dlat,dlon,thislon,thislat, &
                  rlonloc,rlatloc,rlonglob,rlatglob,rad_per_meter!, &
  real(r_kind),allocatable :: drwpol(:,:,:) !tilt,azm,gate
  real(r_kind)             :: mindrwpol,maxdrwpol
  integer(i_kind) :: irid,itilt,iazm,igate,itime,isig

  character(4) this_staid
  logical   :: diagprint
  logical   :: inside
  logical   :: bufrisopen
  logical   :: radar_location
  logical   :: rite_bufr ! decides whether or not to write the bufr info for a particular radar
  integer   :: diagverbose
  integer(i_kind) :: center_t_window

  !---------DEFAULT SETTINGS---------!
  character(6)    :: datatype = 'NEMSIO'       ! Input data format for reading (NEMSIO)
  logical         :: l4denvar = .true.         ! simulate obs based on the needs of 3d/4d file needs
  integer(i_kind) :: maxobrange=250000_i_kind  ! Range (m) *within* which to use observations 
  integer(i_kind) :: minobrange=20000_i_kind   ! Range (m) *outside* of which
  real(r_kind)    :: mintilt=0.0_r_kind        ! Only use tilt(elevation) angles (deg) >= this number 
  real(r_kind)    :: maxtilt=5.5_r_kind        ! Do no use tilt(elevation) angles (deg) <= this number
  real(i_kind)    :: sigma_err=2.0_r_kind      ! observational error standard deviation 
  real(i_kind)    :: mean_err=2.0_r_kind       ! observational error mean 
  integer(i_kind) :: ithin=4_i_kind            ! Gates to skip for ob thinning (must be >=1)
  character(4)    :: staid='KOUN'              ! default station ID to use
  real(r_kind)    :: mindbz=-999_r_kind        ! minimum dbz value needed at a location to create a drw
  real(r_kind),allocatable :: tilts(:)
  integer(i_kind) :: vcpid=212_i_kind          ! default volume coverage pattern (VCP). 
  integer(i_kind) :: azimuths=360_i_kind       ! number of azimuths
  integer(i_kind) :: gatespc=250_i_kind        ! gate spacing (meters)
  integer(i_kind) :: numgates=400_i_kind       ! number of gates (actual radar space)
  integer(i_kind) :: inumgates=0_i_kind        ! number of gates (to be written when obs are missing minimizes the size of the drw files)
  integer(i_kind) :: ntime=1_i_kind            ! number of times from nc file
  character(20)   :: network="nexrad"          ! default station ID to use
  integer(i_kind) :: nelv=1_i_kind             ! number of elvation tilts
  logical         :: use_dbz=.true.            ! check for dbz at obs location?
  logical         :: use_w=.false.             ! decided to use vertical velocity in observation operator. 
  logical         :: gen_ob_err=.true.         ! logical for generating observation errors 
  logical         :: check_err=.false.         ! check the simulated error mean and standard dev.
  logical         :: rand_err=.false.          ! .false. for reproducible results. errors are still "random" 
  logical         :: test_random_number_gen=.false. ! test the random number generator for simulated ob err. 
  !----------------------------------------------!

  !----------------NEMSIO VARS-------------------!
  real(r_kind),allocatable  ::     ges_u(:,:,:) ! (nlon,nlat,nsig)
  real(r_kind),allocatable  ::     ges_v(:,:,:) ! (nlon,nlat,nsig)
  real(r_kind),allocatable  ::     ges_w(:,:,:) ! (nlon,nlat,nsig)
  real(r_kind),allocatable  ::   ges_dbz(:,:,:) ! (nlon,nlat,nsig)
  !----------------------------------------------!

  ! ALL NEMSIO ADDITIONS HERE BEFORE SORTING 
  type(nemsio_gfile),allocatable,dimension(:) :: gfile
  real(r_kind) :: kap1,kapr
  integer(i_kind) :: nlonsin,nlatsin,nlevsin,idvc
  integer(i_kind) :: idate(4)
  integer(i_kind) :: nfhour
  integer(i_kind) :: ilat
  integer(i_kind) :: ilon
  integer(i_kind) :: k
  integer(i_kind) :: azmspc
  real(r_single),allocatable,dimension(:,:,:) :: nems_vcoord
  real(r_single),allocatable,dimension(:    ) :: lats
  real(r_single),allocatable,dimension(:    ) :: lons
  real(r_single),allocatable,dimension(:    ) :: rlats
  real(r_single),allocatable,dimension(:    ) :: rlons
  real(r_kind),  allocatable,dimension(:,:  ) :: ges_ps
  real(r_kind),  allocatable,dimension(:    ) :: ak
  real(r_kind),  allocatable,dimension(:    ) :: bk
  real(r_kind),  allocatable,dimension(:,:,:) :: ges_prsi
  real(r_kind),  allocatable,dimension(:,:,:) :: ges_prsl
  real(r_kind),  allocatable,dimension(:,:,:) :: ges_lnprsl
  real(r_kind),  allocatable,dimension(:,:,:) :: geop_hgtl
  real(r_kind),  allocatable,dimension(:,:,:) :: ges_t
  real(r_kind),  allocatable,dimension(:,:,:) :: ges_tv
  real(r_kind),  allocatable,dimension(:,:,:) :: ges_q
  real(r_kind),  allocatable,dimension(:,:  ) :: ges_z
  real(r_kind),  allocatable,dimension(:,:  ) :: work
  real(r_kind),  allocatable,dimension(:    ) :: zges 
  real(r_kind),              dimension(nsig) :: hges 
  real(r_kind),              dimension(nsig) :: prsltmp 
  real(r_kind),              dimension(nlon*nlat) :: u1d
  real(r_kind),              dimension(nlon*nlat) :: v1d
  real(r_kind),              dimension(nlon*nlat) :: w1d
  real(r_kind),              dimension(nlon*nlat) :: dbz1d
  real(r_kind),              dimension(nlon*nlat) :: sz1d
  real(r_kind),              dimension(nlon*nlat) :: sp1d
  real(r_kind),              dimension(nlon*nlat) :: t1d
  real(r_kind),              dimension(nlon*nlat) :: q1d
  real(r_kind) :: dlonm1,dlonp1,dlatm1,dlatp1
  integer(i_kind) :: bufrcount


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
  real(r_kind)  :: hdr(14)
  real(r_kind),allocatable :: obs(:,:)
  character(8) :: chdr
  character(8) :: subset
  equivalence (hdr(1),chdr)

  !---------L2RWBUFR CSV TABLE VARS------------!
  character(10) :: message_type
  character(50) :: cdummy

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
  character(len=180)  :: filename3
  character(len=180)  :: filename4
  character(len=180)  :: filename5
  character(len=180)  :: filename6
  character(len=180)  :: filename7
  character(len=180)  :: filename8
  character(len=180)  :: filename9
  character(len=180),allocatable,dimension(:)  :: filename
  character(len=180)  :: namelist_atmfxxx
  character(len=180)  :: fcsthr

  !-------------HEIGHT OF INTERFACE VARS--------------!
  real(r_kind),parameter :: thousand = 1000.0_r_kind
  integer(i_kind) i,j
  real(r_kind) dz,rdog
  real(r_kind),dimension(nsig+1):: height
  real(r_kind) cmpr, x_v, rl_hm, fact, pw, tmp_K, tmp_C, prs_sv, prs_a, ehn_fct, prs_v


  namelist/drw/l4denvar,datatype,ntime,network,staid,ithin,mintilt,maxtilt,maxobrange,minobrange,&
               azimuths,use_dbz,mindbz,gatespc,diagprint,diagverbose,radarcsv,vcpid,use_w

  namelist/simoberr/gen_ob_err,sigma_err,mean_err,check_err,rand_err,test_random_number_gen

  namelist/nio/datapath,filename3,filename4,filename5,filename6,filename7,filename8,filename9
     
!--------------------------------------------------------------------------------------!
!                            END OF ALL DECLARATIONS
!                            !
!--------------------------------------------------------------------------------------!

  !--Set up the constants module used here
  call init_constants_derived
  call init_constants(.false.)    !initialize regional constants
  mindat=ione
  diagprint=.false.
  diagverbose=0
  if(l4denvar) ntime=7 !set here just in case it wasn't in namelist...
  if(.not.l4denvar) ntime=1 !set here just in case it wasn't in namelist...

  !----READ NAMELIST----
  call get_command_argument(1,fcsthr)
  2000 format(a22,a3)
  write(namelist_atmfxxx,2000) "./simnml/namelist.atmf",trim(fcsthr)
  namelist_atmfxxx=trim(namelist_atmfxxx)
  open(11,file=namelist_atmfxxx)
  read(11,drw)
  read(11,simoberr)
  if(datatype == 'NEMSIO') read(11,nio)
  !----READ NAMELIST----

  !----RANDOM NUMBER GENERATOR
  if(test_random_number_gen) then
     oberr=0.0_r_kind
     sum_err=0.0_r_kind
     var_err=0.0_r_kind
     std_err=0.0_r_kind
     mu=0.0_r_kind
     itt=144000
     do it=1,itt
        temp=0_r_kind
        temp1=0_r_kind
        temp2=0_r_kind
        if(rand_err) then
           call cpu_time(clock)
           seed(1)=int(clock*500000)**2
           call random_seed(put=seed)
           call random_number(temp1)
           call cpu_time(clock)
           seed(1)=int(clock*600000)**2
           call random_seed(put=seed)
           call random_number(temp2)
        else
           seed(1)=(abs(cos(it))* (10**8) )
           if(itt <= 100) write(6,*) seed 
           call random_seed(put=seed)
           call random_number(temp1)
           seed(1)=(abs(sin(it))* (10**6) )
           if(itt <= 100) write(6,*) seed 
           call random_seed(put=seed)
           call random_number(temp2)
        end if
        temp = sigma_err * sqrt( -2*log(temp1) ) * sin(2*PI*temp2) + mean_err
        if(itt <= 100) write(6,*) temp1,temp2,temp
        mu = mu + temp
        var = var + (temp - mean_err)**2
     enddo
     mu = mu/itt
     mu2=sqrt(var/(itt-1))
     write(6,*) "mean and stdev=",mu,mu2
     write(6,*) "mean and stdev=",mean_err,sigma_err
     stop 
   end if ! end test of random number generator
 
  !----CREATE VOLUME COVERAGE PATTERN (VCP)----start 
  !----SIMPLIFICATION NOTES:
  !    Most precipitation vcps have multiple scans at each layer.
  !    Further, there are SAILS1,2,&3 which jump back to 0.5 every so often.
  !    All simulated obs are assumed to be valid at the background timestamp.
  !    I cannot include this, because I do not have minutely output to do so.
  if(vcpid == 215) then
     nelv=14
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
  allocate(filename(ntime))
  do itime=1,ntime
     if(l4denvar) then
        if(itime==1) filename(itime)=trim(datapath) // trim(filename3)
        if(itime==2) filename(itime)=trim(datapath) // trim(filename4)
        if(itime==3) filename(itime)=trim(datapath) // trim(filename5)
        if(itime==4) filename(itime)=trim(datapath) // trim(filename6)
        if(itime==5) filename(itime)=trim(datapath) // trim(filename7)
        if(itime==6) filename(itime)=trim(datapath) // trim(filename8)
        if(itime==7) filename(itime)=trim(datapath) // trim(filename9)
    end if
    if(.not.l4denvar) then
        if(itime==1) filename(itime)=trim(datapath) // trim(filename6)
    end if
  end do

  if(diagprint .and. diagverbose >= 1) then
     write(6,*) "----NEMSIO FILE SPECS--------"
     write(6,*) "nlat",nlat
     write(6,*) "nlon",nlon
     write(6,*) "nsig",nsig
     1001 format (a16,f5.2,a3)
     write(6,1001) "grid spacing = ~",360.0*111.0/nlon,"-km"

     write(6,*) "----Namelist settings----"
     write(6,*) "radarcsv:   ",radarcsv
     write(6,*) "diagprint:  ",diagprint
     write(6,*) "diagverbose:",diagverbose
     write(6,*) "staid       ",staid
     write(6,*) "tilts       ",tilts
     write(6,*) "azimuths    ",azimuths
     write(6,*) "numgates    ",numgates
     write(6,*) "ithin       ",ithin
     write(6,*) "gatespc     ",gatespc
     write(6,*) "maxobrange  ",maxobrange
     write(6,*) "minobrange  ",minobrange
     write(6,*) "mintilt     ",mintilt
     write(6,*) "maxtilt     ",maxtilt
     write(6,*) "mindbz      ",mindbz
     write(6,*) "datapath:   ",datapath
     do itime=1,ntime
        write(6,*) "nemsio file(",itime,"):  ",filename(itime)
     enddo
  end if
  if(modulo(360,azimuths) /= 0) stop "ERROR: !!! AZIMUTHS IN NAMELIST IS NOT A FACTOR OF 360-DEG. !!!"

  if(diagprint .and. diagverbose >= 3) then
     write(6,*) "deg2rad",deg2rad
     write(6,*) "rearth",rearth
     if(deg2rad>0) write(6,*) "(derived constants succesfully initialized)"
     if(rearth==6370.e03_r_kind) write(6,*) "(regional consts sucess init)"
  end if

  if (minobrange >= maxobrange) then
     STOP 'MINIMUM OB RANGE >= MAXIMUM OB RANGE. PROGRAM STOPPING NOW drwsim.f90'
  end if
  if (ithin < 1) then
     STOP 'ithin MUST BE >=1! CHECK NAMELIST AND RESET.'
  end if


!------ OPEN GLOBAL RADAR LIST ---------------------------
  write(6,*) "Reading Global Radar List"
  numradars=0

  !-Initial read to count the number of radars in the file list.
  open(41,file=trim(radarcsv))
  read(41,*) cdummy !read 1st line which is just a header.
  do !ii=1,numradars
     read(41,*,iostat=io) cdummy 
     if(io == 0) then
       numradars=numradars+1
     elseif(io < 0) then
       exit !End of file
     end if
  end do
  close(41)
  write(6,*) "Number of radars = ",numradars

  !-Second read to store the radar metadata
  allocate(dfid(numradars),dflat(numradars),dflon(numradars),dfheight(numradars))
  open(40,file=trim(radarcsv))
  read(40,*) cdummy !read 1st line which is just a header.
  do ii=1,numradars
     read(40,*) dfid(ii),dflat(ii),dflon(ii),dfheight(ii)
     dfid(ii)=trim(dfid(ii))

     if(diagprint .and. diagverbose >= 1) then
        1002 format(a5,1x,f15.8,1x,f15.8,1x,f15.8)
        write(6,1002) dfid(ii),dflat(ii),dflon(ii),dfheight(ii)
     end if
  end do
  close(40)
  write(6,*) "Done: Reading Global Radar List"

!------ OPEN NEMSIO FILE FOR READING -----------------------
!  if(datatype=='NEMSIO') then
  write(6,*) 'Open and Read NEMSIO files'
  allocate(work(nlon,nlat))
  work    = zero
  if(l4denvar) ntime=7
  if(.not.l4denvar) ntime=1 !set here just in case it wasn't in namelist...
  allocate(gfile(ntime))

  !   do itime=1,ntime
  ndata=izero
  nradars=izero

  !--Call Timer
  call date_and_time(values=time_array_0)
  start_time = time_array_0(5)*3600 + time_array_0(6)*60 + time_array_0(7) + time_array_0(8)*0.001
  write(6,*) 'STARTING RADAR WIND SIMULATION:'

  bufrcount=0
  call nemsio_init(iret=iret)
  center_t_window=(ntime-1)/2 + 1
  call nemsio_open(gfile(center_t_window),filename(center_t_window),'READ',iret=iret)
  call nemsio_getheadvar(gfile(center_t_window),'idate',idate,iret)
  call nemsio_getheadvar(gfile(center_t_window),'nfhour',nfhour,iret)
  call newdate(idate,nfhour,iadate)
  call nemsio_close(gfile(center_t_window))
  call w3ai15(iadate(1),yyyy,1,4,'')
  call w3ai15(iadate(2),  mm,1,2,'')
  call w3ai15(iadate(3),  dd,1,2,'')
  call w3ai15(iadate(4),  hh,1,2,'')
  call w3ai15(iadate(5),  mn,1,2,'')
  cdate=trim(yyyy)//trim(mm)//trim(dd)//trim(hh)

  loopOVERtime: do itime=1,ntime
     call nemsio_init(iret=iret)
     write(6,*) "filename(",itime,")=",filename(itime)
     call nemsio_open(gfile(itime),filename(itime),'READ',iret=iret)
     call nemsio_getheadvar(gfile(itime),'idate',idate,iret)
     call nemsio_getheadvar(gfile(itime),'nfhour',nfhour,iret)
     call nemsio_getfilehead(gfile(itime),iret=iret,dimx=nlonsin,dimy=nlatsin,dimz=nlevsin,idvc=idvc)

     call newdate(idate,nfhour,iadate)
     write(6,*) 'iadate=',iadate

     ! rlats rlons
     write(6,*) 'Reading lat lon'
     if(allocated(lats) ) deallocate(lats )
     if(allocated(lons) ) deallocate(lons )
     if(allocated(rlats)) deallocate(rlats)
     if(allocated(rlons)) deallocate(rlons)
     allocate(lats(nlon*nlat))
     allocate(lons(nlon*nlat))
     allocate(rlats(nlat)    )
     allocate(rlons(nlon)    )
     call nemsio_getfilehead(gfile(itime),iret=iret,lon=lons,lat=lats)
     do ilon=1,nlon
        rlons(ilon)=lons(ilon)*deg2rad
     end do
     do ilat=1,nlat
        rlats(ilat)=lats(ilat*nlon-(nlon-1))*deg2rad 
     end do
     !lats indicies 1-3072 are all the same and are 89.9
     !lons indicies 1-3072 are increase from 0 to 359.88

     ! read 3d fields (u,v,w,dbz)
     !itime=1
     if(allocated(ges_u)   )  deallocate(ges_u  )
     if(allocated(ges_v)   )  deallocate(ges_v  )
     if(allocated(ges_w)   )  deallocate(ges_w  )
     if(allocated(ges_dbz) )  deallocate(ges_dbz)
     if(allocated(ges_t)   )  deallocate(ges_t  )
     if(allocated(ges_tv)  )  deallocate(ges_tv )
     if(allocated(ges_q)   )  deallocate(ges_q  )
     allocate(ges_u(nlon,nlat,nsig)  )
     allocate(ges_v(nlon,nlat,nsig)  )
     allocate(ges_w(nlon,nlat,nsig)  )
     allocate(ges_dbz(nlon,nlat,nsig))
     allocate(ges_t(nlon,nlat,nsig)  )
     allocate(ges_tv(nlon,nlat,nsig) )
     allocate(ges_q(nlon,nlat,nsig)  )
     ges_u(:,:,:)=0.0_r_kind ! initialize
     ges_v(:,:,:)=0.0_r_kind ! initialize
     ges_w(:,:,:)=0.0_r_kind ! initialize
     ges_dbz(:,:,:)=-999.0_r_kind ! initialize
     ges_q(:,:,:)=0.0_r_kind !initialize
     ges_t(:,:,:)=0.0_r_kind !initialize
     ges_tv(:,:,:)=0.0_r_kind !initialize
    

     write(6,*) 'Reading u,v,w,q,t,dbz level by level'
     do isig=1,nsig

        ! u
        call nemsio_readrecv(gfile(itime),'ugrd','mid layer',isig,u1d,iret=iret)
        if(iret/=0) stop "ERROR: ugrd"
        ges_u(:,:,isig)=reshape(u1d(:),(/size(work,1),size(work,2)/))

        ! v
        call nemsio_readrecv(gfile(itime),'vgrd','mid layer',isig,v1d,iret=iret)
        if(iret/=0) stop "ERROR: vgrd"
        ges_v(:,:,isig)=reshape(v1d(:),(/size(work,1),size(work,2)/))

        ! w
        if(use_w) then
           call nemsio_readrecv(gfile(itime),'dzdt','mid layer',isig,w1d,iret=iret) !FV3=dzdt; NMMB=w_tot
           if(iret/=0) stop "ERROR: dzdt"
           ges_w(:,:,isig)=reshape(w1d(:),(/size(work,1),size(work,2)/))
        endif
!
        ! dbz
        if(use_dbz) then !even if this is false, dbz will be initialized to missing everywhere 
           call nemsio_readrecv(gfile(itime),'dbz','mid layer',isig,dbz1d,iret=iret)
           if(iret/=0) stop "ERROR: dbz"
           ges_dbz(:,:,isig)=reshape(dbz1d(:),(/size(work,1),size(work,2)/))
        endif

        ! q 
        call nemsio_readrecv(gfile(itime),'spfh','mid layer',isig,q1d,iret=iret)
        if(iret/=0) stop "ERROR: spfh"
        ges_q(:,:,isig)=reshape(q1d(:),(/size(work,1),size(work,2)/))

        ! t => convert to tv later
        call nemsio_readrecv(gfile(itime),'tmp','mid layer',isig,t1d,iret=iret)
        if(iret/=0) stop "ERROR: tmp"
        ges_t(:,:,isig)=reshape(t1d(:),(/size(work,1),size(work,2)/))

     end do

     ! read sfc height
     if(itime==1) then ! surface height is constant at all times.
        write(6,*) 'Reading surface height'
        allocate(ges_z(nlon,nlat))
        ges_z=0.0_r_kind ! initialize
        call nemsio_readrecv(gfile(itime),'hgt','sfc',1,sz1d,iret=iret)
        ges_z(:,:)=reshape(sz1d(:),(/size(work,1),size(work,2)/))
     end if

     ! read sfc pressure
     write(6,*) 'Reading surface pressure'
     if( allocated(ges_ps) ) deallocate(ges_ps)
     allocate(ges_ps(nlonsin,nlatsin))
     ges_ps=0.0_r_kind ! initialize
     call nemsio_readrecv(gfile(itime),'pres','sfc',1,sp1d,iret=iret)
     sp1d=0.001_r_kind*sp1d
     ges_ps(:,:)=reshape(sp1d(:),(/size(work,1),size(work,2)/))

     ! ak bk --> interface pressure/height calculation
     write(6,*) 'Reading ak bk'
     if(allocated(nems_vcoord)) deallocate(nems_vcoord)
     if(allocated(ges_prsi)   ) deallocate(ges_prsi   )
     if(allocated(ges_prsl)   ) deallocate(ges_prsl   )
     if(allocated(ges_lnprsl) ) deallocate(ges_lnprsl )
     if(allocated(ak)         ) deallocate(ak         )
     if(allocated(bk)         ) deallocate(bk         )
     if(allocated(zges)       ) deallocate(zges       )
     if(allocated(geop_hgtl)  ) deallocate(geop_hgtl  )
     allocate(nems_vcoord(nsig+1,3,2)   )
     allocate(ges_prsi(nlon,nlat,nsig+1))
     allocate(ges_prsl(nlon,nlat,nsig)  )
     allocate(ges_lnprsl(nlon,nlat,nsig))
     allocate(ak(nsig+1)                )
     allocate(bk(nsig+1)                )
     allocate(zges(nsig)                )
     allocate(geop_hgtl(nlon,nlat,nsig ))

     zges = 0.0_r_kind ! initialize
     call nemsio_getfilehead(gfile(itime),iret=iret,vcoord=nems_vcoord)
     if ( idvc == 2 ) then      ! hybrid coordinate
        ak = nems_vcoord(1:nsig+1,1,1) ! convert to mb
        bk = nems_vcoord(1:nsig+1,2,1)
     endif

     ! This is how you compute the pressure at interfaces (not with delz/delp)
     ! guess_grid.F90 @ 1068
     do isig=1,nsig+1
        ges_prsi(:,:,isig)=0.001_r_kind*ak(isig)+(bk(isig)*ges_ps) ! pressure at interfaces
     enddo
     deallocate(ak,bk)
!guess_grids.F90 @ 1140 subroutine load_prsges
     kap1=rd_over_cp+one
     kapr=one/rd_over_cp
     do j=1,nlat
        do i=1,nlon
           do isig=1,nsig
              !load mid-layer pressure by using phillips vertical interpolation
              ges_prsl(i,j,isig)=((ges_prsi(i,j,isig)**kap1-ges_prsi(i,j,isig+1)**kap1)/&
                    (kap1*(ges_prsi(i,j,isig)-ges_prsi(i,j,isig+1))))**kapr
           end do
        end do
     enddo
!guess_grids.F90 @ 1392 subroutine load_geop_hgt     
     do i=1,nlon
        do j=1,nlat
           !--Compute virtual temperature from dry temp and specific humididty.
           do k=1,nsig
              ges_tv(i,j,k) = ges_t(i,j,k) * (one + fv * (ges_q(i,j,k)))
           enddo
        enddo
     enddo
     !-Compute geopotential heights at mid layer
     write(6,*) "Computing geopotential height at mid layer",itime,"/",ntime,": started"
     rdog = rd/grav
     do j=1,nlat     !LIPPI Flipped the order of nlon/nlat
        do i=1,nlon
           k=1
           fact     = one + fv * ges_q(i,j,k)
           pw       = eps + ges_q(i,j,k)*( one - eps )
           tmp_K    = ges_tv(i,j,k) / fact
           tmp_C    = tmp_K - t0c
           prs_sv   = exp(psv_a*tmp_K**2 + psv_b*tmp_K + psv_c + psv_d/tmp_K)  !Pvap sat, eq A1.1 (Pa)
           prs_a   = thousand * exp(half*(log(ges_prsi(i,j,k)) + log(ges_prsl(i,j,k))))     ! (Pa)
           ehn_fct = ef_alpha + ef_beta*prs_a + ef_gamma*tmp_C**2 ! enhancement factor (eq. A1.2)
           prs_v   = ges_q(i,j,k) * prs_a / pw   ! vapor pressure (Pa)
           rl_hm   = prs_v / prs_sv    ! relative humidity
           x_v     = rl_hm * ehn_fct * prs_sv / prs_a     ! molar fraction of water vapor (eq. A1.3)
          ! Compressibility factor (eq A1.4 from Picard et al 2008)
           cmpr = one - (prs_a/tmp_K) * (cpf_a0 + cpf_a1*tmp_C + cpf_a2*tmp_C**2 &
                           + (cpf_b0 + cpf_b1*tmp_C)*x_v + (cpf_c0 + cpf_c1*tmp_C)*x_v**2 ) &
                           + (prs_a**2/tmp_K**2) * (cpf_d + cpf_e*x_v**2)

           h        = rdog * ges_tv(i,j,k)
           dz       = h * cmpr * log(ges_prsi(i,j,k)/ges_prsl(i,j,k))
           height(k)= ges_z(i,j) + dz

           do k=2,nsig
              fact     = one + fv * half * (ges_q(i,j,k-1)+ges_q(i,j,k))
              pw       = eps + half * (ges_q(i,j,k-1)+ges_q(i,j,k))*( one - eps )
              tmp_K    = half * (ges_tv(i,j,k-1)+ges_tv(i,j,k)) / fact
              tmp_C    = tmp_K - t0c
              prs_sv   = exp(psv_a*tmp_K**2 + psv_b*tmp_K + psv_c + psv_d/tmp_K) !Pvap sat, eq A1.1 (Pa)
              prs_a   = thousand * exp(half*(log(ges_prsl(i,j,k-1)) + log(ges_prsl(i,j,k))))     ! (Pa)
              ehn_fct = ef_alpha + ef_beta*prs_a + ef_gamma*tmp_C**2 ! enhancement factor (eq. A1.2)
              prs_v   = half * (ges_q(i,j,k-1)+ges_q(i,j,k)) * prs_a / pw   ! vapor pressure (Pa)
              rl_hm   = prs_v / prs_sv    ! relative humidity
              x_v     = rl_hm * ehn_fct * prs_sv / prs_a     ! molar fraction of water vapor (eq. A1.3)
             ! Compressibility factor (eq A1.4 from Picard et al 2008)
              cmpr = one - (prs_a/tmp_K) * (cpf_a0 + cpf_a1*tmp_C + cpf_a2*tmp_C**2 &
                              + (cpf_b0 + cpf_b1*tmp_C)*x_v + (cpf_c0 + cpf_c1*tmp_C)*x_v**2 ) &
                              + (prs_a**2/tmp_K**2) * (cpf_d + cpf_e*x_v**2)

              h         = rdog * half * (ges_tv(i,j,k-1)+ges_tv(i,j,k))
              dz        = h * cmpr * log(ges_prsl(i,j,k-1)/ges_prsl(i,j,k))
              height(k) = height(k-1) + dz
           enddo
           do k=1,nsig
              geop_hgtl(i,j,k) = height(k) - ges_z(i,j)
           enddo
        enddo
     enddo 


     write(6,*) "Computing geopotential height at mid layer",itime,"/",ntime,": done"
     write(6,*) idate,nfhour
     write(6,*) 'Done: Reading NEMSIO files'
     deallocate(lats,lons)
     deallocate(ges_prsi,ges_prsl,ges_tv,ges_t,ges_q)

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

     if(check_err) allocate(ob_err(nelv*360*numgates*numradars))
     loopOVERradars: do irid=1,numradars 
        allocate(drwpol(nelv,360,numgates))
        drwpol=-999.0_r_kind !Initialize/Reset the drw polar field
        mu=0.0_r_kind
        mindrwpol=huge_single
        maxdrwpol=-huge_single
        radar_location=.true. ! logical to only compute radar x,y once later in loop - preset to true.
        this_staid=adjustl(trim(dfid(irid)))
        ifKGRK: if(this_staid==trim(adjustl(staid)) .or. trim(adjustl(staid))=='all') then
           stahgt=nint(dfheight(irid)*100.0_r_kind)/100.0_r_kind
           rlon0=nint(dflon(irid)*100.0_r_kind)/100.0_r_kind*deg2rad !round to nearest 100th and convert to rads
           rlat0=nint(dflat(irid)*100.0_r_kind)/100.0_r_kind*deg2rad
           clat0=cos(rlat0)
           slat0=sin(rlat0)
           loopOVERtilts:    do itilt=1,nelv
              if(tilts(itilt) > maxtilt .or. tilts(itilt) < mintilt) then
                 cycle
              end if
              bufrisopen=.false.   !Initialize bufr file as closed.
              thistilt=tilts(itilt)
              thistiltr=thistilt*deg2rad
              celev0=cos(thistiltr)
              selev0=sin(thistiltr)
              loopOVERazimuths: do iazm=azmspc,360,azmspc !360=>90;  90=>0; 180=>270;  270=>180
                 1000 format(a5,1x,i4,i2.2,i2.2,i2.2,&
                          3x,a6,1x,i3,a2,i3,1x,a4,&  !3x,a6,1x,a4,&
                          3x,a5,1x,i2,a2,i2,1x,a1,f4.1,a1,f4.1,a1,&
                          3x,a4,1x,i3,a4,&
                          3x,a4,1x,i3,1x,i2,a7,i5,a3) !3x,a4,1x,i3,&
                          !3x,a6,1x,i2,a1,i5,a9)
                 if(iazm == 360) then !print last - less print should speed it up.
                    write(6,1000),"Date:",iadate(1),iadate(2),iadate(3),iadate(4),&
                                  "Radar:",irid,"of",numradars,adjustl(trim(dfid(irid))),&
                                  "Tilt:",itilt,"of",nelv,"(",tilts(itilt),"/",tilts(nelv),")",&
                                  "Azm:",iazm,"/360",&
                                  "VCP:",vcpid,azmspc,"[deg] x",gatespc,"[m]" !&
                                  !"ObRes:",azmspc,"x",gatespc,"[deg x m]"
                 endif
                 thisazimuth=90.0_r_kind-float(iazm) ! 90-azm to be consistent with l2rwbufr
                 if(thisazimuth>=r360) thisazimuth=thisazimuth-r360
                 if(thisazimuth<zero) thisazimuth=thisazimuth+r360
                 thisazimuthr=thisazimuth*deg2rad
                 loopOVERgates: do igate=1,numgates     
                    inside=.false. ! is our ob location inside the bounds? preset to false, then check.
                    if(igate*gatespc >= minobrange .and. igate*gatespc <= maxobrange) inside=.true.
                    ifinside: if(inside) then
! read_l2bufr_mod.f90 @ 686 subroutine radar_bufr_read_all
                       !--Find observation height using method from read_l2bufr_mod.f90 
                       thisrange=(igate)*gatespc
                       aactual=(rearth+stahgt)
                       a43=aactual*four_thirds
                       b   = thisrange*(thisrange+two*aactual*selev0)
                       c   = sqrt(aactual*aactual+b)
                       ha  = b/(aactual+c)
                       epsh=(thisrange*thisrange-ha*ha)/(r8*aactual)
                       h=ha-epsh
                       thishgt=stahgt+h
                       dpres=thishgt !store the absolute ob height (m) in dpres.                   
                       !-Get corrected tilt angle @ 715
                       celev=celev0
                       selev=selev0
                       celev=a43*celev0/(a43+h)
                       selev=(thisrange*thisrange+h*h+two*a43*h)/(two*thisrange*(a43+h))
                       !corrected_tilt=atan2(selev,celev)*rad2deg
                       !thistilt=corrected_tilt
                       !thistiltr=thistilt*deg2rad   
                       gamma=half*thisrange*(celev0+celev)
                       gamma=thisrange

                       !-Get earth lat lon of superob @ 729
                       rad_per_meter=one/rearth
                       rlonloc=rad_per_meter*gamma*cos(thisazimuthr)
                       rlatloc=rad_per_meter*gamma*sin(thisazimuthr)
                       call invtllv(rlonloc,rlatloc,rlon0,clat0,slat0,rlonglob,rlatglob)
                       thislat=rlatglob*rad2deg
                       thislon=rlonglob*rad2deg

                       if(thislon>=r360) thislon=thislon-r360
                       if(thislon<zero) thislon=thislon+r360
                       dlat=thislat
                       dlon=thislon
                       !--Find grid relative location of the ob. !!call tll2xy(thislon,thislat,dlon,dlat)
                       call grdcrd1(dlat,rlats*rad2deg,nlat,-1) !lats are in descending order
                       call grdcrd1(dlon,rlons*rad2deg,nlon, 1)

! read_l2bufr_mod.f90 @ 740 subroutine radar_bufr_read_all
                       !--Get corrected azimuth
                       !clat1=cos(rlatglob)
                       !caz0=cos(thisazimuthr)
                       !saz0=sin(thisazimuthr)
                       !cdlon=cos(rlonglob-rlon0)
                       !sdlon=sin(rlonglob-rlon0)
                       !caz1=clat0*caz0/clat1
                       !saz1=saz0*cdlon-caz0*sdlon*slat0
                       !corrected_azimuth=atan2(saz1,caz1)*rad2deg
                       !delazmmax=max(min(abs(corrected_azimuth-thisazimuth-r720),&
                       !                  abs(corrected_azimuth-thisazimuth-r360),&
                       !                  abs(corrected_azimuth-thisazimuth     ),&
                       !                  abs(corrected_azimuth-thisazimuth+r360),&
                       !                  abs(corrected_azimuth-thisazimuth+r720)),delazmmax)
                       !thisazimuth=corrected_azimuth
                       !thisazimuthr=thisazimuth*deg2rad

                       !--Interpolate surface height to grid relative ob location (dlon,dlat). 
                       dlonm1=floor(dlon)-1
                       dlonp1=floor(dlon)+1
                       dlatm1=floor(dlat)-1
                       dlatp1=floor(dlat)+1
                       zsges=0.0_r_kind
                       prsltmp=0.0_r_kind
                       hges=0.0_r_kind

                       ! surface height (zsges): ges_z(nlon,nlat) => zsges
                       ! setuprw.f90 @ 359
                       call tintrp2a_single_level_sliced(ges_z(dlonm1:dlonp1,dlatm1:dlatp1),&
                                                   zsges,dlon-floor(dlon)+2,dlat-floor(dlat)+2)

                       ! geopotential height at mid layers: geop_hgtl(nlon,nlat,nsig,ntime) => hges(nsig)
                       call tintrp2a_sliced(geop_hgtl(dlonm1:dlonp1,dlatm1:dlatp1,:),&
                                                   hges(:),dlon-floor(dlon)+2,dlat-floor(dlat)+2,nsig)
                      
                       !--Remove terrain height from ob absolute height and reject if below ground.
                       if(zsges>=dpres) then
                         !write(6,*) 'zsges =',zsges,'is greater than dpres',dpres,'. Rejecting ob.'
                         cycle
                       end if
                       dpres=dpres-zsges

!**********************************CONVERT GEOP HEIGHT TO GEOM HEIGHT****!
                       sin2  = sin(thislat*deg2rad)*sin(thislat*deg2rad)
                       termg = grav_equator * ((one+somigliana*sin2)/sqrt(one-eccentricity*eccentricity*sin2))
                       termr = semi_major_axis /(one + flattening + grav_ratio - two*flattening*sin2)
                       termrg = (termg/grav)*termr
                       do k=1,nsig
                          zges(k) = (termr*hges(k)) / (termrg-hges(k))  ! eq (23)
                       end do
                        
                    !  Convert observation height (in dpres) from meters to grid relative units.
                    !  Save the observation height in zob for later use.
                       zob = dpres
                       call grdcrd1(dpres,zges,nsig,1) ! get grid coordinates of dpres from zges
                             
!**********************************CONVERT GEOP HEIGHT TO GEOM HEIGHT****!

                       !--Interpolate guess dbz to observation location - cycle if below threshold.
                       call tintrp3(ges_dbz(:,:,:),dbzgesin,dlon,dlat,dpres)
                       dbzCheck: if(dbzgesin >= mindbz .or..not. use_dbz) then
                          !--Interpolate guess wind to observation location                                  
                          call tintrp3(ges_u(:,:,:),ugesin,dlon,dlat,dpres)
                          call tintrp3(ges_v(:,:,:),vgesin,dlon,dlat,dpres)
                          if(use_w) call tintrp3(ges_w(:,:,:),wgesin,dlon,dlat,dpres)

                          !--Convert guess u,v,w wind components to radial value               
                          cosazm  = cos(thisazimuthr)! cos(azimuth angle)                       
                          sinazm  = sin(thisazimuthr)! sin(azimuth angle)                       
                          costilt = cos(thistiltr)   ! cos(tilt angle)
                          if(use_w) sintilt = sin(thistiltr)   ! sin(tilt angle)
                          cosazm_costilt = cosazm*costilt
                          sinazm_costilt = sinazm*costilt
                          !-------------WIND FORWARD MODEL-----------------------------------------!
                          !drwpol(itilt,iazm,igate) = ugesin*cosazm_costilt +vgesin*sinazm_costilt +wgesin*sintilt
                          drwpol(itilt,iazm,igate) = ugesin*cosazm_costilt +vgesin*sinazm_costilt
                          if(use_w) drwpol(itilt,iazm,igate) = drwpol(itilt,iazm,igate) +wgesin*sintilt
                          ndata=ndata+1
                          !-------------ADD OB ERR-------------------------------------------------!
                          if(gen_ob_err) then
                             temp1=0_r_kind
                             temp2=0_r_kind
                             oberr=0_r_kind
                             if(rand_err) then
                                call cpu_time(clock)
                                seed(1)=int(clock*500000)**2
                                call random_seed(put=seed)
                                call random_number(temp1)
                                call cpu_time(clock)
                                seed(1)=int(clock*600000)**2
                                call random_seed(put=seed)
                                call random_number(temp2)
                             else
                                 seed(1)=(abs(cos(real(ndata)))* (10**8) )
                                 call random_seed(put=seed)
                                 call random_number(temp1)
                                 seed(1)=(abs(sin(real(ndata)))* (10**6) )
                                 call random_seed(put=seed)
                                 call random_number(temp2)
                             end if
                             oberr = sigma_err * sqrt( -2*log(temp1) ) * sin( 2*PI*temp2 ) + mean_err 
                             drwpol(itilt,iazm,igate) = drwpol(itilt,iazm,igate) + oberr

                             if(check_err) then
                                ob_err(ndata) = oberr
                                sum_err = sum_err + ob_err(ndata)
                             end if

                          end if !generate ob error

                          !round to nearest 10th since this is automatically done when writing to bufr.
                          drwpol(itilt,iazm,igate) = nint(drwpol(itilt,iazm,igate)*10.0_r_kind)/10.0_r_kind

                          if(diagprint .and. diagverbose >= 1 .and. drwpol(itilt,iazm,igate) /= -999) then
                             if(drwpol(itilt,iazm,igate) < mindrwpol) then
                                mindrwpol=drwpol(itilt,iazm,igate)
                             end if
                             if(drwpol(itilt,iazm,igate) > maxdrwpol) then
                                maxdrwpol=drwpol(itilt,iazm,igate)
                             end if
                          end if 

                       end if dbzCheck
                    end if ifinside
                 end do loopOVERgates
              end do loopOVERazimuths
           end do loopOVERtilts

           if(check_err) then
             mu = sum_err / ndata
             do n=1,ndata
                var_err = var_err + ( ob_err(n) - mu )**2
             end do
             var_err = var_err / (ndata-1)
             std_err = sqrt( var_err )
             
             write(6,*) 'mean and stadard deviation of oberr=',mu,std_err,ndata,var_err
             write(6,*) 'what the mean and std dev should be=',mean_err,sigma_err 
             deallocate(ob_err)
             stop
           end if

           if(diagprint .and. diagverbose >= 1) write(6,*) "min/max drw: ",mindrwpol,maxdrwpol
           if(mindrwpol <= maxdrwpol) rite_bufr=.true.
           if(mindrwpol >  maxdrwpol) rite_bufr=.false.


           !-------------BUFFERIZE--------------------------------------------------!
           !    At this point we have observations from all radars at every scan
           ! angle at a single time. We will put this information in its own bufr
           ! file hence this is contained within loopOVERtime.
           !
           ifrite_bufr: if(rite_bufr) then !if all obs are missing, skip this radar
              write(6,*)"Writing bufr file for ",trim(dfid(irid))
              hdstr='SSTN CLON CLAT HSMSL HSALG ANEL YEAR MNTH DAYS HOUR MINU SECO QCRW ANAZ'
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
              if(diagprint .and. diagverbose >= 2) write(6,*) message_type
              write(6,*) 'iadate',iadate
              !call w3ai15(iadate(1),yyyy,1,4,'')
              !call w3ai15(iadate(2),  mm,1,2,'')
              !call w3ai15(iadate(3),  dd,1,2,'')
              !call w3ai15(iadate(4),  hh,1,2,'')
              !call w3ai15(iadate(5),  mn,1,2,'')
              !cdate=trim(yyyy)//trim(mm)//trim(dd)//trim(hh)
              subset=trim(adjustl(message_type))
              chdr   = dfid(irid)       !SSTN - RADAR STATION IDENTIFIER -- uses same memory location as hdr(1)
              hdr(2) = dflon(irid)      !CLON - LONGITUDE (COARSE ACCURACY)
              hdr(3) = dflat(irid)      !CLAT - LATITUDE (COARSE ACCURACY)
              hdr(4) = dfheight(irid)   !SELV - HEIGHT OF STATION
              hdr(5) = 00 
             !hdr(6) - tilt loop below.
              hdr(7) = iadate(1)  !YEAR - YEAR
              hdr(8) = iadate(2)  !MNTH - MONTH
              hdr(9) = iadate(3)  !DAYS - DAY
              hdr(10) = iadate(4) !HOUR - HOUR 
              hdr(11)= 00               !MINU - MINUTE
              hdr(12)= 00               !SECO - SECONDS
              hdr(13)= 1                !QCRW - QUALITY MARK FOR WINDS ALONG RADIAL LINE
             !hdr(14)- azm loop below.
              bufrtilt: do itiltbufr=1,nelv
                 intdate=iadate(1)*1000000 + iadate(2)*10000 +iadate(3)*100 + iadate(4) ! int(yyyymmddhh)
                 hdr(6) = tilts(itiltbufr) 

                 if(.not.bufrisopen) then !open a new message for each radar 
                    write(6,*) "intdate",intdate
                    write(6,*) "cdate",cdate
                    bufrfilename='./simbufr/'//trim(adjustl(network))//'_'//trim(cdate)//'_fv3.t'//trim(hh)//'z_drw.bufr'
                    write(6,*) "bufr file name is:",bufrfilename
                    open(unit=11,file='l2rwbufr.table',status='old',action='read',form='formatted')

                    if(bufrcount == 0) then
                       open(unit=10,file=trim(bufrfilename),status='unknown',action='write',form='unformatted')
                       call openbf(10,'OUT',11)
                    else ! APPEND MESSAGES AFTER THE FIRST IS WRITTEN
                       open(unit=10,file=trim(bufrfilename),status='old',    action='write',form='unformatted')
                       call openbf(10,'APN',11)
                    endif
                    bufrisopen=.true.
                    bufrcount=bufrcount+1
                 end if
   
                 call openmb(10,trim(subset),intdate)
                 bufrazm: do iazmbufr=360/azimuths,360,360/azimuths
                    iazmbufr90=90-iazmbufr
                    if(iazmbufr90>=r360) iazmbufr90=iazmbufr90-r360
                    if(iazmbufr90< zero) iazmbufr90=iazmbufr90+r360
                    allocate(obs(3,numgates))
                    obs=-999.0_r_kind ! Initialize as missing values
                    hdr(14)=float(iazmbufr90)
                    inumgates=0
                    bufrgate: do igatebufr=1,numgates
                       obs(1,igatebufr) = igatebufr !DISTANCE (FROM ANTENNA TO GATE CENTER) IN UNITS OF 250M
                       obs(2,igatebufr) = drwpol(itiltbufr,iazmbufr90,igatebufr) !DOPPLER MEAN RADIAL VELOC 
                       obs(3,igatebufr) = 1.0_r_kind                       !DOPPLER VELOCITY SPECTRAL WIDTH
                    end do bufrgate
                    ! encode radial velocity
                    call ufbint(10,hdr,14,1,iret,trim(hdstr))
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
              nradars=nradars+1
           else  !ifrite_bufr
              write(6,*)"Nothing to write for radar id: ",trim(dfid(irid))
           end if ifrite_bufr
        end if ifKGRK
        deallocate(drwpol) !must always be deallocated before processing next radar
     end do loopOVERradars 
  end do loopOVERtime

  !-Call Timer
  call date_and_time(values=time_array_1)
  finish_time = time_array_1(5)*3600 + time_array_1(6)*60 + time_array_1(7) + time_array_1(8)*0.001
  total_time=finish_time-start_time
  hrs =int(        total_time/3600.0       )
  mins=int(    mod(total_time,3600.0)/60.0 )
  secs=int(mod(mod(total_time,3600.0),60.0))
  write(6,'(a14,i2.2,a1,i2.2,a1,i2.2)')"Elapsed time:   ",hrs,":",mins,":",secs
  write(6,*) "Elapsed time (s) =", total_time
  write(6,*) "end of program"
  write(6,*) 'numradars =',nradars
  write(6,*) 'numgates =',numgates
  write(6,*) 'numelevs =',nelv
  write(6,*) 'numazims =',azimuths
  write(6,*) 'numtimes =',ntime
  write(6,*) 'ndata    =',ndata
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
