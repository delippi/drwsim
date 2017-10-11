#!/usr/bin/python
from timeit import default_timer as timer
tic0=timer()
import sys
sys.path.append('./modules')
import drawRadarCoverage
from inspect import currentframe, getframeinfo
#import math
#from math import atan2,degrees
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.basemap import Basemap, cm
import ncepbufr
import ncepy
import math
import numpy as np
from numpy import exp, abs, angle
from netCDF4 import Dataset
from netcdftime import utime
import os
from datetime import datetime,timedelta
import pandas as pd
import pdb
import pygrib #https://github.com/jswhit/pygrib
import progressbar
from progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, \
                        FileTransferSpeed, FormatLabel, Percentage, \
                        ProgressBar, ReverseBar, RotatingMarker, \
                        SimpleProgress, Timer
import scipy.ndimage
import time
import cProfile
import gsimods
import interp
"""
Author     : Donald Lippi
Email      : donald.e.lippi@noaa.gov
Description: This program is intended for the use in an observing system simulation
             experiment (OSSE). This program simulates Doppler radial velocity
             observations given model output (either grib2 or netcdf) and a radar
             database provided in comma separated value format (csv - in MS excel).
             This program only creates DRWs where dBZ exists and writes the
             observations in the bufr format which is the format used in the GSI
             to process level 2 radial winds. This program can also upscale the
             resolution of the model data in order to create full resolution DRW
             observation data (360-degrees, 250-meter gates out to 100-km).
"""

fix='./fix/'
out='./out/'

#-Constants

izero=int(0); ione=int(1)
zero=0; half=0.50; one=1.00; two=2.00; three=3.00; four=4.00; eight=8.00
fourthirds=four/three; pihalf=np.pi*half
huge=999999
sign_pole=-one
erad=6378136.62 # value in GSI-constants - erad=6378137.  # value in NCEPy
rad_per_meter=1./erad
grav_equator = 9.7803253359
grav_polar   = 9.8321849378
semi_minor_axis=6356.7523142e3
semi_major_axis=6378.1370e3
earth_omega     = 7.292115e-5
grav_constant   = 3.986004418e14
grav_ratio = (earth_omega*earth_omega * semi_major_axis*semi_major_axis * semi_minor_axis) / grav_constant
somigliana=(semi_minor_axis/semi_major_axis) * (grav_polar/grav_equator) - one
eccentricity_linear = math.sqrt(semi_major_axis**2 - semi_minor_axis**2)
eccentricity = eccentricity_linear / semi_major_axis
flattening = (semi_major_axis-semi_minor_axis)/semi_major_axis
grav = 9.81



# NOTES FOR DEVELOPEMENT:
#========================================================================================
#1. Consider parallelizizng the code. You'll have to loop over 900+ radars
#   and if each radar takes about 48 seconds to process (only including one tilt
#   with few observations), that'll take nearly 12 hours.
#2. I don't see a dbz field in the fv3 netcdf file, so I can't use it as a
#   constraint for the location of drws.
#3. Do I need to rotate the winds in the fv3 file? I do for the nmmb.
#4. Improve vertical interpolation. Right now I'm getting some stratified looking winds.
#DONE 5. Add a loop over tilts. Process 0.5-19.5 for each call to radar. 
#6. Need to use dist along beam as dist from radar rather than horizontal dist.
#DONE 7. Put all bufr messages from all radars and all tilts in same bufr file.
#========================================================================================


def main():
######### USER DEFINED SETTINGS ################################################# 
    #tilts=np.array(np.arange(0.5,20,0.5)) # 0.5 to 19.5 by 0.5 increments      #
    tilts=np.array([0.5])                                                       #
                    # radar elevation angles to process.                        #
    mindbz=0.       # only compute drws where >= mindbz.                        #
    staid='KGRK'    # station id for debugging.                                 #
    makefigs=True   # would you like to make figures? Helpful in debugging.     #
    netcdf=True     # is the model data in netcdf format?                       #
    netcdfpath='/scratch4/NCEPDEV/meso/save/Donald.E.Lippi/data'                #
    maxobrange=100000                                                           #
    azimuths=360                    #360 azimuths                               #
    ithin=4                         #radar data is very hires 1=full res        #
    gatespc=250.*ithin              #drws are typically in 250 m gates.         #
    gates=int(maxobrange/gatespc)   #100km/250m = 400 gates.                    #
    minobrange=0                                                           #
##### END OF USER DEFINED SETTINGS ##############################################
    tic1=timer()

    print("maxobrange "+str(maxobrange))
    print("minobrange "+str(minobrange))
    print("ithin "     +str(ithin))
    print("gatespc "   +str(gatespc))
    print("gates "     +str(gates))
    print("azimuths "  +str(azimuths))
    print("tilts "     +str(tilts))


    # Read in the microsoft excel comma separated value (csv) file which contains radar database. 
    df=pd.read_csv(fix+'small_list_for_radar_sim_dev.csv')

    if(netcdf): # if the model data is in netcdf format...
       dir=netcdfpath
       nesteddata3d = os.path.join(dir,'nggps3d.nest02.nc')
       nesteddata2d = os.path.join(dir,'nggps2d.nest02.nc')
       nestedgrid   = os.path.join(dir,'grid_spec.nest02.nc')
       print('Reading {:s}'.format(nesteddata3d))
       fnd3d = Dataset(nesteddata3d,'r')
       print('Reading {:s}'.format(nesteddata2d))
       fnd2d = Dataset(nesteddata2d,'r')
       print('Reading {:s}'.format(nestedgrid))
       fng   = Dataset(nestedgrid,'r')
       lons  = fng.variables['grid_lont'][:,:]
       lats  = fng.variables['grid_latt'][:,:]
       #-compute ll of grid centers. This loop stays outside of tll2xy for efficiency.
       t1=timer()
       midlats=np.zeros(shape=(len(lats[:,0])-1,len(lats[0,:])-1))
       midlons=np.zeros(shape=(len(lons[:,0])-1,len(lons[0,:])-1))
       for j in xrange(len(midlats[0,:])-1):
           for i in xrange(len(midlats[:,0])-1):
               midlats[i,j]=(lats[i,j]+lats[i,j+1]+lats[i+1,j]+lats[i+1,j+1])/4
               midlons[i,j]=(lons[i,j]+lons[i,j+1]+lons[i+1,j]+lons[i+1,j+1])/4
       t2=timer()
       print("mid points took "+str(t2-t1)+" seconds")

       # Grab the cycledate
       times = fnd2d.variables['time'][:]
       cdftime = utime(getattr(fnd2d.variables['time'],'units'))
       cycledate=roundTime(cdftime.num2date(times[0]),roundTo=60.*60.)
       yyyy=cycledate.year; mm=cycledate.month; dd=cycledate.day; cyc=cycledate.hour
       fhr=12; hh=fhr/3 #output is every 3 hours
       yyyymmdd=int(str(yyyy)+str(mm).zfill(2)+str(dd).zfill(2))

       # Get the needed fields of shape (time,pfull,grid_yt,grid_xt)=(25,63,1440,1728) for a 72hr fcst.
       ges_u  = fnd3d.variables['ucomp'][hh,:,:,:]
       ges_v  = fnd3d.variables['vcomp'][hh,:,:,:]
       ges_w  = fnd3d.variables['w'][hh,:,:,:]
       ges_z  = 0
       geop_hgtl=fnd3d.variables['delz'][hh,:,:,:]
       dbznc  = fnd2d.variables['PRATEsfc'][hh,:,:] # Don't see dbz in the fv3 output yet.
       pcoord = fnd3d.variables['pfull'] # 63 levels
       gridspacing=3000. # I don't know how to get this from the netcdf file yet...


    toc=timer()
    print("step 1 took "+str(toc-tic1)+" seconds out of "+str(toc-tic0))

    count=0; irid=-1
    for rid in xrange(len(df.ID)):
       if(strip(df.ID[rid])==staid):
          count+=1; irid+=1
          # define an observation space of missing values.
          drwpol=np.zeros(shape=(len(tilts),azimuths,int(gates))); drwpol.fill(-999.)
          
          deltiltmax, deltiltmin, deldistmax, deldistmin, delazmmax=-huge,huge,-huge,huge,-huge

          stahgt=df.Height[rid]
          rlon0=np.radians(df.Lon[rid])
          rlat0=np.radians(df.Lat[rid])
          clat0=np.cos(rlat0)
          slat0=np.sin(rlat0)
          itilt=-1
          for thistilt in tilts: # Loop for simulating drws.
             itilt+=1
             thistiltr=np.radians(thistilt)
             celev0=np.cos(thistiltr)
             selev0=np.sin(thistiltr)

             #widgetstring=str(strip(df.ID[rid]))+' '+str(thistilt)
             #widgets=[widgetstring+': [',Percentage(),' ',Timer(),'] ', Bar(),' (', ETA(), ') ']
             #bar = ProgressBar(widgets=widgets)
             iazm=-1
             #for thisazimuth in bar(xrange(azimuths)): #thetas=360
             for thisazimuth in (xrange(azimuths)): #thetas=360
              iazm+=1
              if(thisazimuth == 60):
                thisazimuthr=np.radians(thisazimuth)
                igate=-1
                widgetstring=str(strip(df.ID[rid]))+' '+str(thistilt)+' '+str(thisazimuth)
                widgets=[widgetstring+': [',Percentage(),' ',Timer(),'] ', Bar(),' (', ETA(), ') ']
                bar = ProgressBar(widgets=widgets)
                for gate in bar(xrange(int(gates))):         #gates=400 or 100km
                 igate+=1
                 if(gate >= 20 and gate <= 22):
                   inside=False
                   if(gate*gatespc >= minobrange and gate*gatespc <= maxobrange): inside=True
                   if(inside):
                      #--Find observation height using method from read_l2bufr_mod.f90 
                      thisrange=gate*gatespc
                      aactual=(erad+stahgt)
                      a43=aactual*fourthirds
                      b   = thisrange*(thisrange+two*aactual*selev0)
                      c   = np.sqrt(aactual*aactual+b)
                      ha  = b/(aactual+c)
                      epsh=(thisrange**2-ha**2)/(eight*aactual)
                      h=ha-epsh
                      thishgt=stahgt+h


                      dpres=thishgt #store the absolute ob height (m) in dpres.                   


                      #--Find observation location using method from read_l2bufr_mod.f90

                      #-Get corrected tilt angle
                      celev=celev0
                      selev=selev0
                      celev=a43*celev0/(a43+h)
                      selev=(thisrange*thisrange+h*h+two*a43*h)/(two*thisrange*(a43+h))
 
                      gamma=half*thisrange*(celev0+celev) 

                      #-Get earth lat lon of ob
                      rlonloc=rad_per_meter*gamma*np.cos(thisazimuthr)
                      rlatloc=rad_per_meter*gamma*np.sin(thisazimuthr)

                      rlonglob,rlatglob=gsimods.invtllv(rlonloc,rlatloc,rlon0,clat0,slat0)

                      #-Find grid relative location of the radar
                      radar_lat=df.Lat[rid]
                      radar_lon=df.Lon[rid]
                      if(radar_lon>=360): radar_lon=radar_lon-360
                      if(radar_lon<zero): radar_lon=radar_lon+360
                      radar_x,radar_y=tll2xy(lats,lons,midlats,midlons,radar_lat,radar_lon)

                      #-Find grid relative location of the ob.
                      thislat=np.degrees(rlatglob)
                      thislon=np.degrees(rlonglob)
                      if(thislon>=360): thislon=thislon-360
                      if(thislon<zero): thislon=thislon+360
                      dlat,dlon=tll2xy(lats,lons,midlats,midlons,thislat,thislon)

                      #zsges=tintrp2a_single_level(gues_z,dlat,dlon)
                      zsges=100.
                      dpres=dpres-zsges
                      if(dpres<zero): exit("obs is below surface...")
                      nsig=63
                      hges=tintrp2a(geop_hgtl,dlat,dlon,nsig) 

                      #Convert geopotential height at layer midpoints to geometric height using
                      #equations (17, 20, 23) in MJ Mahoney's note "A discussion of various
                      #measures of altitude" (2001).  Available on the web at
                      #http://mtp.jpl.nasa.gov/notes/altitude/altitude.html
                      #termg, termr, termrg, zges  = equation 17, 21, 1st term denom 23, 23
                      thislat=np.radians(thislat); thislon=np.radians(thislon)
                      sin2 = np.sin(thislat)**2
                      termg = grav_equator*((one+somigliana*sin2)/np.sqrt(one-eccentricity**2*sin2))
                      termr = semi_major_axis / (one + flattening + grav_ratio - two*flattening*sin2)
                      termrg= (termg/grav)*termr
                      zges=np.zeros(shape=(nsig))
                      for n in range(nsig):
                          zges[n] = (termr*hges[n]) / (termrg-hges[n]) # eq (23)

                      #Convert observation height (in dpres) from meters to grid relative
                      #units.  Save the observation height in zob for later use.
                      zob = dpres
                      dpres=grdcrd(dpres,1,zges,nsig,1) 
                      #Interpolate guess u and v to observation location 
                      ugesin=interp.interp_util.tintrp3(ges_u,dlat,dlon,dpres)
                      vgesin=interp.interp_util.tintrp3(ges_v,dlat,dlon,dpres)
                      wgesin=interp.interp_util.tintrp3(ges_w,dlat,dlon,dpres)
                      #Convert guess u,v wind components to radial value consident with obs
                      cosazm  = np.cos(thisazimuthr) 
                      sinazm  = np.sin(thisazimuthr) 
                      costilt = np.cos(thistiltr) 
                      sintilt = np.sin(thistiltr) 
                      cosazm_costilt=cosazm*costilt
                      sinazm_costilt=sinazm*costilt
                      #-------------WIND FORWARD MODEL-----------------------------------------!
                      drwpol[itilt,iazm,igate] = ugesin*cosazm_costilt + vgesin*sinazm_costilt + wgesin*sintilt

                      if(dpres < zero or dpres > nsig): exit("something went wrong... ")
                      
          #------------BUFFERIZE---------------------------------------------------!
          # Mnemonics for getting data from the bufr file.
          print("Writing bufr file for "+str(strip(df.ID[rid])))
          hdstr='SSTN CLON CLAT SELV ANEL YEAR MNTH DAYS HOUR MINU QCRW ANAZ'
          obstr='DIST125M DMVR DVSW'                     #NL2RW--level 2 radial wind.
          l2rwdf=pd.read_csv(fix+'l2rwbufr.table.csv')
          date=yyyymmdd
          cyc=00
          fhr=fhr
          message_type=strip(l2rwdf.MNEMONIC[fhr+cyc])
          del_anel=0.25
          # create a message with radial velocity data.
          idate=int(str(yyyymmdd)+str(fhr).zfill(2))
          if(count==1): # Only open new file at the sign of first message.
             bufr = ncepbufr.open(out+str(idate)+'_fv3.t'+str(fhr)+'z.drw.bufr','w',table=fix+'l2rwbufr.table') # bufr file for reading.
             bufrisopen=True
          subset=strip(l2rwdf.MNEMONIC[fhr+cyc])
          # set header
          hdr    = bufr.missing_value*np.ones(len(hdstr.split()),np.float)
          hdr[0] = np.fromstring(strip(df.ID[rid])+'    ',dtype=np.float)[0] #SSTN - RADAR STATION IDENTIFIER (SHORT)
          hdr[1] = df.Lon[rid]      #CLON - LONGITUDE (COARSE ACCURACY)
          hdr[2] = df.Lat[rid]      #CLAT - LATITUDE (COARSE ACCURACY)
          hdr[3] = df.Height[rid]   #SELV - HEIGHT OF STATION
        # hdr[4] goes in tilt loop below.
          hdr[5] = yyyy       #YEAR - YEAR
          hdr[6] = mm         #MNTH - MONTH
          hdr[7] = dd         #DAYS - DAY
          hdr[8] = fhr        #HOUR - HOUR
          hdr[9] = 00         #MINU - MINUTE
          hdr[10]= 1          #QCRW - QUALITY MARK FOR WINDS ALONG RADIAL LINE
          itiltbufr=-1
          for tilt in tilts: # Loop for writing simulated drws to bufr file.
             itiltbufr+=1
             hdr[4] = tilt       #ANEL - ANTENNA ELEVATION ANGLE
             bufr.open_message(subset,idate) # open a new message for each tilt.
             # set obs for radial velocity
             for iazmbufr in xrange(azimuths):
                obs = bufr.missing_value*np.ones((len(obstr.split()),int(gates)),np.float)
                hdr[11]= np.float(iazmbufr)  #ANAZ - ANTENNA AZIMUTH ANGLE
                igatebufr=-1          # counter
                r=0.                  # meters
                rend=maxobrange       # meters - 400 gates for 100km
                while(r < rend):
                    igatebufr+=1
                    obs[0,igatebufr] = igatebufr   #DISTANCE (FROM ANTENNA TO GATE CENTER) IN UNITS OF 250M 
                    obs[1,igatebufr] = drwpol[itiltbufr,iazmbufr,igatebufr] #DOPPLER MEAN RADIAL VELOC
                    obs[2,igatebufr] = 1.                                   #DOPPLER VELOCITY SPECTRAL WIDTH
                    r+=gatespc 
                # encode radial velocity
                bufr.write_subset(hdr,hdstr)
                bufr.write_subset(obs,obstr,end=True) # end subset
             # close bufr message and bufr file.
             bufr.close_message()
             if(tilt==tilts[-1]): # close if no more tilts to process.
                bufr.close(); bufrisclosed=True
                print("Done writing bufr file.")# for "+str(strip(df.ID[rid])))



def grdcrd(d,nd,x,nx,flg):
#   Treat "normal" case in which nx>1
    d=np.array([d])
    for id in range(nd):
       if(nx > ione):
          if(flg==ione): 
             if(d[id]<=x[0]): ix=ione
             else:            ix=isrchf(nx-ione,x,d,flg)-ione 
             if(ix==nx-ione): ix=ix-ione
       elif(flg==-ione):
             if(d[id]>=x[0]): ix=ione
             else:            ix=isrchf(nx-ione,x,d,flg)-ione
       d[id]=float(ix)+(d-x[ix])/(x[ix+ione]-x[ix])
    return(d)

def isrchf(nx1,x,y,flg):
    if(flg==ione):
       for k in range(nx1):
           if(y<=x[k]): 
               isrchf=k
               break
    else:
       for k in range(nx1):
           if(y>=x[k]): 
               isrchf=k
               break 
    isrchf=nx1+ione
    if(nx1<=izero): isrchf=izero
    return(isrchf)

def tintrp3(fin,dxin,dyin,dzin): #tintrp3(f,gout,dxin,dyin,dzin)
#!   input argument list:
#!     f        - input interpolator
#!     dx,dy,dz - input x,y,z-coords of interpolation points (grid units)
#!     n        - number of interpolatees
#!   output argument list:
#!     g        - output interpolatees
    f=fin; nlevs,nlat,nlon=np.shape(f); dx=dxin; dy=dyin; dz=dzin; g=np.zeros(shape=(1)); n=1
    for i in range(n):
        ix1=int(dx)
        iy1=int(dy)
        iz =int(dz)
        delx=dx-ix1
        dely=dy-iy1
        delz=dz-iz
        ixp=ix1+1
        iyp=iy1+1
        izp=iz +1
        delxp=one-delx
        delyp=one-dely
        delzp=one-delz
        ix=ix1
        iy=iy1
        g[i] =(f[iz ,ix ,iy ]*delxp*delyp*delz  \
             + f[iz ,ixp,iy ]*delx *delyp*delz  \
             + f[iz ,ix ,iyp]*delxp*dely *delz  \
             + f[iz ,ixp,iyp]*delx *dely *delz  \
             + f[izp,ix ,iy ]*delxp*delyp*delzp \
             + f[izp,ixp,iy ]*delx *delyp*delzp \
             + f[izp,ix ,iyp]*delxp*dely *delzp \
             + f[izp,ixp,iyp]*delx *dely *delzp )
    return(g)

def tintrp2a(fin,dxin,dyin,nlevs):#hges=tintrp2a(geop_hgtl,dlat,dlon,nsig)
#   input argument list:
#     f        - input interpolator
#     dx,dy    - input x,y,z-coords of interpolation points (grid units)
#     obstime  - time to interpolate to
#     gridtime - grid guess times to interpolate from
#     n        - number of interpolatees
#     nlevs    - number of vertical levels over which to perform the 
#                2-d intrpolation 
#   output argument list:
#     g        - output interpolatees
#    print("Check tintrp2a for correctness...")
    f=fin; nlevs,nlat,nlon=np.shape(f); dx=dxin; dy=dyin; g=np.zeros(shape=(nlevs,1)); n=1
    istart=0; jstart=0; lon1=nlon
    for i in range(n):
        ix1=int(dx)
        iy1=int(dy)
        delx=dx-ix1
        dely=dy-iy1
        ixp=ix1+1
        iyp=iy1+1
        delxp=one-delx
        delyp=one-dely
        ix=ix1
        iy=iy1
        for k in range(nlevs):
            g[k,i] =(f[k,ix,iy]*delxp*delyp+f[k,ixp,iy]*delx*delyp \
                   + f[k,ix,iyp]*delxp*dely+f[k,ixp,iyp]*delx*dely) 
    return g

def tintrp2a_single_level(fin,dxin,dyin):
#    print("tintrp2a_single_level: Make sure you get gues_z from nggps2d at some point.")
    f=fin; nlat,nlon=np.shape(f); nlevs=1; dx=dxin; dy=dyin; g=np.zeros(shape=(nlat,nlon)); n=1
    ix1=int(dx)
    iy1=int(dy)
    delx=dx-ix1
    dely=dy-y1
    ixp=ix1+1
    iyp=iy1+1
    delxp=one-delx
    delyp=one-dely
    ix=ix1
    iy=iy1
    g =(f[ix,iy]*delxp*delyp+f[ixp,iy]*delx*delyp + f[ix,iyp]*delxp*dely+f[ixp,iyp]*delx*dely)
    return g


def tll2xy(lats,lons,midlats,midlons,lat,lon):
  diag=False
  gcc=ncepy.gc_dist(midlats,midlons,lat,lon)
  # index of grid centers closet to radar.
  # xc, yc give the bottom right index of grid box...
  xc,yc=np.unravel_index(gcc.argmin(), gcc.shape)
  gc=ncepy.gc_dist(lats,lons,lat,lon)
  x,y=np.unravel_index(gc.argmin(), gc.shape)
  x1,y1=xc+0,yc+0 #bottom right
  x2,y2=xc+0,yc+1 #bottom left    point of grid box that contains the lat/lon
  x3,y3=xc+1,yc+1 #top left
  x4,y4=xc+1,yc+0 #top right
  dlat =ncepy.gc_dist(lats[x1,y1],lons[x1,y1],lats[x4,y4],lons[x4,y4])
  dlon =ncepy.gc_dist(lats[x1,y1],lons[x1,y1],lats[x2,y2],lons[x2,y2])
  dlatr=ncepy.gc_dist(lats[x1,y1],lons[x1,y1],lat        ,lons[x4,y4])
  dlonr=ncepy.gc_dist(lats[x1,y1],lons[x1,y1],lats[x2,y2],lon        )
  idxx=x1+dlatr/dlat
  idxy=y1+dlonr/dlon
  if(diag):
     print(str(midlats[xc,yc])+", "+str(midlons[xc,yc]))
     print(str(lats[x1,y1])+", "+str(lons[x1,y1]))
     print(str(lats[x2,y2])+", "+str(lons[x2,y2]))
     print(str(lats[x3,y3])+", "+str(lons[x3,y3]))
     print(str(lats[x4,y4])+", "+str(lons[x4,y4]))
     print(str(lat)        +", "+str(lon))
  return(idxx,idxy)#,dlat,dlon)


#def invtllv(ALM,APH,TLMO,CTPH0,STPH0):
##   input arg list:
##     alm  -- input eart longitude
##     aph  -- input earth latitude
##     tlmo -- input earth longitude of rotated grid origin (radrees)
##     ctph0-- cos(earth lat of rotated grid origin)
##     stph0-- sin(earth lat of rotated grid origin) 
##
##   output arg list
##     tlm  -- rotated grid longitude
##     tph  -- rotated grid latitude
#
#    RELM=ALM
#    SRLM=np.sin(RELM)
#    CRLM=np.cos(RELM)
#    SPH=np.sin(APH)
#    CPH=np.cos(APH)
#    CC=CPH*CRLM
#    ANUM=CPH*SRLM
#    DENOM=CTPH0*CC-STPH0*SPH
#    TLM=TLMO+np.arctan2(ANUM,DENOM)
#    TPH=np.arcsin(CTPH0*SPH+STPH0*CC)
#    return(TLM,TPH)

# HELPER FUNCTIONS ##############################

def roundTime(dt=None, roundTo=60):
   """Round a datetime object to any time laps in seconds
   dt : datetime.datetime object, default now.
   roundTo : Closest number of seconds to round to, default 1 minute.
   Author: Thierry Husson 2012 - Use it as you want but don't blame me.
   """
   if dt == None : dt = datetime.datetime.now()
   seconds = (dt.replace(tzinfo=None) - dt.min).seconds
   rounding = (seconds+roundTo/2) // roundTo * roundTo
   return dt + timedelta(0,rounding-seconds,-dt.microsecond)

def strip(csv):
    return (" ".join(str(csv).split()))


if __name__ == "__main__":
    tic=timer()
    cProfile.run('main()',out+'main.profile')
    #main()
    toc=timer()
    time=toc-tic
    hrs=int(time/3600)
    mins=int(time%3600/60)
    secs=int(time%3600%60)
    print("Total elapsed time: "+str(toc-tic)+" seconds.")
    print("Total elapsed time: "+str(hrs).zfill(2)+":"+str(mins).zfill(2)+":"+str(secs).zfill(2))

