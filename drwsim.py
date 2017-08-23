#!/usr/bin/python
import drawRadarCoverage
#import math
#from math import atan2,degrees
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.basemap import Basemap, cm
import ncepbufr
import ncepy
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
    factor=int(10)  # factor for upscaling the native resolution of model res.  #
    resample=True   # if False, does no resample to get upscaled drws.          # 
#    tilts=np.array(np.arange(0.5,20,0.5)) # 0.5 to 19.5 by 0.5 increments     #
    tilts=np.array([0.5])
                    # radar elevation angles to process.                        #
    mindbz=0.       # only compute drws where >= mindbz.                        #
    staid='KGRK'    # station id for debugging.                                 #
    makefigs=True   # would you like to make figures? Helpful in debugging.     #
    grib=False      # is the model data in grib2 format?                        #
    netcdf=True     # is the model data in netcdf format?                       #
    gribfile  ='/home/Donald.E.Lippi/obssim/namrr.t00z.conusnest.hiresf18.tm00.grib2'
    netcdfpath='/scratch4/NCEPDEV/meso/save/Donald.E.Lippi/data'
##### END OF USER DEFINED SETTINGS ##############################################
    tic=time.clock()
    if(grib): # if the model data is in grib2 format...
       grbs = pygrib.open(gribfile)
       #print_grib_inv(grbs)
       # Get the needed fields.
       ugrbs  =grbs.select(name='U component of wind',typeOfLevel='isobaricInhPa')
       vgrbs  =grbs.select(name='V component of wind',typeOfLevel='isobaricInhPa')
       wgrbs  =grbs.select(name='Vertical velocity'  ,typeOfLevel='isobaricInhPa')
       # NEED TO GET TEMPERATURE AND SEA LEVEL PRESSURE FOR USE IN HEIGHT CALCULATION
       dbzgrbs=grbs.select(name='Derived radar reflectivity',typeOfLevel='heightAboveGround')
       dbz=dbzgrbs[0].values
       pcoord=np.arange(10000,102500,2500) # probably only uses 1000mb - 
       gridspacing=dbzgrbs[0].DxInMetres/1000.
       lats,lons=dbzgrbs[0].latlons()
       yyyymmdd=ugrbs[0].dataDate; hh=ugrbs[0]['forecastTime']; fhr=hh
       yyyy=ugrbs[0].year; mm=ugrbs[0].month; dd=ugrbs[0].day; cyc=ugrbs[0].hour
 
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
       #print_nc_variables(fnd3d,fnd2d,fng); exit()
       lons  = fng.variables['grid_lont'][:,:] # these are much different than grib
       lats  = fng.variables['grid_latt'][:,:]
       # Grab the cycledate
       times = fnd2d.variables['time'][:]
       cdftime = utime(getattr(fnd2d.variables['time'],'units'))
       cycledate=roundTime(cdftime.num2date(times[0]),roundTo=60.*60.)
       yyyy=cycledate.year; mm=cycledate.month; dd=cycledate.day; cyc=cycledate.hour
       fhr=12; hh=fhr/3 #output is every 3 hours
       yyyymmdd=int(str(yyyy)+str(mm).zfill(2)+str(dd).zfill(2))
       # Get the needed fields.
       # These fields are of shape (time, pfull, grid_yt, grid_xt) = (25, 63, 1440, 1728) for a 72hr fcst.
       unc    = fnd3d.variables['ucomp']
       vnc    = fnd3d.variables['vcomp']
       wnc    = fnd3d.variables['w']
       pcoord = fnd3d.variables['pfull'] # 63 levels
       #dbznc = fnd2d.variables[''] # Don't see this in the fv3 output yet.
       # Let's approximate dbz for now.
       pratesfc=fnd2d.variables['PRATEsfc'] # surface precipitation rate
       # https://vlab.ncep.noaa.gov/web/wdtd/-/surface-precipitation-rate-spr-?selectedFolder=668041
       dbz=pratesfc[hh,:,:]  # Z = 300*R**1.4
       gridspacing=3. # I don't know how to get this from the netcdf file yet...

    # Read in the microsoft excel comma separated value (csv) file which contains radar database. 
    df=pd.read_csv('small_list_for_radar_sim_dev.csv')
    #print("number of messages = "+str(grbs.messages))

    #radius=100; radius=radius+2.*gridspacing; dx=int(np.ceil(radius/gridspacing)); dy=dx
    radius=100; dx=int(np.ceil(radius/gridspacing)); dy=dx
    if(resample):
       factor=factor
       thetas=360                #360 azimuths
       gate=250                  #drws are typically in 250 m gates.
       ngates=radius*1000/gate   #100km/250m = 400 gates.
    else:
       factor=1                  #if we're not resampling, reset factor to 1.
       thetas=360                #360 azimuths
       gate=gridspacing*1000     #grid spacing in meters (3000-m)
       ngates=np.ceil(radius*1000/gate)       

    toc=time.clock()
    print("step 1 took "+str(toc-tic)+" seconds")

    nummessages=0
    for rid in xrange(len(df.ID)):
       maxnummessages=1 # equal to the number of radars being processed.
       #if(strip(df.ID[rid])=='KAMA' or strip(df.ID[rid])==staid):
       if(strip(df.ID[rid])==staid):
          nummessages=nummessages+1
          tic=time.clock()
          # get coord of radar
          x,y=ncepy.find_nearest_ij(lats,df.Lat[rid],lons,df.Lon[rid])
          # create subset for faster processing. 
          dbzsub,lat,lon=dbz[x-dx:x+dx,y-dy:y+dy],lats[x-dx:x+dx,y-dy:y+dy],lons[x-dx:x+dx,y-dy:y+dy]
          dbzsub,lat,lon=upscale(dbzsub,factor),upscale(lat,factor),upscale(lon,factor)
          drw=np.zeros(shape=(len(tilts),len(dbzsub),len(dbzsub))); drw.fill(-999)#; drw=np.ma.masked_array(drw,dbzsub.mask)
          drwpol=np.zeros(shape=(len(tilts),thetas,ngates)); drwpol.fill(-999.)
          twodxfactor=2*dx*factor
          toc=time.clock()
          print("step 2 took "+str(toc-tic)+" seconds for #"+str(rid)+' '+str(strip(df.ID[rid])))
          print("upscale factor = "+str(factor)+" and subset dimensions of "+str(2*dx*factor)+" by "+str(2*dx*factor))
          tic=time.clock()
          usub=[]; vsub=[]; wsub=[]
          if(grib):
             for level in xrange(len(ugrbs)): # putting this before the a,b loop significanly sped things up!
                usub.append(upscale(ugrbs[level].values[x-dx:x+dx,y-dy:y+dy],factor))
                vsub.append(upscale(vgrbs[level].values[x-dx:x+dx,y-dy:y+dy],factor))
                wsub.append(upscale(wgrbs[level].values[x-dx:x+dx,y-dy:y+dy],factor))
          if(netcdf):
             for level in xrange(len(unc[hh,:,0,0])): 
                usub.append(upscale(unc[hh,level,x-dx:x+dx,y-dy:y+dy],factor))
                vsub.append(upscale(vnc[hh,level,x-dx:x+dx,y-dy:y+dy],factor))
                wsub.append(upscale(wnc[hh,level,x-dx:x+dx,y-dy:y+dy],factor))
          toc=time.clock()
          print("step 3 took "+str(toc-tic)+" seconds for #"+str(rid)+' '+str(strip(df.ID[rid])))
          alpha=0
          for tilt in tilts: # Loop for simulating drws.
             widgetstring=str(strip(df.ID[rid]))+' '+str(tilt)
             widgets=[widgetstring+': [',Percentage(),' ',Timer(),'] ', Bar(),' (', ETA(), ') ']
             bar = ProgressBar(widgets=widgets)
             for a in bar(xrange(2*dx*factor)):
                for b in xrange(2*dy*factor):
                   if(dbzsub[b,a] >= mindbz): # if dbz exists.
                      # check if obs are within range (100-km); get distance from radar to compute height.
                      dist,iswithinrange=check_iswithinrange(df.Lat[rid],df.Lon[rid],lat[b,a],lon[b,a],ngates)
                      dist250=int(np.round(dist/gate)) #dist250=int(np.round(dist/250.))
                      # beam width at dist.
                      beamwidth=1.
                      if(iswithinrange): # check if horizontal dist is within range.
                        # compute compass bearing
                        radarlatlon=(df.Lat[rid],df.Lon[rid]); obslatlon=(lat[b,a],lon[b,a]) # must be tuple
                        azm=int(np.round(calculate_initial_compass_bearing(radarlatlon,obslatlon))) # round, integer
                        if(azm>=360): azm=azm-360
#                        tilt=tilts[0]
                        # four thirds height and corrected tilt?
                        tiltrad=np.radians(tilt)
                        height=fourthirdsheight(dist,df.Height[rid],tiltrad)
                        # find nearest pressure v-coord and slice u,v,w
                        nearest,idx=height2nearestpressure(height,pcoord)      
                        u=usub[idx]; v=vsub[idx]; w=wsub[idx]
                        azmrad=np.radians(azm)
                        costilt=np.cos(tiltrad)
                        cosazm=np.cos(azmrad)
                        sintilt=np.sin(tiltrad)
                        sinazm=np.sin(azmrad)
                        if(makefigs and tilt==tilts[0]):
                           drw[alpha,b,a]=u[b,a]*costilt*cosazm+v[b,a]*costilt*sinazm+w[b,a]*sintilt
                        drwpol[alpha,azm,dist250]=u[b,a]*costilt*cosazm+v[b,a]*costilt*sinazm+w[b,a]*sintilt
                      else:
                        drw[alpha,b,a]=-999.
             alpha=alpha+1
          # Mnemonics for getting data from the bufr file.
          print("Writing bufr file for "+str(strip(df.ID[rid])))
          hdstr='SSTN CLON CLAT SELV ANEL YEAR MNTH DAYS HOUR MINU QCRW ANAZ'
          obstr='DIST125M DMVR DVSW'                     #NL2RW--level 2 radial wind.
          l2rwdf=pd.read_csv('l2rwbufr.table.csv')
          date=yyyymmdd
          cyc=00
          fhr=fhr
          message_type=strip(l2rwdf.MNEMONIC[fhr+cyc])
          del_anel=0.25
          # create a message with radial velocity data.
          idate=int(str(yyyymmdd)+str(fhr).zfill(2))
          if(nummessages==1): # we want to put all messages from all radars and all tilts in one bufr file.
             bufr = ncepbufr.open(str(idate)+'_fv3.t'+str(fhr)+'z.drw.bufr','w',table='l2rwbufr.table') # bufr file for reading.
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
          alpha=0
          for tilt in tilts: # Loop for writing simulated drws to bufr file.
             hdr[4] = tilt       #ANEL - ANTENNA ELEVATION ANGLE
             bufr.open_message(subset,idate) # open a new message for each tilt.
             # set obs for radial velocity
             for phi in xrange(thetas):
                obs = bufr.missing_value*np.ones((len(obstr.split()),ngates),np.float)
                hdr[11]= np.float(phi)  #ANAZ - ANTENNA AZIMUTH ANGLE
                i=0                   # counter
                r=0.                  # meters
                rend=radius*1000      # meters - 400 gates for 100km
                while(r < rend):
                    obs[0,i] = r/gate            #DISTANCE (FROM ANTENNA TO GATE CENTER) IN UNITS OF 125M 
                    obs[1,i] = drwpol[alpha,phi,r/gate]  #DOPPLER MEAN RADIAL VELOCITY 
                    obs[2,i] = 1.                #DOPPLER VELOCITY SPECTRAL WIDTH
                    r=r+gate; i=i+1
                # encode radial velocity
                bufr.write_subset(hdr,hdstr)
                bufr.write_subset(obs,obstr,end=True) # end subset
             # close bufr message and bufr file.
             bufr.close_message()
             if(tilt==tilts[-1] and nummessages==maxnummessages): # close if no more tilts to process.
                bufr.close(); bufrisclosed=True
                print("Done writing bufr file.")# for "+str(strip(df.ID[rid])))
             alpha=alpha+1

          fignamestr=str(idate)+'_'+str(strip(df.ID[rid]))+'_'+str(tilt)
          if(makefigs and len(tilts)==1):
             # Set up basemap
             print("Making figures for "+str(strip(df.ID[rid])))
             dom='CONUS'
             proj='lcc'
             llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat,res=ncepy.corners_res(dom,proj=proj)
             fig = plt.figure(figsize=(12,8))
             ax = fig.add_axes([0.1,0.1,0.8,0.8])
             lat_1=25.0  # True latitude for the LCC projection
             lon_0=-95.0 # Reference longitude for the LCC projection
             m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
                      rsphere=(6378137.00,6356752.3142),\
                      resolution=res,projection=proj,\
                      lat_1=lat_1,lon_0=lon_0,ax=ax)
             m.drawcoastlines(linewidth=1.25)
             m.drawstates(linewidth=1.25)
             m.drawcountries(linewidth=1.25)
             m.drawcounties(linewidth=0.2)
             # DRW on basemap 
             draw_radar_coverage(df,rid,m)
             cmap=drw_colormap()
             clevs=np.arange(-40,41,1)
             cs = m.contourf(lon,lat,drw[tilt,:,:],clevs,cmap=cmap,latlon=True,extend='both')
             plt.title('drw')
             plt.savefig('./drw'+fignamestr+'.png',bbox_inches='tight')
             print("Fig 1 (drw) done")
             # DBZ on basemap 
             draw_radar_coverage(df,rid,m)
             cmap=ncepy.mrms_radarmap()
             clevs = [5.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,65.,70.,75.]
             cs = m.contourf(lon,lat,dbzsub,clevs,cmap=cmap,latlon=True,extend='both')
             #cs = m.contourf(lons,lats,dbz,clevs,cmap=cmap,latlon=True,extend='both')
             plt.title('dbz')
             plt.savefig('./dbz'+fignamestr+'.png',bbox_inches='tight')
             print("Fig 2 (dbz) done")
             # Polar plot 
             fig = plt.figure(figsize=(8,8)) # 8" x 8" seems plenty large.
             ax = fig.add_subplot(111,polar=True) # we would like it to be a polar plot.
             ax.set_theta_zero_location("N") # set theta zero location to point North.
             ax.set_theta_direction(-1) # set theta to increase clockwise (-1).
             theta,r = np.meshgrid(np.arange(0,thetas),np.arange(1,ngates+1)) # create meshgrid
             theta=np.radians(theta) # convert theta from degrees to radians.
             cmap=drw_colormap()
             mesh = ax.pcolormesh(theta,r,drwpol[tilt,:,:].T,shading='flat',cmap=cmap,vmin=-40,vmax=40) # plot the data.
             ax.grid(True)
             cbar = fig.colorbar(mesh,shrink=0.85,pad=0.10,ax=ax) # add a colorbar.
             plt.title('drw')
             plt.savefig('./drw_polar'+fignamestr+'.png',bbox_inches='tight')
             print("Fig 3 (drw polar) done")
             print("Done making figures for "+str(strip(df.ID[rid])))
          if(makefigs and len(tilts)!=1):
             exit("Try changing len of tilts to 1 or running with makefigs=False")
    if(bufrisopen):
       bufr.close()
       print("bufr file is now closed")
    
###############################################################################
# Below are functions used by the main program:                               #
#     calculate_initial_compass_bearing                                       #
#     check_iswithinrange                                                     #
#     draw_radar_coverage                                                     #
#     drw_colormap                                                            #
#     find_nearest                                                            #
#     fourthirdsheight                                                        #
#     height2nearestpressure                                                  #
#     make_colormap                                                           #
#     print_grib_inv                                                          #
#     roundTime                                                               #
#     strip                                                                   #
#     upscale                                                                 # 
###############################################################################

    #reference: https://gist.github.com/jeromer/2005586
def calculate_initial_compass_bearing(pointA, pointB):
    if (type(pointA) != tuple) or (type(pointB) != tuple):
        raise TypeError("Only tuples are supported as arguments")
    lat1 = np.radians(pointA[0])
    lat2 = np.radians(pointB[0])
    diffLong = np.radians(pointB[1] - pointA[1])
    x = np.sin(diffLong) * np.cos(lat2)
    y = np.cos(lat1) * np.sin(lat2) - (np.sin(lat1)
            * np.cos(lat2) * np.cos(diffLong))
    initial_bearing = np.arctan2(x, y)
    # Now we have the initial bearing but math.atan2 return values
    # from -180 to + 180 which is not what we want for a compass bearing
    # The solution is to normalize the initial bearing as shown below
    initial_bearing = np.degrees(initial_bearing)
    compass_bearing = (initial_bearing + 360) % 360
    return compass_bearing

def check_iswithinrange(lat1,lon1,lat2,lon2,ngates):
    d = ncepy.gc_dist(lat1,lon1,lat2,lon2)
    gate=100000./ngates
    dist=np.round(d/gate)*gate #round to nearest 250 for 250m gates.
    if(dist < 100000.):
       return(dist,True)
    else:
       return(dist,False)

def draw_radar_coverage(df,i,m):
    color='black'; radius=100
    glon1 = df.Lon[i]
    glat1 = df.Lat[i]
    X = []
    Y = []
    for azimuth in xrange(0,360):
       glon2, glat2, baz = drawRadarCoverage.shoot(glon1,glat1,azimuth,radius)
       X.append(glon2)
       Y.append(glat2)
    X.append(X[0])
    Y.append(Y[0])
    X,Y = m(X,Y)
    plt.plot(X,Y,color)

def drw_colormap():
    c = mcolors.ColorConverter().to_rgb
    cmap = make_colormap(
           [c('deepskyblue'),c('navy')    ,0.20, # light blue to dark blue
            c('#02ff02')    ,c('#003500') ,0.47, # bright green to dark green
            c('#809e80')    ,c('white')   ,0.50, # gray with green tint to white
            c('white')      ,c('#9e8080') ,0.53, # white to gray with red tint
            c('#350000')    ,c('#ff0000') ,0.80, # dark red to bright red
            c('salmon')     ,c('yellow')])       # salmon to yellow
    cmap.set_under('#999999')
    cmap.set_over('purple')
    return(cmap)

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx],idx

def fourthirdsheight(thisrange,stahgt,thistiltr):
    #4/3rds height
    if(True): # I don't think this is what I want to do...
       fourthirds=4./3.
       half=0.5
       r2=2.0
       r8=8.0
       erad=6371000    #m
       aactual=(erad+stahgt)
       a43 =aactual*fourthirds
       celev0=np.cos(thistiltr)
       selev0=np.sin(thistiltr)    
       b   = thisrange*(thisrange+r2*aactual*selev0)
       c   = np.sqrt(aactual*aactual+b)
       ha  = b/(aactual+c)
       epsh=(thisrange**2-ha**2)/(r8*aactual)
       h=ha-epsh
       thishgt=stahgt+h
    return(thishgt)

def height2nearestpressure(h,pcoord):
    # Hypsometric Eq. h = ((p0/p)^(1/5.257)-1)*(T+273.15)/0.0065
    #                 p = p0/(  (h*0.0065)/(T+273.15) + 1  )**5.257
    psfc=pcoord[-1]
    T=15
    ptop = psfc/((((h*0.0065)/(T+273.15)) + 1)**5.257)
    idx = (np.abs(pcoord-ptop)).argmin()
    #hft=h*3.28084 # convert height to feet.
    #Pa=(1013.25-hft/30.)*100 # in lower atm, 1hPa drop off per 30ft., and convert to Pa 
    #idx = (np.abs(pcoord-Pa)).argmin()
    return(pcoord[idx],idx)

def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)

def print_grib_inv(grbs):
    filename="./grib_inventory.txt"
    f=open(filename,"w")
    grbs.seek(0)
    for grb in grbs:
        #print(grb)
        f.write(str(grb)+"\n")
    f.close()

def print_nc_variables(fnd3d,fnd2d,fng):
    CGREEN='\033[32m'; CRED='\033[31m'; CORANGE='\033[33m'; CEND='\033[0m'
    print(CGREEN+"3D variables:\n"+str(fnd3d.variables.keys())+"\n")
    print(CRED+"2D variables:\n"+str(fnd2d.variables.keys())+"\n")
    print(CORANGE+"Grid specs  :\n"+str(fng.variables.keys())+CEND+"\n")

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

def upscale(arr,factor):
    out=scipy.ndimage.interpolation.zoom(input=arr, zoom=(factor), order = 3)
    return out

if __name__ == "__main__":
    tic=time.clock()
    cProfile.run('main()','main.profile')
    toc=time.clock()
    time=toc-tic
    hrs=int(time/3600)
    mins=int(time%3600/60)
    secs=int(time%3600%60)
    print("Total elapsed time: "+str(toc-tic)+" seconds.")
    print("Total elapsed time: "+str(hrs).zfill(2)+":"+str(mins).zfill(2)+":"+str(secs).zfill(2))
