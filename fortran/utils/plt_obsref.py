#!/usr/bin/python

import pygrib
import numpy as np
import matplotlib
matplotlib.use('Agg')   #Necessary to generate figs when not running an Xserver (e.g. via PBS)
import ncepy
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import time
import sys 
import  definitions as d

if __name__ == '__main__':


  t1a = time.clock()
  print("Starting...")

#  grbs = pygrib.open('/scratch2/portfolios/NCEPDEV/meso/noscrub/Jacob.Carley/POWER/RETRO/nwpower_retro_fcst_AugCTL/namrr.2004080912/namrr.t12z.conusnest.hiresf06.tm00')  
  
  # Read forecast valid time grib file from command line
  valpdy=sys.argv[1]  
  valcyc=sys.argv[2]
  gribfile=sys.argv[3]
  domid=sys.argv[4]
  grbs=pygrib.open(gribfile)
  
# Uncomment this below to print the contents of the grib file to the screen
#  for grb in grbs:
#    print(grb) 

  # Get the lats and lons
  lats, lons = grbs[1].latlons()
  
  # Grib grid projection info   -  not in the obs dbz files
  #LOV=grbs[1]['LoVInDegrees']
  #true_lat=grbs[1]['Latin1InDegrees']

  #Get the date/time and forecast hour
  fhr=grbs[1]['stepRange'] # Forecast hour
  # Pad fhr with a 0
  if int(fhr) < 10:
    fhr='0'+fhr
  cyctime=grbs[1].dataTime #Cycle (e.g. 1200 UTC)
 
  if fhr==0: cyctime=cytime+'00' 
  #Pad with a zero and convert to a string
  if cyctime < 1000:
    grbtime='0'+repr(cyctime)
  else:
    grbtime=repr(cyctime) 
  
  date=grbs[1].dataDate    #PDY


  # Specify some plotting domains which have the regions pre-set in ncepy
#  domains=['CONUS','NW','NC','NE','SW','SC','SE','Great_Lakes']
  domains=['SC']
  proj='lcc'
  llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat,res=ncepy.corners_res('SC',proj=proj)
  #llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat=-99.38,29.25,-95.38,32.18
  # Start reading fields to be plotted


  # read refc dBZ
  dbz_msgs = grbs.select(name='Maximum/Composite radar reflectivity')[0] # This is actually REFC dBZ
  #Now get the values from this msg
  dbz = dbz_msgs.values
  

  t2a=time.clock()
  t3a=round(t2a-t1a, 3)
  print(repr(t3a)+" seconds to read all gribs msgs!")

  ###################################################
  #       START PLOTTING FOR EACH DOMAIN            #
  ###################################################

  for dom in domains:
    
    t1dom=time.clock()
    print('Working on '+dom)

    # create figure and axes instances
    fig = plt.figure(figsize=(11,11))
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    
    # create LCC basemap instance and set the dimensions
    #llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat,res=ncepy.corners_res(dom)    
    if(dom == 'SC1'):
       llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat=-101.0,29.0,-94.0,33.0
    if(dom == 'SC2'):
       llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat=-98.0,29.5,-95.0,31.5
    if(dom == 'SC3'):
       llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat=-99.38,29.25,-95.38,32.18
    if(dom == 'SC4'):
       llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat=-104.0,28.0,-92.0,35.0

    Lon0=262.5 #grbs[1]['LoVInDegrees']
    Lat0=38.5  #grbs[1]['LaDInDegrees']
    Lat1=38.5  #grbs[1]['Latin1InDegrees'] 
    Lat2=38.5  #grbs[1]['Latin2InDegrees']
    gribproj=grbs[1]['gridType']
    #rearth=grbs[1]['radius']
    rearth=6371229
    m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
   	      #rsphere=(6378137.00,6356752.3142),\
   	      rsphere=rearth,\
   	      resolution=res,projection='lcc',\
   	      #lat_1=25.0,lon_0=-95.0,ax=ax)
   	      #lat_1=30.0,lon_0=-100.0,ax=ax)
   	      lat_1=Lat1,lat_2=Lat2,lat_0=Lat0,lon_0=Lon0,ax=ax)

    #parallels = np.arange(-80.,90,10.)
    #parallels = np.arange(30.,90,2.)
    #meridians = np.arange(0.,360.,2.)
    parallels = np.arange(26.,90,4.)
    meridians = np.arange(0.,360.,4.)
    m.drawmapboundary(fill_color='#7777ff')
    m.fillcontinents(color='#ddaa66', lake_color='#7777ff', zorder = 0)
    m.drawcoastlines(linewidth=1.25)
    m.drawstates(linewidth=1.25)
    m.drawcountries(linewidth=1.25)
    m.drawparallels(parallels,labels=[1,0,0,1])
    m.drawmeridians(meridians,labels=[1,0,0,1])
    # DRAW 100-KM RADAR COVERAGE CIRLCES AROUND EACH RADAR
    #d.drawRadarCoverage(m)


    print('Working on REFC for '+dom) 
    
    if dom != 'CONUS':  
      # Draw the the counties if not CONUS
      # Note that drawing the counties can slow things down!
      m.drawcounties(linewidth=0.2, color='k')
      skip=25
    else:
      skip=45
    
    #  Map/figure has been set up here (bulk of the work), save axes instances for
    #     use again later   
    keep_ax_lst = ax.get_children()[:]


    #Now plot REFC dBZ
    clevs = [5.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,65.,70.,75.]
    mycmap = ncepy.mrms_radarmap()
    cs = m.contourf(lons,lats,dbz,clevs,cmap=mycmap,latlon=True,extend='max')
    cs.set_clim(5,75) 
    cbar = m.colorbar(cs,location='bottom',pad="5%",ticks=clevs)
    cbar.ax.tick_params(labelsize=8.5) 
    cbar.set_label('dBZ')
    plt.title(domid+' Column Max Reflectivity \n'+repr(date)+' '+grbtime+'Z')
   
    plt.savefig('./refc_'+dom+'_'+domid+'_'+repr(date)+grbtime+'.png',bbox_inches='tight')
    t2 = time.clock()
    #t3=round(t2-t1, 3)
    #print(repr(t3)+" seconds to plot refc for: "+dom)
      
    t3dom=round(t2-t1dom, 3)
    print(repr(t3dom)+" seconds to plot ALL for: "+dom)
    plt.clf()
    t3all=round(t2-t1a,3)
    print(repr(t3all)+" seconds to run everything!")

  
  
