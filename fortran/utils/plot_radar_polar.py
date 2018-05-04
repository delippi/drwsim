from __future__ import print_function
import ncepbufr
import matplotlib
matplotlib.use('Agg')
import ncepy, sys
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import progressbar
from progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, \
                        FileTransferSpeed, FormatLabel, Percentage, \
                        ProgressBar, ReverseBar, RotatingMarker, \
                        SimpleProgress, Timer
import time
import matplotlib.colors as mcolors
import pdb
import subprocess

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

def main():
######### USER DEFINED SETTINGS #################################################
   #OBS_FILE='../run/2017080712_fv3.t12z.drw.bufr'
   #date=2017080712      # date and time of the observations                    #
   #OBS_FILE='../run/2017101500_fhr06_fv3.t00z.drw.bufr'
   OBS_FILE='../run/2017101500_fv3.t00z.drw.bufr'
   date=2017101500      # date and time of the observations                    #
   STAID='KGRK'         # station id you want to plot.                         #
                        # tilt angle you want to plot.                         #
   anel_list=[0.5]#,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
   field='RADIAL WIND'  # you can plot 'RADIAL WIND' or 'REFLECTIVITY'         #
   with open('/scratch4/NCEPDEV/meso/save/Donald.E.Lippi/PhD-globalOSSE/obssim/fortran/run/namelist','r') as searchfile:
        for line in searchfile:
            if "ithin" in line:
                ithin=line
                ps1=subprocess.Popen(("echo",ithin), stdout=subprocess.PIPE)
                ithin=subprocess.check_output(("cut","-f","2","-d","="),stdin=ps1.stdout).strip("\n")
                ps2=subprocess.Popen(("echo",ithin), stdout=subprocess.PIPE)
                ithin=int(subprocess.check_output(("cut","-f","1","-d",","),stdin=ps2.stdout).strip("\n"))
                break
        for line in searchfile:
            if "azimuths" in line:
                azimuths=line
                ps1=subprocess.Popen(("echo",azimuths), stdout=subprocess.PIPE)
                azimuths=subprocess.check_output(("cut","-f","2","-d","="),stdin=ps1.stdout).strip("\n")
                ps2=subprocess.Popen(("echo",azimuths), stdout=subprocess.PIPE)
                azimuths=float(subprocess.check_output(("cut","-f","1","-d","."),stdin=ps2.stdout).strip("\n"))
                break
   gatespc=250.*ithin
##### END OF USER DEFINED SETTINGS ##############################################
   for anel0 in anel_list:
    tic = time.clock()
    # get message type from date.
    if(field=='RADIAL WIND'):
       message_type='NC00'+str(int(str(date)[-2:])+6010) #Ex: NC006022 if for 12z and NC006023 is for 13z.
    elif(field=='REFLECTIVITY'):
       message_type='NC00'+str(int(str(date)[-2:])+6040) #Ex: NC006052 if for 12z and NC006053 is for 13z.

    # Mnemonics for getting data from the bufr file.
    hdstr= 'SSTN CLON CLAT SELV ANEL YEAR MNTH DAYS HOUR MINU QCRW ANAZ' #PRFR' # MGPT'
    obstr= 'DIST125M DMVR DVSW' # PRFR'                     #NL2RW--level 2 radial wind.
    obstr2='STDM SUPLON SUPLAT HEIT RWND RWAZ RSTD' #RWSOB--radial wind super ob.

    #1. INITIALIZING SOME BASIC LISTS AND GETTING INPUT DATA.
    i=0; sids=[]; lons=[]; lats=[]; l2rw=[]; anel=[]; anaz=[]; dist125m=[]; ymdhm=[]; radii=[]; PRF=[]
    #OBS_FILE=sys.argv[1]; message_type=sys.argv[2]; date=sys.argv[3]; STAID=sys.argv[4]

    del_anel=0.25
     
    #2. READ PREPBUFR FILE.
    b='false' # used for breaking the loop.
    bufr = ncepbufr.open(OBS_FILE) # bufr file for reading.
    bufr.dump_table('../run/l2rwbufr.table') # dump table to file.
    while bufr.advance() == 0: # loop over messages.
        print(bufr.msg_counter, bufr.msg_type, bufr.msg_date)
        if(bufr.msg_type == message_type):
            while bufr.load_subset() == 0: # loop over subsets in message.
                hdr = bufr.read_subset(hdstr).squeeze() # parse hdstr='SSTN CLON ... etc.'
                station_id = hdr[0].tostring() # convert SSTN to string.
                station_id=station_id.strip()  # remove white space from SSTN.
                if(station_id == STAID): # comes from input. used for picking a single SSTN.
                    if(hdr[4] >= anel0-del_anel and hdr[4] <= anel0+del_anel): # read an elevation angle.
                        obs = bufr.read_subset(obstr).squeeze() # parse obstr='DIST125M DMVR DVSW'
                        ###print(bufr.msg_counter,'SSTN,CLON,CLAT,ANAL,ANAZ : ' \
                        ###      ,station_id,hdr[1],hdr[2],hdr[4],hdr[11])
                        i=i+1  #number of observations counter.
                        sids.append(station_id) # station ids
                        #lons.append(360+hdr[1]) # longitudes
                        #lats.append(hdr[2]) # latitudes
                        l2rw.append(obs[1]) # level 2 radial winds
                        anel.append(hdr[4]) # elevation angles
                        anaz.append(hdr[11]) # azimuthal angles
                        #PRF.append(obs[2]) #Pulse repitition frequency
                        #dist125m.append(obs[0]) # distances in units of 125 m
                        #radii.append(obs[0]*125) #distances in units of 1 m
                        #radii.append(obs[0]) #distances in units of 1 m
                        radii.append(obs[0]*gatespc) #distances in units of 1 m
                        ymdhm.append(int(hdr[9]))
                        #print(ymdhm)
                        if(anaz[-1] >= 360.0): # we would like to break the loop now.
                            b='true'
        if(b == 'true'): break # stop reading after all anaz's read.
    bufr.close()
    try:
       MM=str(ymdhm[-1])
    except IndexError: 
       exit("You have an empty bufr file...")
    print(np.min(anel),np.max(anel))
    print(np.min(ymdhm),np.max(ymdhm))
    toc = time.clock() # check how long reading the data in took.
    sec = str(toc-tic)
    print('time it took to run: '+sec+' seconds.')

    #3. FIND THE MAX/MIN VALUES OF THE DIST125M AND THETA ARRAYS.
    maxRadii=100000 # set max of max distances to about 100 mi.
    minRadii=0 # min of min distances.

    #4. INITIALIZE THE POLAR PLOT AS ALL MISSING DATA (-999).
    fig = plt.figure(figsize=(8,8)) # 8" x 8" seems plenty large.
    ax = fig.add_subplot(111,polar=True) # we would like it to be a polar plot.
    ax.set_theta_zero_location("N") # set theta zero location to point North.
    ax.set_theta_direction(-1) # set theta to increase clockwise (-1).
    theta,r = np.meshgrid(anaz,np.arange(minRadii,maxRadii,gatespc)) # create meshgrid
    theta=np.radians(theta) # convert theta from degrees to radians.
    rw = np.zeros(shape=(len(theta),len(r[0]))) # initialize the polar plot with all 0's.
    rw.fill(-999) # change all values to missing data (-999).
    #PRF_new=np.zeros(shape=(len(theta),len(r[0])))
    #PRF_new.fill(-999)
    r=r.T; rw=rw.T#; PRF_new.T # transpose r and rw for later manipulation.
    
    #5. POPULATE THE EMPTY RW ARRAY WITH ACTUAL VALUES.
    #for i in range(len(anaz)): # for every azimuth angle ...
    widgetstring=STAID+' '+str(anel0)
    widgets=[widgetstring+': [',Percentage(),' ',Timer(),'] ', Bar(),' (', ETA(), ') ']
    bar = ProgressBar(widgets=widgets)
    for i in bar(xrange(int(360./azimuths))): # for every azimuth angle ...
        #print(i,'/',len(anaz))
        for j in xrange(len(radii[i])): # ... loop over every observation distance from radar ...
            for k in xrange(len(r[i])): # ... and loop over an equally 125m spaced array ...
                if(radii[i][j] == r[i][k]): # ... if the observation dist = dist125m ...
                    rw[i][k]=l2rw[i][j] # ... assign the value of the observation to that index.
                    #PRF_new[i][k]=PRF[i][j]
                    break # speeds things up by about 50-60%.
    rw[np.isnan(rw)] = -999 # the masked values are converted to nan. set nan to -999.
    print("max value is: "+str(np.max(rw))) # check that it is not nan.

    #6. FINISH MAKING THE POLAR PLOT WITH THE FILLED IN VALUES.
    if(False):
       cmap = plt.cm.jet # use the jet colormap.
       cmap.set_under('white') # set the -999 values to white.
    else:
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
    mesh = ax.pcolormesh(theta,r.T,rw.T,shading='flat',cmap=cmap,\
                         vmin=-40,vmax=40) # plot the data.
    cbar = fig.colorbar(mesh,shrink=0.85,pad=0.10,ax=ax) # add a colorbar.
    cbar.set_label('$m/s$') # radial wind data is in units of meters per second.
    plt.title('Doppler Velocity \n Station ID: '+STAID\
              +'  Scan Angle: '+str(anel[0])\
              +'  Date: '+str(date)+MM.zfill(2),fontsize=15,y=1.12) # add a useful title.
    ax.grid(True)
    plt.show() # make the plot.
    plt.savefig('../figs/'+str(ithin*250)+'_'+STAID+'_'+str(anel[0])+'_'+str(date)+MM.zfill(2)+'.png'\
                ,bbox_inches='tight') # save figure.

    #7. CALCULATE SOME STATS
    toc = time.clock() # check to see how long it took to run.
    sec = str(toc-tic)
    print('time it took to run: '+sec+' seconds.')

if __name__ == "__main__":
    main()
