import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def main():
    fill_between=True    # True for filling between 0.5 and 20 degree tilts
    both_directions=True # True for showing radar coverage to the left and right of station.
    proportional=False   # True for km x km plot, otherwise km x kft plot
    two=2.; r8=8.; rearth=6370000. #Constants
    deg2rad=3.14149265/180.00000000; m2ft=3.28084; m2km=0.001; ft2kft=0.001 # Conversion factors
    range_max=100000.
    if(both_directions): ranges=np.linspace(-1*range_max,range_max,501)
    else: ranges=np.linspace(0.,range_max,501)
    title_fontsize=19.; xy_axis_fontsize=16.
    if(proportional): # Makes an 8x8 plot with axes units km x km
       if(fill_between): tilts=np.array([0.5,20])
       else: tilts=np.array([0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20])
       ymax=140.; line_label_fontsize=9.; fig = plt.figure(figsize=(8,8))
       plt.gca().set_aspect('equal',adjustable='box')
    else: # Makes a 12x8 plot with axes km x kft
       if(fill_between): tilts=np.array([0.5,20])
       else: tilts=np.array([0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,12,14,16,18,20])
       ymax=70.; line_label_fontsize=12.; fig = plt.figure(figsize=(12,8))
    rads=tilts*deg2rad
    thishgt=np.zeros(shape=(len(tilts),len(ranges))) # Initialized values to 0
    i=-1
    for tilt in rads:
        selev0=np.sin(tilt)
        celev0=np.cos(tilt)
        i=i+1; j=0
        for thisrange in ranges:
            thisrange=np.abs(thisrange)
            b=thisrange*(thisrange+two*rearth*selev0)
            c=np.sqrt(rearth*rearth+b)
            ha=b/(rearth+c)
            epsh=(thisrange*thisrange-ha*ha)/(r8*rearth)
            h=ha-epsh
            thishgt[i,j]=h
            j=j+1
    for t in range(len(tilts)):
        if(proportional): x=ranges*m2km; y=thishgt[t,:]*m2km; yunits='km'; svfigstr='kmkm'
        else: x=ranges*m2km; y=thishgt[t,:]*m2ft*ft2kft; yunits='kft'; svfigstr='kmkft'
        label=str(tilts[t])+'$^\circ$'
        if(tilts[t]>=5.0):
           label=str(int(tilts[t]))+'$^\circ$'
        #plt.plot(x,y,color=c)
        plt.plot(x,y)
        if(y[-1] < ymax):
           plt.text(x[-1],y[-1],s=label,fontsize=line_label_fontsize,color='r',weight='bold') 
        else:
           z=np.argmax(y>=ymax)
           plt.text(x[z],ymax,s=label,fontsize=line_label_fontsize,color='r',weight='bold') 
    if(proportional and fill_between): plt.fill_between(x,thishgt[0,:]*m2km,thishgt[1,:]*m2km,facecolor='lightgray')
    elif(fill_between): plt.fill_between(x,thishgt[0,:]*m2ft*ft2kft,thishgt[1,:]*m2ft*ft2kft,facecolor='lightgray')
    plt.grid(b=True)
    plt.ylim(ymax=ymax,ymin=0)
    plt.xlabel('Range (km)',fontsize=xy_axis_fontsize,weight='bold')
    plt.ylabel('Height Above Radar Level ('+yunits+')',fontsize=xy_axis_fontsize,weight='bold')
    plt.title('Range vs. Height using the 4/3$^{rds}$ Rule \n',fontsize=title_fontsize,weight='bold')
   # plt.show()
    plt.savefig('./beam_height_'+svfigstr+'.png',bbox_inches='tight') 

if __name__ == "__main__":
   main()


