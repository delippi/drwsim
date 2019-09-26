#!/bin/ksh
#SBATCH -J getNatureFiles 
#SBATCH -t 1:00:00
#SBATCH -n 1
#SBATCH --partition=service 
#SBATCH -q batch
#SBATCH -A fv3-cpu
#SBATCH -o getNatureFiles096-119_0927.log 


instructions='''Before running this script you should have:
/NCEPDEV/emc-meso/5year/Donald.E.Lippi/rw_FV3GFS/NATURE-2018092300-2018100700

rh2018:
201809/  201810/  

rh2018/201809:
20180923/  20180924/  20180925/  20180926/  20180927/  20180928/  20180929/  20180930/  

rh2018/201809/20180923:
gfs.t00z.20180923.atm.nemsio.tar        gfs.t00z.20180923.pgrb2.0p25.tar        
gfs.t00z.20180923.atm.nemsio.tar.idx    gfs.t00z.20180923.pgrb2.0p25.tar.idx    
gfs.t00z.20180923.master.grb2.tar       gfs.t00z.20180923.sfc.nemsio.tar        
gfs.t00z.20180923.master.grb2.tar.idx   gfs.t00z.20180923.sfc.nemsio.tar.idx    

...

rh2018/201810:
20181001/  20181002/  20181003/  20181004/  20181005/  20181006/  20181007/

...

rh2018/201810/20181007:
gfs.t00z.20181007.atm.nemsio.tar        gfs.t00z.20181007.pgrb2.0p25.tar        
gfs.t00z.20181007.atm.nemsio.tar.idx    gfs.t00z.20181007.pgrb2.0p25.tar.idx    
gfs.t00z.20181007.master.grb2.tar       gfs.t00z.20181007.sfc.nemsio.tar        
gfs.t00z.20181007.master.grb2.tar.idx   gfs.t00z.20181007.sfc.nemsio.tar.idx 
'''

#This script downloads the known truth model output to generate simulated obs from.

valtime=2018092700 #20180923-20180930
pdy=`echo $valtime | cut -c 1-8`
cyc=`echo $valtime | cut -c 9-10`
year=`echo $valtime | cut -c 1-4`
yymm=`echo $valtime | cut -c 1-6`

cd /scratch1/NCEPDEV/stmp2/Donald.E.Lippi/fv3gfs_dl2rw
hsi "cd /NCEPDEV/emc-meso/5year/Donald.E.Lippi/rw_FV3GFS/NATURE-2018092300-2018100700/rh$year/$yymm/20180923/00/$pdy; cget -R gfs.t00z.${pdy}.atm.nemsio"

echo done
