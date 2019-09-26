#!/bin/ksh
export ndate=/home/Donald.E.Lippi/bin/ndate

pdy_m0=20190904
cyc_m0=00


valtime_m3=`${ndate} -3 $valtime_m0`
valtime_m2=`${ndate} -2 $valtime_m0`
valtime_m1=`${ndate} -1 $valtime_m0`
valtime_m0=`${ndate} +0 $valtime_m0`
valtime_p1=`${ndate} +1 $valtime_m0`
valtime_p2=`${ndate} +2 $valtime_m0`
valtime_p3=`${ndate} +3 $valtime_m0`


pdy_m3=`echo ${valtime_m3} | cut -c 1-8`
pdy_m2=`echo ${valtime_m2} | cut -c 1-8`
pdy_m1=`echo ${valtime_m1} | cut -c 1-8`
pdy_m0=`echo ${valtime_m0} | cut -c 1-8`
pdy_p1=`echo ${valtime_p1} | cut -c 1-8`
pdy_p2=`echo ${valtime_p2} | cut -c 1-8`
pdy_p3=`echo ${valtime_p3} | cut -c 1-8`

cyc_m3=`echo ${valtime_m3} | cut -c 9-10`
cyc_m2=`echo ${valtime_m2} | cut -c 9-10`
cyc_m1=`echo ${valtime_m1} | cut -c 9-10`
cyc_m0=`echo ${valtime_m0} | cut -c 9-10`
cyc_p1=`echo ${valtime_p1} | cut -c 9-10`
cyc_p2=`echo ${valtime_p2} | cut -c 9-10`
cyc_p3=`echo ${valtime_p3} | cut -c 9-10`

PATH1="/gpfs/hps3/emc/meso/save/Donald.E.Lippi/rw_tools/nexrad_1p00_250/simbufr/"
nexrad_L2RWBUFRm3="nexrad_${pdy_m3}_fv3.t${cyc_m3}z_drw.bufr"
nexrad_L2RWBUFRm2="nexrad_${pdy_m2}_fv3.t${cyc_m2}z_drw.bufr"
nexrad_L2RWBUFRm1="nexrad_${pdy_m1}_fv3.t${cyc_m1}z_drw.bufr"
nexrad_L2RWBUFRm0="nexrad_${pdy_m0}_fv3.t${cyc_m0}z_drw.bufr"
nexrad_L2RWBUFRp1="nexrad_${pdy_p1}_fv3.t${cyc_p1}z_drw.bufr"
nexrad_L2RWBUFRp2="nexrad_${pdy_p2}_fv3.t${cyc_p2}z_drw.bufr"
nexrad_L2RWBUFRp3="nexrad_${pdy_p3}_fv3.t${cyc_p3}z_drw.bufr"
sha1sum $nexrad_L2RWBUFRm3 $nexrad_L2RWBUFRm2 $nexrad_L2RWBUFRm1 $nexrad_L2RWBUFRm0
sha1sum $nexrad_L2RWBUFRp1 $nexrad_L2RWBUFRp2 $nexrad_L2RWBUFRp3


echo $valtime_m3
echo $valtime_m2
echo $valtime_m1
echo $valtime_m0
echo $valtime_p1
echo $valtime_p2
echo $valtime_p3
