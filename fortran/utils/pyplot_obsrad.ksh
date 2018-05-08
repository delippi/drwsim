#!/bin/ksh -l

export ndate=/home/Donald.E.Lippi/bin/ndate
exp=obsrefd    


# 00z 5/20/2013

PDY=20171015
cyc=00 #starting hour 
maxhr=01 #max number of hours to pull


hr=00
typeset -Z2 hr

while [ $hr -le $maxhr ]; do
echo "\n >>>Making plot for ${PDY}${cyc} ${hr} max hour $maxhr"
pyplotdir=${TROOT}/stmp4/${USER}/pyplot_work_${exp}.${PDY}_${hr}
mkdir -p $pyplotdir

figout=${TROOT}/stmp4/${USER}/pyplot_nest_${exp}.${PDY}
mkdir -p $figout
cd ${pyplotdir}


domid=Obs
valtime=`${ndate} +${hr} ${PDY}${cyc}`
valpdy=`echo ${valtime} | cut -c 1-8`
valcyc=`echo ${valtime} | cut -c 9-10`
valyr=`echo ${valtime} | cut -c 1-4`
valmon=`echo ${valtime} | cut -c 5-6`
if [ ! -s refd3d.t${valcyc}z.grb2f00 ]; then
  htar -xvf /NCEPPROD/hpssprod/runhistory/rh${valyr}/${valyr}${valmon}/${valpdy}/com_hourly_prod_radar.${valpdy}.save.tar ./refd3d.t${valcyc}z.grb2f00
fi

gbfile=refd3d.t${valcyc}z.grb2f00

if [ ! -d ${TROOT}/ptmp/$USER/output/pyplot ]; then
mkdir -p ${TROOT}/ptmp/$USER/output/pyplot
fi

export gbfile=${gbfile}
export valcyc=${valcyc}
export valpdy=${valpdy}
export domid=${domid}

cd ${pyplotdir}
cp /home/$USER/plotting/python/plt_obsrad/plt_obsref.py .
cp /home/Donald.E.Lippi/plotting/python/pynemsio/pyscripts/definitions.py .

python plt_obsref.py ${valpdy} ${valcyc} ${gbfile} ${domid}
rm -f definitions.py

mv *png ${figout}/
echo "Your plots are here:"
echo ${figout}
#rm -rf ${pyplotdir} 

(( hr = hr + 1 ))


done


exit
