#!/bin/ksh
#PBS -N arch_obs 
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=1
#PBS -q service 
#PBS -A fv3-cpu
#PBS -o arch_obs.log 
#PBS -j oe

export ndate=/home/Donald.E.Lippi/bin/ndate

PDY=20180911
CYC=06
valtime=${PDY}${CYC}
cycles=1 #28 #how many cycles? 1=06z; 2=06z,12z; 3=06z,12z,18z; etc.
network="nexrad"


cycle=0
while [[ $cycle -le  $cycles ]]; do
   
   (( FH=6*$cycle ))
   (( cycle=cycle+1 ))
   valtime=`${ndate} +$FH ${PDY}${CYC}`

   valcyc=`echo $valtime | cut -c 9-10`

   HPSSDIR="/NCEPDEV/emc-meso/5year/Donald.E.Lippi/rw_FV3GFS/obs"
   simbufr="/scratch4/NCEPDEV/fv3-cam/save/Donald.E.Lippi/PhD-globalOSSE/obssim/fortran/run/simbufr"
   cd $simbufr
   #2018091106_fv3.t06z_drw.bufr
   obs=${network}_${valtime}_fv3.t${valcyc}z_drw.bufr
   echo ""
   echo ${obs} $FH

   htar -cvf ${HPSSDIR}/${obs}.tar ${obs} 
done
