#!/bin/ksh
#PBS -N drwsim
#PBS -l walltime=0:30:00
#PBS -l nodes=1:ppn=20
#PBS -q debug
#PBS -A blender
#PBS -o log_drwsim.out
#PBS -j oe

src=../src

cd $src
make clean
make
cd -

exefile=drwsim.x
namelist=namelist
setup=setup.ksh
csvfile=radar_sim_dev.csv
ln -sf $src/${exefile} .
ln -sf ../fix/${csvfile} .
ln -sf ../fix/l2rwbufr.table.csv .
ln -sf ../fix/l2rwbufr.table .


FHMAX_GFS=72 #end hour
FHMIN_GFS=72   #start hour
FHOUT_GFS=3   #increment hour by
fcsthr=$FHMIN_GFS
typeset -Z3 fcsthr

CDUMP="gfs"
CYC="00"
PDY="20180911"
EXP="NATURE-2018091100-2018091800"

while [[ $fcsthr -le $FHMAX_GFS ]]; do
   #filename="gfs.t00z.atmf${fcsthr}.nemsio" #Get filename change
   filename="${CDUMP}.t${CYC}z.atmf${fcsthr}.nemsio" #Get filename change
   cp $namelist namelist.atmf${fcsthr}
   sed -i "s/@filename@/${filename}/g" namelist.atmf${fcsthr}
   sed -i "s/@CDUMP@/${CDUMP}/g"       namelist.atmf${fcsthr}
   sed -i "s/@PDY@/${PDY}/g"           namelist.atmf${fcsthr}
   sed -i "s/@CYC@/${CYC}/g"           namelist.atmf${fcsthr}
   sed -i "s/@EXP@/${EXP}/g"           namelist.atmf${fcsthr}

   cp $setup setup_atmf${fcsthr}.ksh 
   sed -i "s/@namelist.atmfxxx@/namelist.atmf${fcsthr}/g" setup_atmf${fcsthr}.ksh
   sed -i "s/@atmfxxx@/atmf${fcsthr}/g"                   setup_atmf${fcsthr}.ksh
   sed -i "s/@fcsthr@/${fcsthr}/g"                        setup_atmf${fcsthr}.ksh
   work=`pwd`
   sed -i "s#@work@#${work}#g"                            setup_atmf${fcsthr}.ksh
   ksh ./setup_atmf${fcsthr}.ksh
   (( fcsthr=fcsthr+${FHOUT_GFS} ))
done
