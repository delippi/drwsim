#!/bin/ksh

instructions='''Before running this script you should run ./tools/getNatureFiles.ksh'''
echo $instructions
sleep 3

export ndate=/home/Donald.E.Lippi/bin/ndate

# Check for existence of sim directories which contain bufr, log, nml, and setup files.
# If they don't already exist, create them. Otherwise continues.
dirs="bufr log nml setup"
for dir in $dirs; do
   if [[ ! -e ./sim${dir} ]];then
      mkdir ./sim${dir}
   fi
done

# Compile the simulation code.
print -n "Compile drwsim? (y/n): ";read var; print ""
if [[ $var == 'y' ]]; then
   if [[ `hostname | cut -c 1` == "h" ]]; then; machine="hera"; fi
   if [[ `hostname | cut -c 1` == "t" ]]; then; machine="theia"; fi
   src=../src
   cd $src
   ln -sf Makefile.$machine Makefile
   make clean
   make
   cd -
   exefile=drwsim.x
   ln -sf $src/${exefile} .
fi

namelist=namelist
setup=setup.ksh
csvfile=radar_sim_dev.csv
ln -sf ../fix/${csvfile} .
ln -sf ../fix/l2rwbufr.table.csv .
ln -sf ../fix/l2rwbufr.table .

# Create the jobs for simulating the observations now.
# EXAMPLE:
# FHMIN_GFS=6, FHMAX_GFS=18, FHOUT_GFS=1, FHOUT_GDAS=6
# This example computes 3 simulated drw bufr files: one for 06z, 12z, and 18z.

#23,24,25,26, 27, 28, 29, 30 (September 23-30)
#000,024,048,072,096,120,144,168
#023,047,071,095,119,143,167,171
FHMIN_GFS=120 #start hour (6, 12, 18, 24, 30, 36, etc). 
FHMAX_GFS=143 #end  hour (6, 12, 18, 24, 30, etc). FHMIN = FHMAX if simulating obs for single cyc (1 bufr file).
FHOUT_GFS=1   #increment hour by (1; hourly 4DEnVar). Probably don't want to change this.
FHOUT_GDAS=1  #increment for gdas analysis (6; analysis every 6 hours @ 00, 06, etc.). Don't change.
fcsthr=$FHMIN_GFS
typeset -Z3 fcsthr
typeset -Z3 fcsthrm0

#The next 4 lines are for the NATURE RUN. There is only one NATURE, the config used
#forecast only (gfs), started at 00z on 20180911, and is named NATURE-SDATE-EDATE.
CDUMP="gfs"
valtime=2018092300 #this is where the experiment starts.
valtime=`${ndate} +$FHMIN_GFS $valtime` #add the first x-hours given by FHMIN_GFS to original start time.
PDY=`echo $valtime | cut -c 1-8`
CYC=`echo $valtime | cut -c 9-10`
EXP="NATURE-2018092300-2018100700" #This should be the NATURE RUN

while [[ $fcsthr -le $FHMAX_GFS ]]; do
   valcyc=`echo $valtime | cut -c 9-10`
   valpdy=`echo $valtime | cut -c 1-8`
   fcsthrm0=$fcsthr
   filename_tm00="${CDUMP}.t${CYC}z.atmf${fcsthrm0}.nemsio" #Get filename change
   
   cp $namelist ./simnml/namelist.atmf${fcsthr}
   sed -i "s/@filename_tm00@/${filename_tm00}/g" ./simnml/namelist.atmf${fcsthr}
   sed -i "s/@CDUMP@/${CDUMP}/g"       ./simnml/namelist.atmf${fcsthr}
   sed -i "s/@PDY@/${valpdy}/g"        ./simnml/namelist.atmf${fcsthr}
   sed -i "s/@CYC@/${CYC}/g"           ./simnml/namelist.atmf${fcsthr}
   sed -i "s/@EXP@/${EXP}/g"           ./simnml/namelist.atmf${fcsthr}

   cp $setup ./simsetup/setup_atmf${fcsthr}.ksh 
   sed -i "s/@namelist.atmfxxx@/namelist.atmf${fcsthr}/g" ./simsetup/setup_atmf${fcsthr}.ksh
   sed -i "s/@atmfxxx@/atmf${fcsthr}/g"                   ./simsetup/setup_atmf${fcsthr}.ksh
   sed -i "s/@fcsthr@/${fcsthr}/g"                        ./simsetup/setup_atmf${fcsthr}.ksh
   work=`pwd`
   sed -i "s#@work@#${work}#g"                            ./simsetup/setup_atmf${fcsthr}.ksh
   sbatch ./simsetup/setup_atmf${fcsthr}.ksh
   #ksh ./simsetup/setup_atmf${fcsthr}.ksh
   (( fcsthr=fcsthr+${FHOUT_GDAS} ))
   valtime=`${ndate} +$FHOUT_GFS $valtime` #now update by FHOUT_GFS.
done
