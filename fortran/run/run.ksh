
# Check for existence of sim directories which contain bufr, log, nml, and setup files.
# If they don't already exist, create them. Otherwise continues.
dirs="bufr log nml setup"
for dir in $dirs; do
   if [[ ! -e ./sim${dir} ]];then
      mkdir ./sim${dir}
   fi
done

# Compile the simulation code.
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


# Create the jobs for simulating the observations now.
# EXAMPLE:
# FHMIN_GFS=6, FHMAX_GFS=18, FHOUT_GFS=1, FHOUT_GDAS=6
# This example computes 3 simulated drw bufr files: one for 06z, 12z, and 18z.
FHMIN_GFS=6   #start hour (6, 12, 18, 00). 
FHMAX_GFS=18  #end   hour (6, 12, 18, 00). FHMIN = FHMAX if simulating obs for single cycle (1 bufr file).
FHOUT_GFS=1   #increment hour by (1; hourly 4DEnVar). Probably don't want to change this.
FHOUT_GDAS=6  #increment for gdas analysis (6; analysis every 6 hours @ 00, 06, etc.). Don't change.
fcsthr=$FHMIN_GFS
typeset -Z3 fcsthr
typeset -Z3 fcsthrm3
typeset -Z3 fcsthrm2
typeset -Z3 fcsthrm1
typeset -Z3 fcsthrm0
typeset -Z3 fcsthrp1
typeset -Z3 fcsthrp2
typeset -Z3 fcsthrp3

CDUMP="gfs"
CYC="00"
PDY="20180911"
EXP="NATURE-2018091100-2018091800" #This should be the NATURE RUN

while [[ $fcsthr -le $FHMAX_GFS ]]; do
   (( fcsthrm3=fcsthr - 3 ))
   (( fcsthrm2=fcsthr - 2 ))
   (( fcsthrm1=fcsthr - 1 ))
   fcsthrm0=$fcsthr
   (( fcsthrp1=fcsthr + 1 ))
   (( fcsthrp2=fcsthr + 2 ))
   (( fcsthrp3=fcsthr + 3 ))
   #filename="gfs.t00z.atmf${fcsthr}.nemsio" #Get filename change
   filename3="${CDUMP}.t${CYC}z.atmf${fcsthrm3}.nemsio" #Get filename change
   filename4="${CDUMP}.t${CYC}z.atmf${fcsthrm2}.nemsio" #Get filename change
   filename5="${CDUMP}.t${CYC}z.atmf${fcsthrm1}.nemsio" #Get filename change
   filename6="${CDUMP}.t${CYC}z.atmf${fcsthrm0}.nemsio" #Get filename change
   filename7="${CDUMP}.t${CYC}z.atmf${fcsthrp1}.nemsio" #Get filename change
   filename8="${CDUMP}.t${CYC}z.atmf${fcsthrp2}.nemsio" #Get filename change
   filename9="${CDUMP}.t${CYC}z.atmf${fcsthrp3}.nemsio" #Get filename change
   cp $namelist ./simnml/namelist.atmf${fcsthr}
   sed -i "s/@filename3@/${filename3}/g" ./simnml/namelist.atmf${fcsthr}
   sed -i "s/@filename4@/${filename4}/g" ./simnml/namelist.atmf${fcsthr}
   sed -i "s/@filename5@/${filename5}/g" ./simnml/namelist.atmf${fcsthr}
   sed -i "s/@filename6@/${filename6}/g" ./simnml/namelist.atmf${fcsthr}
   sed -i "s/@filename7@/${filename7}/g" ./simnml/namelist.atmf${fcsthr}
   sed -i "s/@filename8@/${filename8}/g" ./simnml/namelist.atmf${fcsthr}
   sed -i "s/@filename9@/${filename9}/g" ./simnml/namelist.atmf${fcsthr}
   sed -i "s/@CDUMP@/${CDUMP}/g"       ./simnml/namelist.atmf${fcsthr}
   sed -i "s/@PDY@/${PDY}/g"           ./simnml/namelist.atmf${fcsthr}
   sed -i "s/@CYC@/${CYC}/g"           ./simnml/namelist.atmf${fcsthr}
   sed -i "s/@EXP@/${EXP}/g"           ./simnml/namelist.atmf${fcsthr}

   cp $setup ./simsetup/setup_atmf${fcsthr}.ksh 
   sed -i "s/@namelist.atmfxxx@/namelist.atmf${fcsthr}/g" ./simsetup/setup_atmf${fcsthr}.ksh
   sed -i "s/@atmfxxx@/atmf${fcsthr}/g"                   ./simsetup/setup_atmf${fcsthr}.ksh
   sed -i "s/@fcsthr@/${fcsthr}/g"                        ./simsetup/setup_atmf${fcsthr}.ksh
   work=`pwd`
   sed -i "s#@work@#${work}#g"                            ./simsetup/setup_atmf${fcsthr}.ksh
   qsub ./simsetup/setup_atmf${fcsthr}.ksh
   #ksh ./simsetup/setup_atmf${fcsthr}.ksh
   (( fcsthr=fcsthr+${FHOUT_GDAS} ))
done
