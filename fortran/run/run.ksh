
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
FHMAX_GFS=24  #end hour
FHMIN_GFS=0   #start hour
FHOUT_GFS=1   #increment hour by
fcsthr=$FHMIN_GFS
typeset -Z3 fcsthr

CDUMP="gfs"
CYC="00"
PDY="20180911"
EXP="NATURE-2018091100-2018091800" #This should be the NATURE RUN

while [[ $fcsthr -le $FHMAX_GFS ]]; do
   #filename="gfs.t00z.atmf${fcsthr}.nemsio" #Get filename change
   filename="${CDUMP}.t${CYC}z.atmf${fcsthr}.nemsio" #Get filename change
   cp $namelist ./simnml/namelist.atmf${fcsthr}
   sed -i "s/@filename@/${filename}/g" ./simnml/namelist.atmf${fcsthr}
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
   (( fcsthr=fcsthr+${FHOUT_GFS} ))
done
