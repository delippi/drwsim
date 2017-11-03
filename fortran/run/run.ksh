#!/bin/ksh
#PBS -N drwsim
#PBS -l walltime=0:30:00
#PBS -l nodes=1:ppn=20
#PBS -q debug
#PBS -A blender
#PBS -o log_drwsim.out
#PBS -j oe


cd /scratch4/NCEPDEV/meso/save/Donald.E.Lippi/PhD-globalOSSE/obssim/fortran/run
cd ../src
make clean
make
cd -

exefile=drwsim.x
namelist=./namelist
FIXFILES="small_list_for_radar_sim_dev.csv 
          l2rwbufr.table.csv 
          l2rwbufr.table"
if ! [[ -L $exefile ]]; then
   echo "Linking $exefile" 
   ln -s ../${exefile}
else
   echo "$exefile is already linked"
fi

for fixfile in $FIXFILES; do
   if ! [[ -L $fixfile ]]; then
      echo "Linking $fixfile"
      ln -s ../fix/${fixfile}
   else
      echo "$fixfile is already linked"
   fi
done

./$exefile < $namelist
