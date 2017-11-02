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

cp ../fix/small_list_for_radar_sim_dev.csv .
cp ../fix/l2rwbufr.table.csv .
cp ../fix/l2rwbufr.table .

EXECNAME=./drwsim.x
namelist=./namelist

$EXECNAME < $namelist
