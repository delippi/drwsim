#!/bin/ksh
#PBS -N drwsim
#PBS -l walltime=0:30:00
#PBS -l nodes=1:ppn=20
#PBS -q debug
#PBS -A blender
#PBS -o log_drwsim.out
#PBS -j oe

#cd /scratch4/NCEPDEV/meso/save/Donald.E.Lippi/PhD-globalOSSE/obssim/fortran/run
#src=/scratch4/NCEPDEV/meso/save/Donald.E.Lippi/PhD-globalOSSE/obssim/practice
src=../src

cd $src
make clean
make
cd -

exefile=drwsim.x
namelist=namelist
csvfile=radar_sim_dev.csv

ln -sf $src/${exefile} .
ln -sf ../fix/${csvfile} .
ln -sf ../fix/l2rwbufr.table.csv .
ln -sf ../fix/l2rwbufr.table .

$exefile < $namelist

