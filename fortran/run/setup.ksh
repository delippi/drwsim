#!/bin/ksh
#PBS -N drwsim.@atmfxxx@
#PBS -l walltime=1:30:00
#PBS -l nodes=1:ppn=20
#PBS -q batch
#PBS -A fv3-cpu
#PBS -o ./simlog/log_drwsim.@atmfxxx@.out
#PBS -j oe

work="@work@"
cd $work
 
#exefile=drwsim.x
#namelist=@namelist.atmfxxx@
#$exefile < $namelist

drwsim.x @fcsthr@ < ./simnml/@namelist.atmfxxx@
