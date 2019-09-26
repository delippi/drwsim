#!/bin/ksh
#SBATCH -J drwsim.@atmfxxx@
#SBATCH -t 0:20:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -q batch 
#SBATCH -A fv3-cpu
#SBATCH -o ./simlog/log_drwsim.@atmfxxx@.out

work="@work@"
cd $work
 
#exefile=drwsim.x
#namelist=@namelist.atmfxxx@
#$exefile < $namelist

drwsim.x @fcsthr@ < ./simnml/@namelist.atmfxxx@
