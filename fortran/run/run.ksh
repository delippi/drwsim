#!/bin/ksh

cd ../src
make clean
make
cd -

#cp ../src/drwsim.x .
#ln -s ../src/drwsim.x .

cp ../fix/small_list_for_radar_sim_dev.csv .
cp ../fix/l2rwbufr.table.csv .
cp ../fix/l2rwbufr.table .

EXECNAME=./drwsim.x
namelist=./namelist

$EXECNAME < $namelist
