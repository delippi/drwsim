#!/bin/ksh

rm -f *.pyf *.so *.mod *.o

#name='kinds'
#f2py -m ${name} -h ${name}.pyf ${name}.f90 constants.f90 interp.f90 grdcrd.f90 nc2rw.f90
#gfortran -shared -fPIC -o kinds.so kinds.f90
#gfortran -shared -fPIC -o constants.so constants.f90 
#gfortran -c constants.f90 kinds.o
#grdcrd
name='grdcrd'
f2py -m ${name} -h ${name}.pyf ${name}.f90 constants.f90
f2py -c ${name}.pyf            ${name}.f90 constants.f90
#interp
#name='interp'
#f2py -m ${name} -h ${name}.pyf ${name}.f90 constants.f90
#f2py -c ${name}.pyf            ${name}.f90 constants.f90
#nc2rw
#name='nc2rw'
#f2py -m ${name} -h ${name}.pyf ${name}.f90 constants.f90 interp.f90 
#f2py -c ${name}.pyf            ${name}.f90 constants.f90 interp.f90 








#f2py -m ${name} -h ${name}.pyf kinds.f90 ${name}.f90 
#f2py -c ${name}.pyf kinds.f90 ${name}.f90 

#f2py -m ${name} -h ${name}.pyf ${name}.f90 
#f2py -c ${name}.pyf ${name}.f90 


#f2py -m fmods -h fmods.pyf fmods.f90 kinds.f90 
#f2py -c fmods.pyf fmods.f90 kinds.f90
#f2py -m interp_util -h interp_util.pyf kinds.f90 constants.f90 interp_grid_to_ob_util.f90
#f2py -c interp_util.pyf kinds.f90 constants.f90 interp_grid_to_ob_util.f90



