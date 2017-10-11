#!/bin/ksh

rm -f *.pyf *.so *.mod *.o

name='gsimods'
f2py -m ${name} -h ${name}.pyf ${name}.f90
f2py -c ${name}.pyf            ${name}.f90



name='interp'
f2py -m ${name} -h ${name}.pyf constants.f90 ${name}.f90 
f2py -c --fcompiler=intelem ${name}.pyf            constants.f90 ${name}.f90 
#f2py -c -m ${name} constants.f90 ${name}.f90


