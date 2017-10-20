
#INCnetcdf=/apps/netcdf/3.6.3-intel/include

W3EMC_VER = 2.0.5
W3NCO_VER = 2.0.6
BUFR_VER = 10.2.5

CORELIB=/scratch3/NCEPDEV/nwprod/lib

INCnetcdf=$(NETCDF)/include
#INCw3    = $(CORELIB)/w3emc/v$(W3EMC_VER)/incmod/w3emc_v$(W3EMC_VER)_d
INCw3    = /scratch3/NCEPDEV/nwprod/lib/w3nco/v2.0.6/incmod/w3nco_v2.0.6_d 
#W3EMC_LIBd=$(CORELIB)/w3emc/v$(W3EMC_VER)/libw3emc_v$(W3EMC_VER)_d.a
W3NCO_LIBd=/scratch3/NCEPDEV/nwprod/lib
BUFR_LIBd=/scratch3/NCEPDEV/nwprod/lib

LIBnetcdf=$(NETCDF)/lib -lnetcdf

LIBS=	-L$(NETCDF)/lib -lnetcdf \
	-L$(W3NCO_LIBd) -lw3nco_v2.0.6_d \
	-L$(BUFR_LIBd)  -lbufr_v10.2.5_d_64

CF	= mpif90 -f90=ifort -I $(INCnetcdf) -I $(INCw3)
CFFLAG	= mpif90 -f90=ifort -I $(INCnetcdf) -I $(INCw3)
execname= drwsim.x

drwsim.x: kinds.o constants.o grdcrd.o interp_grid_to_ob_util.o drwsim.f90
	$(CFFLAG) kinds.o constants.o grdcrd.o interp_grid_to_ob_util.o -o $(execname) drwsim.f90 $(LIBS)

interp_grid_to_ob_util.o: kinds.o constants.o grdcrd.o interp_grid_to_ob_util.f90
	$(CF) kinds.o constants.o grdcrd.o -c interp_grid_to_ob_util.f90 $(LIBS)

grdcrd.o: kinds.o constants.o grdcrd.f90
	$(CF) kinds.o constants.o -c grdcrd.f90

constants.o: kinds.o constants.f90
	$(CF) kinds.o -c constants.f90

kinds.o: kinds.f90
	$(CF) -c kinds.f90

clean:
	rm -f *.o *.x *.mod

