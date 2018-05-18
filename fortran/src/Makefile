BACIO_VER = 2.0.1
BUFR_VER = 10.2.5
NEMSIO_VER = 2.2.1
SFCIO_VER = 1.0.0
SIGIO_VER =2.0.1
SP_VER = 2.0.2
W3EMC_VER = 2.0.5
W3NCO_VER = 2.0.6

CORELIB=/scratch3/NCEPDEV/nwprod/lib

INCsfcio = $(CORELIB)/sfcio/v$(SFCIO_VER)/incmod/sfcio_v$(SFCIO_VER)_4
INCsigio = $(CORELIB)/sigio/v$(SIGIO_VER)/incmod/sigio_v$(SIGIO_VER)_4
INCnemsio= $(CORELIB)/nemsio/v$(NEMSIO_VER)/incmod/nemsio_v$(NEMSIO_VER)
INCw3    = $(CORELIB)/w3emc/v$(W3EMC_VER)/incmod/w3emc_v$(W3EMC_VER)_d

BACIO_LIB4=$(CORELIB)/bacio/v$(BACIO_VER)/libbacio_v$(BACIO_VER)_4.a
BUFR_LIBd=$(CORELIB)/bufr/v$(BUFR_VER)/libbufr_v$(BUFR_VER)_d_64.a
NEMSIO_LIB=$(CORELIB)/nemsio/v$(NEMSIO_VER)/libnemsio_v$(NEMSIO_VER).a
SFCIO_LIB=$(CORELIB)/sfcio/v$(SFCIO_VER)/libsfcio_v$(SFCIO_VER)_4.a
SIGIO_LIB=$(CORELIB)/sigio/v$(SIGIO_VER)/libsigio_v$(SIGIO_VER)_4.a
SP_LIBd=$(CORELIB)/sp/v$(SP_VER)/libsp_v$(SP_VER)_d.a
W3EMC_LIBd=$(CORELIB)/w3emc/v$(W3EMC_VER)/libw3emc_v$(W3EMC_VER)_d.a
W3NCO_LIBd=$(CORELIB)/w3nco/v$(W3NCO_VER)/libw3nco_v$(W3NCO_VER)_d.a

INCnetcdf  = ${NETCDF}/include
LIBnetcdf  = -L${NETCDF}/lib -lnetcdf


CF	= mpif90 -f90=ifort -warn unused -I $(INCsfcio) -I $(INCsigio) -I $(INCnemsio) \
	-I $(INCnetcdf) -I $(INCw3) -O3 -traceback -convert big_endian -assume byterecl \
	-g -fp-model source 

CFFLAG	= $(CF) 


LIBS	= $(SP_LIBd) $(NEMSIO_LIB) $(BACIO_LIB4) $(SIGIO_LIB) \
          $(SFCIO_LIB) $(BUFR_LIBd) $(W3NCO_LIBd) $(W3EMC_LIBd) $(LIBnetcdf)

execname= drwsim.x

drwsim.x: kinds.o constants.o grdcrd.o interp_grid_to_ob_util_nio.o drwsim_nio.f90
	$(CFFLAG) kinds.o constants.o grdcrd.o interp_grid_to_ob_util_nio.o -o $(execname) drwsim_nio.f90 $(LIBS)

interp_grid_to_ob_util_nio.o: kinds.o constants.o grdcrd.o interp_grid_to_ob_util_nio.f90
	$(CF) -c interp_grid_to_ob_util_nio.f90 #$(LIBS)

grdcrd.o: kinds.o constants.o grdcrd.f90
	$(CF) -c grdcrd.f90

constants.o: kinds.o constants.f90
	$(CF) -c constants.f90

kinds.o: kinds.f90
	$(CF) -c kinds.f90

clean:
	rm -f *.o *.x *.mod

