
ifneq ($(CUSTOM),yes)
FC = ifort
XFLAGS = -O 
LIBS = -L $(NETCDF_ROOT)/lib -lnetcdf -lnetcdff
INC = -I $(NETCDF_ROOT)/include
PPFLAG90 = -fpp
PPFLAG77 = -fpp
endif

ifeq ($(GFORTRAN),yes)
FC = gfortran
XFLAGS = -O2 -mtune=native -march=native -I $(NETCDF_ROOT)/include
PPFLAG90 = -x f95-cpp-input
PPFLAG77 = -x f77-cpp-input
endif

OBJT = smclim.o smread.o setxyz_m.o ccinterp.o readswitch.o jimcc_m.o \
       latltoij_m.o ncwrite.o ncread.o xyzinfo_m.o newmpar_m.o \
       indices_m.o parm_m.o precis_m.o ind_m.o jimco_m.o jim_utils.o \
       nfft_m.o misc.o netcdf_m.o

smclim :$(OBJT)
	$(FC) $(XFLAGS) $(OBJT) $(LIBS) -o smclim

clean:
	rm *.o core *.mod smclim
# This section gives the rules for building object modules.

.SUFFIXES:.f90
.f90.o:
	$(FC) -c $(XFLAGS) $(INC) $(PPFLAG90) $<
.f.o:
	$(FC) -c $(XFLAGS) $(INC) $(PPFLAG77) $<

# Remove mod rule from Modula 2 so GNU make doesn't get confused
%.o : %.mod

smclim.o : ccinterp.o netcdf_m.o
smread.o : ccinterp.o netcdf_m.o
ccinterp.o : ccinterp.f90 setxyz_m.o xyzinfo_m.o latltoij_m.o newmpar_m.o precis_m.o
latltoij_m.o : latltoij_m.f90 xyzinfo_m.o newmpar_m.o precis_m.o
setxyz_m.o : setxyz_m.f90 newmpar_m.o indices_m.o parm_m.o precis_m.o ind_m.o xyzinfo_m.o jimco_m.o jimcc_m.o 
xyzinfo_m.o : xyzinfo_m.f90 precis_m.o
newmpar_m.o : newmpar_m.f90 
precis_m.o : precis_m.f90
indices_m.o : indices_m.f90
parm_m.o : parm_m.f90 precis_m.o 
ind_m.o : ind_m.f90 newmpar_m.o 
jimcc_m.o : jimcc_m.f90 parm_m.o precis_m.o 
jimco_m.o : jimco_m.f90 precis_m.o jim_utils.o nfft_m.o 
jim_utils.o : jim_utils.f90 precis_m.o 
nfft_m.o : nfft_m.f90 precis_m.o 
ncread.o : netcdf_m.o
ncwrite.o : netcdf_m.o
