FF = ifort
XFLAGS = -O -I $(NETCDF_ROOT)/include
LIBS = -L $(NETCDF_ROOT)/lib -L $(HDF5_HOME)/lib -lnetcdf -lnetcdff -lhdf5 -lhdf5_hl
LDFLAGS = 

OBJT = smclim.o smread.o setxyz_m.o ccinterp.o readswitch.o jimcc_m.o \
       latltoij_m.o ncwrite.o ncread.o xyzinfo_m.o newmpar_m.o \
       indices_m.o parm_m.o precis_m.o ind_m.o jimco_m.o jim_utils.o \
       nfft_m.o misc.o

smclim :$(OBJT)
	$(FF) $(XFLAGS) $(OBJT) $(LDFLAGS) $(LIBS) -o smclim

clean:
	rm *.o core *.mod 
# This section gives the rules for building object modules.

.SUFFIXES:.f90
.f90.o:
	$(FF) -c $(XFLAGS) $<
.f.o:
	$(FF) -c $(XFLAGS) $<

smclim.o : ccinterp.o
smread.o : ccinterp.o
ccinterp.o : ccinterp.f90 setxyz_m.o xyzinfo_m.o latltoij_m.o newmpar_m.o
latltoij_m.o : latltoij_m.f90 xyzinfo_m.o newmpar_m.o
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
