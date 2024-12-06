# Set compiler, flags and netcdf library according to machine we are using:



TARGET = lare3d
$(info Machine name: $(shell hostname))

ifeq ($(shell hostname),brillouin.dur.ac.uk)
archie_flag = false
MPIF90 ?= /usr/lib64/openmpi/bin/mpifort
NETCDF = -I /usr/lib64/gfortran/modules
NETCDFLIB = -L/usr/lib64/libnetcdff.so.7 -lnetcdff
FFLAGS = -O3 #-fcheck=all -Wuninitialized -march=native -fimplicit-none -Wall -Wextra -ffast-math -funroll-loops --param max-unroll-times=5

else ifeq ($(shell hostname),login1.ham8.dur.ac.uk)
archie_flag = false
MPIF90 ?= mpif90
FFLAGS = -O3 -Wuninitialized -march=native -fimplicit-none -Wall -Wextra -ffast-math -funroll-loops --param max-unroll-times=5
NETCDF = -I /apps/developers/libraries/netcdf/4.8.1/1/gcc-11.2-openmpi-4.1.1/include
NETCDFLIB = -L/apps/developers/libraries/netcdf/4.8.1/1/gcc-11.2-openmpi-4.1.1/lib  -lnetcdff
MODULEFLAG = -module $(OBJDIR)
else ifeq ($(shell hostname),login2.ham8.dur.ac.uk)
archie_flag = false
MPIF90 ?= mpif90
FFLAGS = -O3 -Wuninitialized -march=native -fimplicit-none -Wall -Wextra -ffast-math -funroll-loops --param max-unroll-times=5
NETCDF = -I /apps/developers/libraries/netcdf/4.8.1/1/gcc-11.2-openmpi-4.1.1/include
NETCDFLIB = -L/apps/developers/libraries/netcdf/4.8.1/1/gcc-11.2-openmpi-4.1.1/lib  -lnetcdff
MODULEFLAG = -module $(OBJDIR)
else
archie_flag = true
MPIF90 ?= mpiifort
FFLAGS = -O2
NETCDF = -I /opt/software/netcdf/intel-2020.4/fortran-4.5.4/include
NETCDFLIB = -L/opt/software/netcdf/intel-2020.4/fortran-4.5.4/lib  -lnetcdff
MODULEFLAG = -I/usr/include -I $(OBJDIR)
endif
# --------------------------------------------------
# Shouldn't need to touch below here
# --------------------------------------------------
SRCFILES = shared_data.o mpi_routines.o openboundary.o mpiboundary.o boundary.o oliver_output.o diagnostics.o  setup.o lagran.o  \
 remap.o xremap.o yremap.o zremap.o observables.o emergence.o strings.o initial_conditions.o welcome.o read_data.o lare3d.o

ifeq ($(archie_flag),true)
$(info Compiling on Archie-West)

SRCDIR = src
OBJDIR = obj
BINDIR = bin
MODULEFLAG = -module
MACHINEFLAGS = $(COPSON) $(MHDCLUSTER)
OPFLAGS = $(QMONO)
FC= mpiifort $(OPFLAGS)
PREPROFLAGS = $(NONMPIIO)

OBJFILES := $(SRCFILES:.f90=.o)
OBJFILES := $(OBJFILES:.F90=.o)

FULLTARGET = $(BINDIR)/$(TARGET)

#vpath %.f90 $(SRCDIR)
#vpath %.o $(OBJDIR)
VPATH = $(SRCDIR):$(OBJDIR)

# Rule to build the fortran files

%.o: %.f90
	@mkdir -p $(BINDIR) $(OBJDIR)
	$(FC) -c $(FFLAGS)  $(NETCDF) $(MODULEFLAG) $(OBJDIR) -o $(OBJDIR)/$@ $<

%.o: %.F90
	@mkdir -p $(BINDIR) $(OBJDIR)
	$(FC) -c $(FFLAGS)  $(NETCDF) $(MODULEFLAG) $(OBJDIR) -o $(OBJDIR)/$@ $(PREPROFLAGS) $<

$(FULLTARGET): $(OBJFILES)
	$(FC) $(FFLAGS) $(MODULEFLAG) $(OBJDIR) -o $@ $(addprefix $(OBJDIR)/,$(OBJFILES)) $(NETCDFLIB)

.PHONEY: clean
clean:
	@rm -rf *~ $(BINDIR) $(OBJDIR) *.pbs.* *.sh.* $(SRCDIR)/*~ *.log

.PHONEY: tidy
tidy:
	@rm -rf $(OBJDIR) *.pbs.* *.sh.* $(SRCDIR)/*~ *.log

.PHONEY: datatidy
datatidy:
	@rm -rf Data/*

.PHONEY: visit
visit:
	@cd VisIT;xml2makefile -clobber l3dv2.xml;make

.PHONEY: visitclean
visitclean:
	@cd VisIT;make clean;rm -f .depend


else

$(info Compiling Locally)

all: main

SDF := SDF/FORTRAN
SDFMOD = $(SDF)/include/sdf.mod
SRCDIR = src
OBJDIR = obj
BINDIR = bin
FC = $(MPIF90)
DATE := $(shell date +%s)
MACHINE := $(shell uname -n)
  $(D)_MACHINE='"$(MACHINE)"'


OBJFILES := $(SRCFILES:.f90=.o)
OBJFILES := $(OBJFILES:.F90=.o)

FULLTARGET = $(BINDIR)/$(TARGET)

VPATH = $(SRCDIR):$(SRCDIR)/core:$(SDF)/src:$(OBJDIR)

-include $(SRCDIR)/COMMIT

ifeq ($(DONE_COMMIT),)
main: commit
else
main: $(FULLTARGET)
endif

# Rule to build the fortran files

%.o: %.f90
	@mkdir -p $(BINDIR) $(OBJDIR)
	$(FC) -c $(FFLAGS) $(NETCDF) -J $(OBJDIR) -o $(OBJDIR)/$@ $<

%.o: %.F90
	@mkdir -p $(BINDIR) $(OBJDIR)
	$(FC) -c $(FFLAGS) $(NETCDF) -J $(OBJDIR) -o $(OBJDIR)/$@ $<

$(FULLTARGET): $(OBJFILES)
	$(FC) $(FFLAGS) -J $(OBJDIR) -o $@ $(addprefix $(OBJDIR)/,$(OBJFILES)) $(NETCDFLIB)


clean:
	@rm -rf $(BINDIR) $(OBJDIR)


#$(OBJFILES): | $(OBJDIR)

$(OBJDIR):
	@mkdir -p $(OBJDIR)

commit: FORCE
	@sh $(SRCDIR)/gen_commit_string.sh && $(MAKE) $(MAKECMDGOALS) DONE_COMMIT=1

FORCE:

.PHONY: commit clean cleanall tidy datatidy visit visitclean main FORCE

endif

shared_data.o:shared_data.F90
oliver_output.o:shared_data.o
mpi_routines.o:mpi_routines.f90 shared_data.o
setup.o:setup.F90 shared_data.o
openboundary.o:openboundary.f90 shared_data.o
mpiboundary.o:mpiboundary.f90 shared_data.o
boundary.o:boundary.f90 openboundary.o mpiboundary.o shared_data.o
xremap.o:xremap.f90 shared_data.o boundary.o
yremap.o:yremap.f90 shared_data.o boundary.o
zremap.o:zremap.f90 shared_data.o boundary.o
diagnostics.o:diagnostics.F90 shared_data.o boundary.o
lagran.o:lagran.F90 shared_data.o boundary.o emergence.o diagnostics.o
remap.o:remap.f90 shared_data.o xremap.o yremap.o zremap.o
observables.o:observables.F90 shared_data.o
emergence.o:emergence.f90 shared_data.o boundary.o
strings.o:strings.f90 shared_data.o
initial_conditions.o:initial_conditions.f90 strings.o shared_data.o
welcome.o:welcome.f90 shared_data.o
read_data.o: read_data.f90 shared_data.o
lare3d.o:lare3d.f90 shared_data.o setup.o boundary.o diagnostics.o lagran.o remap.o mpi_routines.o observables.o emergence.o initial_conditions.o welcome.o


