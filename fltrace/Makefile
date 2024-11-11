# Set compiler, flags and netcdf library according to machine we are using:

$(info Machine name: $(shell hostname))

ifeq ($(shell hostname),brillouin.dur.ac.uk)
archie_flag = false
MPIF90 ?= /usr/lib64/openmpi/bin/mpif90
FFLAGS = -O3 -Wuninitialized -march=native -fimplicit-none -Wall -Wextra -ffast-math -funroll-loops --param max-unroll-times=5
NETCDF = -I /usr/lib64/gfortran/modules
NETCDFLIB = -L/usr/lib64/libnetcdff.so.7 -lnetcdff
MODULEFLAG = -module $(OBJDIR)
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

TARGET = fltrace


# --------------------------------------------------
# Shouldn't need to touch below here
# --------------------------------------------------

SRCFILES = shared_data.f90 grid.f90 fltrace.f90

ifeq ($(archie_flag),true)
$(info Compiling on Archie-West)

FFLAGS += $(MODULEFLAG)
LDFLAGS = $(FFLAGS)


# Set pre-processor defines
DEFINES := $(DEFINE)

# --------------------------------------------------
# Shouldn't need to touch below here
# --------------------------------------------------

all: main

SRCDIR = src
OBJDIR = obj
BINDIR = bin
FC = $(MPIF90)
DATE := $(shell date +%s)
MACHINE := $(shell uname -n)
PREPROFLAGS = $(DEFINES) $(D)_COMMIT='"$(COMMIT)"' $(D)_DATE=$(DATE) \
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
	$(FC) -c $(FFLAGS) $(NETCDF) -module $(OBJDIR) -o $(OBJDIR)/$@ $<

%.o: %.F90
	$(FC) -c $(FFLAGS) $(NETCDF) -module $(OBJDIR) -o $(OBJDIR)/$@ $(PREPROFLAGS) $<

$(FULLTARGET): $(OBJFILES)
	@mkdir -p $(BINDIR)
	$(FC) -o $@ $(addprefix $(OBJDIR)/,$(OBJFILES)) $(LDFLAGS) $(NETCDFLIB)

clean:
	@rm -rf $(BINDIR) $(OBJDIR)

cleanall: tidy

tidy:
	@rm -rf $(OBJDIR) *~ *.pbs.* *.sh.* $(SRCDIR)/*~ *.log
	$(MAKE) -C $(SDF) cleanall

datatidy:
	@rm -rf Data/*

tarball:
	@sh $(SRCDIR)/make_tarball.sh


$(OBJFILES): | $(OBJDIR)

$(OBJDIR):
	@mkdir -p $(OBJDIR)

commit: FORCE
	@sh $(SRCDIR)/gen_commit_string.sh && $(MAKE) $(MAKECMDGOALS) DONE_COMMIT=1

FORCE:

.PHONY: commit clean cleanall tidy datatidy visit visitclean main FORCE

else

$(info Compiling locally or on Hamilton)


all: main

SDF := SDF/FORTRAN
SDFMOD = $(SDF)/include/sdf.mod
SRCDIR = src
OBJDIR = obj
BINDIR = bin
FC = $(MPIF90)
DATE := $(shell date +%s)
MACHINE := $(shell uname -n)
PREPROFLAGS = $(DEFINES) $(D)_COMMIT='"$(COMMIT)"' $(D)_DATE=$(DATE) \
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
	$(FC) -c $(FFLAGS) $(NETCDF) -J $(OBJDIR) -o $(OBJDIR)/$@ $(PREPROFLAGS) $<

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

shared_data.o: shared_data.f90
grid.o: grid.f90 shared_data.o
fltrace.o: fltrace.f90 shared_data.o grid.o


