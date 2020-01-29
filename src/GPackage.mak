# This is an include file for using HydroGrid with LBNL's BoxLib/amrex codes
# It assumes HYDROLIB_HOME is set to the repo stochasticHydroTools/HydroGrid

vpath %.f90 $(HYDROLIB_HOME)/src
vpath %.c $(HYDROLIB_HOME)/src
vpath %.h $(HYDROLIB_HOME)/src
INCLUDE_LOCATIONS += $(HYDROLIB_HOME)/src

myf90sources += Precision.f90
myf90sources += FFTW.f90
myf90sources += VisitWriter.f90
myf90sources += HydroGridModule.f90
myf90sources += HydroGridCInterface.f90

# Note that the RNGs are already included in the Parallel version of BoxLib, but not fParallel
# They are extracted here in C and made separate so there is no inconsistencies
ifdef USE_HG_RNGs
   myf90sources += Random.f90
   myf90sources += NURNGs.f90
   myf90sources += RNG.f90
   myf90sources += RNGEngine.f90
   mycsources += RNGs.c visit_writer.c
   mycheaders += RNGs.h visit_writer.h
endif

mycsources += visit_writer.c
mycheaders += visit_writer.h

# For fParallel-based codes:
f90sources += $(myf90sources) 
csources += $(mycsources)
cheaders += $(mycheaders) # This does not actually seem to exist in fParallel?!?

# For Parallel-based codes:
f90EXE_sources += $(myf90sources)
cEXE_sources   += $(mycsources)
cEXE_headers   += $(mycheaders)
