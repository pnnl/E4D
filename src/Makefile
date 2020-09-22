PETSC_DIR  ?= ~/codes/petsc
PETSC_ARCH ?= arch-linux2-c-opt
INSTALL_PREFIX=/usr

f95 = ${PETSC_DIR}/${PETSC_ARCH}/bin/mpif90
libs = ${PETSC_LIB}


CFLAGS= ${PETSC_CC_INCLUDES}

# default configuration options
ENABLE_OPT=0
ENABLE_DEBUG=0
OPT_FLAGS=-O2
MODE=res

# base compiler flags 

FFLAGS= $(PETSC_FC_INCLUDES)

# evaluate configuration options
ifeq ($(MODE), res)
	FFLAGS += -D resmode 
else ifeq ($(MODE), ini)
	FFLAGS += -D inimode
endif

ifeq ($(ENABLE_PETSC7), 1)
	FFLAGS += -D petsc_7 
endif 

ifeq ($(ENABLE_OPT), 1)
	FFLAGS += $(OPT_FLAGS)
endif
ifeq ($(ENABLE_DEBUG), 1)
	FFLAGS += -g
endif


include ${PETSC_DIR}/lib/petsc/conf/variables

sources := $(wildcard *.F90) $(wildcard *.F) $(wildcard *.f90)
objects := $(sources:%.F90=%.o)
objects := $(objects:%.f90=%.o)
objects := $(objects:%.F=%.o)


#Evaluate Configuration Options


# Release Target (Default)
release: FFLAGS += $(OPT_FLAGS)
release: e4d

# Debug Target
debug: FFLAGS += -g
debug: e4d

# Executable

e4d: $(objects)
	$(f95) $(objects) -o $@ $(libs)

READ_INP.o: private FFLAGS += -D sip

# Generate Dependencies

dependencies.d: makedepend.py $(sources)
	python3 ./makedepend.py $(sources)

# Include generated dependencies
include dependencies.d



# Generic Compilation Rules

%.o: %.f90
	$(f95) $(FFLAGS) -c $< -o $@ 

%.o: %.F
	$(f95) $(FFLAGS) -c $< -o $@ 

%.o: %.F90
	$(f95) $(FFLAGS) -c $< -o $@ 

# Install target
install: e4d
	install e4d $(INSTALL_PREFIX)/bin

# Clean target
clean: 
	-rm $(objects)
	-rm e4d
	-rm *.mod
	-rm dependencies.d

