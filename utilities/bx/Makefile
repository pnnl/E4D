f95 = gfortran

EXODUS = /home/john775/codes/e4d_dev/third_party/exodus-6.02/exodus
NETCDF_LIB = -L/usr/lib64 -lnetcdf 
includes = -I$(EXODUS)/forbind/include -I/usr/include/ -I$(EXODUS)/include
libs = -lm ${EXODUS}/libexoIIv2for.a ${EXODUS}/libexodus.a ${NETCDF_LIB}
objects = subroutines.o exo_subs.o  build_exo.o


bx : $(objects)
	$(f95) $(objects) $(includes) $(libs) -o bx  

build_exo.o : build_exo.f90
	$(f95) -c -static build_exo.f90 -o  build_exo.o

exo_subs.o : exo_subs.f90
	$(f95) $(includes) -c -static exo_subs.f90 -o  exo_subs.o

subroutines.o : subroutines.f90
	$(f95) -c -static subroutines.f90 -o  subroutines.o


