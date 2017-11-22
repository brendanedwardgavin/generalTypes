#your shared BLAS library goes here:
atlaslibs=/usr/lib/libblas.so #$(HOME)/ATLAS/lib/libsatlas.so
all:
	gfortran matmul.f90 $(atlaslibs) -fPIC -g -fbacktrace -O2 -shared -o mymatmul.so
