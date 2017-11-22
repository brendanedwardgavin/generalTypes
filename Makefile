#your shared BLAS library goes here:
atlaslibs=$(HOME)/ATLAS/lib/libsatlas.so
all:
	gfortran matmul.f90 $(atlaslibs) -fPIC -shared -o mymatmul.so
