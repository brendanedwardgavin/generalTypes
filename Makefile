#your shared BLAS library goes here:
blasSharedLib=/usr/lib/libblas.so #$(HOME)/ATLAS/lib/libsatlas.so
mpiSharedLib=/home/brendan/software/openmpi/lib/libmpi_mpifh.so
mpifort=/home/brendan/software/openmpi/bin/mpif90
mklloc=/opt/intel/compilers_and_libraries_2017.1.132/linux/mkl/lib/intel64_lin/libmkl_def.so
all:
	ifort $(mklloc) matmul.f90 -mkl -fPIC -g -O2 -shared -o mymatmul.so
	$(mpifort) mpiwrap.f90 $(mpiSharedLib) -fPIC -shared -o mpiwrap.so
	#-fPIC -g -fbacktrace -O2 -shared -L/local/binaries/gfortran-bindings/install/lib -lfutile-1 -lmpi_f90 -lmpi_f77 -lmpi -o mymatmul_futile.so
	$(mpifort) $(mklloc) matmul_futile.f90 $(mpiSharedLib) -mkl -fPIC -g -O2 -shared -I/home/brendan/software/bigdft/build/install/include -L/home/brendan/software/bigdft/build/install/lib -lfutile-1 -lyaml -ldl -lrt -o mymatmul_futile.so
