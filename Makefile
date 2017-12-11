#your shared BLAS library goes here:
blasSharedLib=/usr/lib/libblas.so #$(HOME)/ATLAS/lib/libsatlas.so
mpiSharedLib=/local/gavin/OpenMPI/lib/libmpi_mpifh.so
mpifort=/local/gavin/OpenMPI/bin/mpif90
all:
	gfortran matmul.f90 $(blasSharedLib) -fPIC -g -fbacktrace -O2 -shared -o mymatmul.so
	$(mpifort) mpiwrap.f90 $(mpiSharedLib) -fPIC -shared -o mpiwrap.so
	-fPIC -g -fbacktrace -O2 -shared -L/local/binaries/gfortran-bindings/install/lib -lfutile-1 -lmpi_f90 -lmpi_f77 -lmpi -o mymatmul_futile.so
	#-fPIC -g -fbacktrace -O2 -shared -L/local/gavin/software/bigdft/build_dir/install/lib -lfutile-1 -o mymatmul_futile.so
