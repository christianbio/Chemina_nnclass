CC=gfortran
CFLAGS = -O3 -ffree-form -ffast-math 
MFLAG = -mcmodel='medium'
all: 
	$(CC) $(CFLAGS) -c parse_text.f
	$(CC) $(CFLAGS) -c nr_common.f95
	$(CC) $(CFLAGS) -c nr_mod.f95
	$(CC) $(CFLAGS) -c eigen_fit_common.f90
	$(CC) $(CFLAGS) -c neural.f95
	$(CC) $(CFLAGS) -L/home/user2/apps/lapack-3.10.1 eigen_fit.f90 neural.o nr_common.o nr_mod.o eigen_fit_common.o parse_text.o -o eigen_fit -llapack -lrefblas
parallel:
	$(CC) $(CFLAGS) -c parse_text.f
	$(CC) $(CFLAGS) -c nr_common.f95
	$(CC) $(CFLAGS) -c nr_mod.f95
	$(CC) $(CFLAGS) -c eigen_fit_common.f90
	$(CC) $(CFLAGS) -c neural.f95
	$(CC) $(CFLAGS) -fopenmp -L/home/user2/apps/lapack-3.10.1 eigen_fit.f90 neural.o nr_common.o nr_mod.o eigen_fit_common.o parse_text.o -o eigen_fit_parallel -llapack -lrefblas
eigen_model:
	$(CC) $(CFLAGS) -c nr_common.f95
	$(CC) $(CFLAGS) -c nr_mod.f95
	$(CC) $(CFLAGS) -c eigen_model_common.f90
	$(CC) $(CFLAGS) -fopenmp eigen_model.f90 nr_common.o nr_mod.o eigen_model_common.o parse_text.o -o eigen_model
findbonds:
	$(CC) $(CFLAGS) -c findbonds_common.f
	$(CC) $(CFLAGS) findbonds.f90 parse_text.o findbonds_common.o -o findbonds
calc_distances:
	$(CC) $(CFLAGS) -c calc_distances_common.f
	$(CC) $(CFLAGS) calc_distances.f90 parse_text.o calc_distances_common.o -o calc_distances
make_spline:
	$(CC) $(CFLAGS) -c make_spline_common.f
	$(CC) $(CFLAGS) make_spline.f nr_common.o nr_mod.o make_spline_common.o -o make_spline
nnmodel:
	$(CC) $(CFLAGS) -c nnmodel_common.f90
	$(CC) $(CFLAGS) nnmodel.f90 nnmodel_common.f90 nr_common.o nr_mod.o parse_text.o permutations.o -o nnmodel
find_limits:
	$(CC) $(CFLAGS) -c find_limits.f90
	$(CC) $(CFLAGS) find_limits.f90 parse_text.o -o find_limits
