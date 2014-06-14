FC=gfortran
# Choose a set of flags for either debugging purposes or else for optimized code
#FFLAGS=-O3 -fdefault-real-8 -fdefault-double-8 -cpp -fbounds-check
#FFLAGS=-Wall -f
FFLAGS=-O3 -fdefault-real-8 -fdefault-double-8 -cpp #-fopenmp

FFTWINC=-I/usr/local/fftw-3.3.3/include/
FFTWLIB=-L/usr/local/fftw-3.3.3/lib/

#OBJS = fftw_mod.o field_model.o gl_integrator.o fundamental_matrix.o get_floquet.o
#OBJS2 = fftw_mod.o field_model.o gl_integrator.o fundamental_matrix.o check_efunction.o
#OBJS3 = fftw_mod.o field_model_sinegordon.o gl_integrator.o initial_field_sg.o

OBJS = fftw_mod.o field_model.o hamiltonian.o symplectic_integrate.o fundamental_matrix.o get_floquet.o
OBJS2 = fftw_mod.o field_model.o hamiltonian.o symplectic_integrate.o fundamental_matrix.o check_efunction.o
OBJS3 = fftw_mode.o field_model.o hamiltonian.o symplectic_integrate.o initial_field.o

floquet: floquet

check_func: check_efunc

bound_states: bound_states

clean:
	rm -f *.o
	rm -f *.mod

floquet: %: $(OBJS)
	$(FC) $(FFLAGS) $(FFTWINC) $^ -o floquet $(FFTWLIB) -llapack -lfftw3 -lm

check_efunc: %: $(OBJS2)
	$(FC) $(FFLAGS) $(FFTWINC) $^ -o check_efunc $(FFTWLIB) -llapack -lfftw3 -lm

bound_states: %: $(OBJS3)
	$(FC) $(FFLAGS) $(FFTWINC) $^ -o bound_states $(FFTWLIB) -llapack -lfftw3 -lm

fundamental_matrix.o: fundamental_matrix.f90
	$(FC) $(FFLAGS) $(FFTWINC) -c $< -o $@ $(FFTWLIB) -llapack -lfftw3 -lm

%.o: %.f90
	$(FC) $(FFLAGS) $(FFTWINC) -c $< -o $@ $(FFTWLIB) -llapack -lfftw3 -lm

fftw_mod.o: fftw_mod.f90
	$(FC) $(FFLAGS) $(FFTWINC) -c $< -o $@ $(FFTWLIB) -llapack -lfftw3 -lm