# AUTHOR
# Adriaan Larmuseau - 1ste master computerwetenschappen
#
# NOTES
# Nieuw in deze Makefile:
# g95 wordt gebruikt

#### MODE ###########################################################

# Select the mode: DEBUG or OPTIMIZE
MODE=OPTIMIZE

#### COMPILER #######################################################

# G95 (GCC 4.1.2 (g95 0.93!) Jun 16 2010)
FC=g95

# Basisvlaggen
FFLAGS_gfortran = -pedantic -Wall -Wunderflow
FFLAGS_g95 = -Wall -Wextra -pedantic

# Debugvlaggen (tijdens testfase)
FFLAGS_DEBUG_gfortran = -g -fbounds-check
FFLAGS_DEBUG_g95 = -fbounds-check -g

# Optimalisatievlaggen (ZELF IN TE VULLEN)
FFLAGS_OPTIMIZE_gfortran = -O4
FFLAGS_OPTIMIZE_g95 = -O3

# Stel de selectie van de juiste vlaggen samen op basis van MODE en FC
FFLAGS = $(FFLAGS_$(MODE)_$(FC)) $(FFLAGS_$(FC))

#### BLAS ###########################################################
# blas -- version 1.2-7
BLAS_FLAGS= -lblas

#### AUTOMATISCHE COMPILATIE *.f95 -> *.o ###########################

.SUFFIXES:
.SUFFIXES: .f95 $(SUFFIXES)

# Lees als: "Een object-bestand kan je genereren uit een
#            Fortran-bestand d.m.v. het volgende commando."
.f95.o:
	$(FC) -c $(FFLAGS) $<

#### TARGETS ########################################################

matrix_product:  zmatrix.o multiplication.o matrix_product.o
	$(FC) -o $@ $^ $(BLAS_FLAGS)

shortest_distance: zmatrix.o shortest_distance.o
	$(FC) -o $@ $^ 

generate: zmatrix.o generate_matrices.o
	$(FC) -o $@ $^

generate_edges : zmatrix.o generate_edges.o
	$(FC) -o $@ $^ 

testconv : zmatrix.o test_conversion.o
	$(FC) -o $@ $^ 

testmult : zmatrix.o multiplication.o test_multiply.o
	$(FC) -o $@ $^ $(BLAS_FLAGS)

testmultedges : zmatrix.o multiplication.o test_multiply_edges.o
	$(FC) -o $@ $^ $(BLAS_FLAGS)

test_distance : zmatrix.o test_distance.o
	$(FC) -o $@ $^ 

test_steps : zmatrix.o test_steps.o
	$(FC) -o $@ $^ 

testread : zmatrix.o test_read.o
	$(FC) -o $@ $^ 

perf_edges : zmatrix.o multiplication.o testperformance_edges.o
	$(FC) -o $@ $^ $(BLAS_FLAGS)

perf : zmatrix.o multiplication.o testperformance.o
	$(FC) -o $@ $^ $(BLAS_FLAGS)

clean:
	rm -f *.o *.mod *.data *.txt

veryclean: clean
	rm -f matrix_product shortest_distance generate_edges generate 
