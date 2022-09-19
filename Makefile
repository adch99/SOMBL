CC=gcc
# BLASDIR=../openblas
# BLASDIR=../../Build/openblas
# IDIRS=-I$(BLASDIR)/include
# LBLAS=-L$(BLASDIR)/lib -Wl,-rpath,$(BLASDIR)/lib -lopenblas -lpthread
LBLAS=-lcblas -lblas
CFLAGS=-Wall -Wextra -g -fdiagnostics-color=always -ffast-math $(IDIRS)
# CFLAGS=-O3 -ffast-math -fopenmp $(IDIRS) 
LFLAGS=-llapacke -lm $(LBLAS)
ERRORLOG=logs/compiler_error.log

_DEPS = utils/utils.c ham_gen/ham_gen.c params/params.c io/io.c
DEPS = $(patsubst %,src/%,$(_DEPS))

OBJ = $(patsubst %.c,build/%.o,$(_DEPS))

TESTS = $(wildcard tests/*.c)

default: exact_diag_simulation

all: exact_diag_simulation calculate_dist_vs_gfuncsq calculate_imbalance tests

bins: exact_diag_simulation calculate_dist_vs_gfuncsq calculate_imbalance

build/%.o: src/%.c
	$(CC) -c -o $@ $^ $(CFLAGS)

exact_diag_simulation: $(OBJ)
	$(CC) -o build/$@ src/$@.c $^ $(CFLAGS) $(LFLAGS)

calculate_dist_vs_gfuncsq: $(OBJ)
	$(CC) -o build/$@ src/$@.c $^ $(CFLAGS) $(LFLAGS)

calculate_imbalance: $(OBJ)
	$(CC) -o build/$@ src/$@.c $^ $(CFLAGS) $(LFLAGS)

output_hamiltonian: $(OBJ)
	$(CC) -o build/$@ src/$@.c $^ $(CFLAGS) $(LFLAGS)

build/tests/%: tests/%.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LFLAGS)

tests: $(patsubst tests/%.c,build/tests/%,$(TESTS))

clean:
	rm -rf build/*/*.o build/exact_diag_simulation
	find data -size 0 -print -delete

