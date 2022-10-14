CC=gcc
# BLASDIR=../openblas
BLASDIR=../../Build/openblas/build/
IDIRS=-I$(BLASDIR)/generated
LBLAS=-L$(BLASDIR)/lib -Wl,-rpath,$(BLASDIR)/lib -lopenblas -lpthread
# LBLAS=-lcblas -lblas
# CFLAGS=-Wall -Wextra -g -fdiagnostics-color=always -ffast-math $(IDIRS) -fopenmp
CFLAGS=-ffast-math -O3 $(IDIRS)
# CFLAGS=-O3 -ffast-math -fopenmp $(IDIRS) 
LFLAGS=-llapacke -llapack -lm $(LBLAS)
ERRORLOG=logs/compiler_error.log

_DEPS = utils/utils.c ham_gen/ham_gen.c params/params.c io/io.c diag/diag.c
DEPS = $(patsubst %,src/%,$(_DEPS))

OBJ = $(patsubst %.c,build/%.o,$(_DEPS))

_EXECS = exact_diag_simulation calculate_dist_vs_gfuncsq \
calculate_imbalance output_hamiltonian gfunc_complexity_test gfunc_complexity_test_2
EXECS = $(patsubst %,build/%,$(_EXECS))

TESTS = $(wildcard tests/*.c)

COUP=0

default: $(EXECS)

all: $(EXECS) tests

bins: $(EXECS)

build/%.o: src/%.c
	$(CC) -c -o $@ $^ $(CFLAGS)

$(EXECS): build/%: src/%.c $(OBJ)
	$(CC) -o $@  $^ $(CFLAGS) $(LFLAGS)

build/run_diag_params_c$(COUP): build/%: src/run_diag_params.c $(OBJ)
	$(CC) -o $@  $^ $(CFLAGS) $(LFLAGS)

build/tests/%: tests/%.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LFLAGS)

tests: $(patsubst tests/%.c,build/tests/%,$(TESTS))

clean:
	rm -rf build/*/*.o build/exact_diag_simulation
	find data -size 0 -print -delete

