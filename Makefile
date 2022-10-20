CC=gcc
# MKLROOT=/opt/intel/mkl
MKLROOT=/apps/intel_2018/mkl
LMKL=-Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
# BLASDIR=../openblas
# BLASDIR=../../Build/openblas/build/
# IDIRS=-I$(BLASDIR)/generated
# IDIRS=-I$(MKLROOT)/include
# LBLAS=-L$(BLASDIR)/lib -Wl,-rpath,$(BLASDIR)/lib -lopenblas -lpthread
# LBLAS=-lcblas -lblas
# CFLAGS=-Wall -Wextra -g -fdiagnostics-color=always -ffast-math $(IDIRS) -fopenmp -ftrapv -fwrapv
CFLAGS=-ffast-math -O3 $(IDIRS)
# CFLAGS=-O3 -ffast-math -fopenmp $(IDIRS) 
# LFLAGS=-llapacke -llapack -lm $(LBLAS)
# LFLAGS=$(LMKL)
LFLAGS=$(LMKL)
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

$(wildcard build/run_diag_params_c*): build/%: src/%.c $(OBJ)
	$(CC) -o $@  $^ $(CFLAGS) $(LFLAGS)

build/tests/%: tests/%.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LFLAGS)

tests: $(patsubst tests/%.c,build/tests/%,$(TESTS))

clean:
	rm -rf build/*/*.o build/exact_diag_simulation
	find data -size 0 -print -delete

