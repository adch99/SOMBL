CC=gcc
# BLASDIR=../openblas
BLASDIR=../../Build/openblas/build/
IDIRS=-I$(BLASDIR)/generated
LBLAS=-L$(BLASDIR)/lib -Wl,-rpath,$(BLASDIR)/lib -lopenblas -lpthread
# LBLAS=-lcblas -lblas
CFLAG=-Wall
CFLAGS += -Wextra
CFLAGS += -g
CFLAGS += -fdiagnostics-color=always
# CFLAGS += -ffast-math
# CFLAGS += -ftrapv -fwrapv
CFLAGS += -O2
# CFLAGS += -Iextern/unity
# CFLAGS += -fopenmp
# CFLAGS += $(IDIRS)

TESTFLAGS=
TESTFLAGS += -DUNITY_INCLUDE_DOUBLE

LFLAGS=
LFLAGS += -llapacke -llapack
LFLAGS += -lm
LFLAGS += $(LBLAS)
ERRORLOG=logs/compiler_error.log

_DEPS = utils/utils.c ham_gen/ham_gen.c params/params.c io/io.c diag/diag.c
DEPS = $(patsubst %,src/%,$(_DEPS))

OBJ = $(patsubst %.c,build/%.o,$(_DEPS))

EXTDEPS = extern/unity/unity.c

_EXECS = exact_diag_simulation calculate_dist_vs_gfuncsq \
calculate_imbalance output_hamiltonian gfunc_complexity_test gfunc_complexity_test_2 check_sigma_gfuncsq
EXECS = $(patsubst %,build/%,$(_EXECS))

# TESTS = $(wildcard tests/*.c)
TESTS = tests/test_gfuncsq.c

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

build/tests/test_gfuncsq: tests/test_gfuncsq.c $(OBJ) $(EXTDEPS)
	$(CC) -o $@ $^ $(TESTFLAGS) $(CFLAGS) $(LFLAGS)

build/tests/%: tests/%.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LFLAGS) 

tests: $(patsubst tests/%.c,build/tests/%,$(TESTS))

clean:
	rm -rf build/*/*.o build/exact_diag_simulation
	find data -size 0 -print -delete

