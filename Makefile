CC=gcc
CFLAGS=-Wall -Wextra -g -fdiagnostics-color=always -ffast-math -fopenmp
# CFLAGS=-O2 -ffast-math -fopenmp
LFLAGS=-llapacke -lm -lgsl -lcblas
ERRORLOG=logs/compiler_error.log

_DEPS = utils/utils.c ham_gen/ham_gen.c params/params.c
DEPS = $(patsubst %,src/%,$(_DEPS))

OBJ = $(patsubst %.c,build/%.o,$(_DEPS))

TESTS = $(wildcard tests/*.c)

default: exact_diag_simulation calculate_dist_vs_gfuncsq

build/%.o: src/%.c
	$(CC) -c -o $@ $^ $(CFLAGS)

exact_diag_simulation: $(OBJ)
	$(CC) -o build/$@ src/$@.c $^ $(CFLAGS) $(LFLAGS)

calculate_dist_vs_gfuncsq: $(OBJ)
	$(CC) -o build/$@ src/$@.c $^ $(CFLAGS) $(LFLAGS)

build/tests/%: tests/%.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LFLAGS)

tests: $(patsubst tests/%.c,build/tests/%,$(TESTS))

clean:
	rm -rf build/*/*.o build/exact_diag_simulation
	find data -size 0 -print -delete

