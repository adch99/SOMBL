CC=gcc
# CFLAGS=-Wall -Wextra -g -fdiagnostics-color=always -ffast-math -fopenmp -pg -O2
# TAU_MAKEFILE=Makefile.tau
# CC=tau_cc.sh
# CFLAGS=-Wall -Wextra -g -fdiagnostics-color=always \
# 	-ffast-math -pg -pedantic -W \
# 	-Wmissing-prototypes -Wstrict-prototypes \
# 	-Wconversion -Wshadow -Wpointer-arith -Wcast-qual \
# 	-Wcast-align -Wwrite-strings -Wnested-externs \
# 	-fshort-enums -fno-common -Dinline= -g -O2 -fopenmp
# CFLAGS=-Wall -Wextra -g -pg  -fopenmp
CFLAGS=-O2 -ffast-math -fopenmp
LFLAGS=-llapacke -lm -lgsl -lblas
# LFLAGS=-llapacke -lm -lgsl -lcblas -pg -lgcov
ERRORLOG=logs/compiler_error.log

default: exactdiag

exactdiag: hamgen utils
	$(CC) $(CFLAGS) src/exact_diag_simulation.c -o build/exact_diag_simulation build/ham_gen/ham_gen.o build/utils/utils.o $(LFLAGS)| tee $(ERRORLOG)

exactdiag1d: hamgen utils
	$(CC) $(CFLAGS) src/1d_exact_diag_simulation.c -o build/1d_exact_diag_simulation build/ham_gen/ham_gen.o build/utils/utils.o $(LFLAGS)| tee $(ERRORLOG)

hamgen: utils
	$(CC) $(CFLAGS) src/ham_gen/ham_gen.c -o build/ham_gen/ham_gen.o -c

utils: 
	$(CC) $(CFLAGS) src/utils/utils.c -o build/utils/utils.o -c

tests: utils hamgen
	$(CC) $(CFLAGS) tests/test_utils.c -o build/tests/test_utils build/utils/utils.o $(LFLAGS)
	$(CC) $(CFLAGS) tests/test_ham_gen.c -o build/tests/test_ham_gen build/utils/utils.o build/ham_gen/ham_gen.o $(LFLAGS)
	$(CC) $(CFLAGS) tests/check_matrix.c -o build/tests/check_matrix build/utils/utils.o build/ham_gen/ham_gen.o $(LFLAGS)
clean:
	rm -rf build/*/*.o build/exact_diag_simulation

