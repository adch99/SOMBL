CC=gcc
CFLAGS=-Wall -Wextra -g -fdiagnostics-color=always
# CFLAGS=-Wall -Wextra -g
# CFLAGS=-O2
LFLAGS=-llapacke -lm -lgsl -lcblas
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

