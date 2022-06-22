CC=gcc
CFLAGS=-Wall -g
LFLAGS=-llapacke -lm
ERRORLOG=logs/compiler_error.log

default: exactdiag

exactdiag: hamgen utils
	$(CC) $(CFLAGS) src/exact_diag_simulation.c -o build/exact_diag_simulation build/ham_gen/ham_gen.o build/utils/utils.o $(LFLAGS)| tee $(ERRORLOG)

hamgen: utils
	$(CC) $(CFLAGS) src/ham_gen/ham_gen.c -o build/ham_gen/ham_gen.o -c

utils: 
	$(CC) $(CFLAGS) src/utils/utils.c -o build/utils/utils.o -c

tests: utils hamgen
	$(CC) $(CFLAGS) tests/test_utils.c -o build/tests/test_utils build/utils/utils.o $(LFLAGS)

clean:
	rm -rf build/*/*.o build/exact_diag_simulation
