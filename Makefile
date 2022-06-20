CC=gcc
CFLAGS=-Wall -g
LFLAGS=-llapacke
ERRORLOG=logs/compiler_error.log

default: exactdiag

exactdiag: hamgen utils
	$(CC) $(CFLAGS) src/exact_diag_simulation.c -o build/exact_diag_simulation $(LFLAGS) build/ham_gen/ham_gen.o | tee $(ERRORLOG)

hamgen: 
	$(CC) $(CFLAGS) src/ham_gen/ham_gen.c -o build/ham_gen/ham_gen.o $(LFLAGS) -c

utils: 
	$(CC) $(CFLAGS) src/utils/utils.c -o build/utils/utils.o $(LFLAGS) -c

clean:
	rm -rf build/*/*.o build/exact_diag_simulation
