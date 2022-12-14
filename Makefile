CC=gcc
IDIRS=
LDIRS=

# -----------------
# MKL Configuration
# -----------------

# On cluster use these lines
# MKLROOT=/apps/intel_2018/mkl
# CMKL=
# LMKL=-Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group

# On local machine use these lines
# MKLROOT=/opt/intel/mkl
CMKL=-DMKL_ILP64 
CMKL += -m64
CMKL += -I"${MKLROOT}/include"
# LMKL += -L${MKLROOT}/lib/intel64
LMKL=
# Static
LMKL += -Wl,--start-group
LMKL += ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a
# LMKL += ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a
LMKL += ${MKLROOT}/lib/intel64/libmkl_intel_thread.a
LMKL += ${MKLROOT}/lib/intel64/libmkl_core.a
LMKL += ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_ilp64.a
LMKL += -Wl,--end-group
# Common lines
# LMKL += -lgomp
LMKL += -liomp5
LMKL += -lpthread
LMKL += -lm
LMKL += -ldl

# ----------------------
# OpenBLAS Configuration
# ----------------------

# On triality use this
# BLASDIR=../openblas
# On local machine use this
BLASDIR=../../Build/openblas/build/

CBLAS += -I$(BLASDIR)/generated
LBLAS=-L$(BLASDIR)/lib 
LBLAS += -Wl,-rpath,$(BLASDIR)/lib
LBLAS += -lopenblas
LBLAS += -lpthread

# -------------------------
# System BLAS Configuration
# -------------------------
LSYSBLAS=-lcblas -lblas

# -------------------
# Unity Configuration
# -------------------
TESTFLAGS=
TESTFLAGS += -DUNITY_INCLUDE_DOUBLE

# --------------------------------------
# Main Compiler and Linker Configuration
# --------------------------------------
CFLAGS=
LFLAGS=
CFLAGS += -Wall
CFLAGS += -Wextra
CFLAGS += -g
CFLAGS += -fdiagnostics-color=always
# CFLAGS += -ffast-math
# CFLAGS += -ftrapv -fwrapv
CFLAGS += -O2
# CFLAGS += -Iextern/unity
# CFLAGS += -fopenmp
CFLAGS += $(IDIRS)

# Use MKL
CFLAGS += $(CMKL)
LFLAGS += $(LMKL)

# Use OpenBLAS
# CFLAGS += $(CBLAS)
# LFLAGS += $(LBLAS)
# LFLAGS += -llapacke -llapack
# LFLAGS += -lm

# Use system BLAS
# LFLAGS += $(LSYSBLAS)
# LFLAGS += -llapacke -llapack
# LFLAGS += -lm


# -----------------
# Files to Compile 
# -----------------

# Function Files
_DEPS = utils/utils.c ham_gen/ham_gen.c params/params.c io/io.c diag/diag.c
DEPS = $(patsubst %,src/%,$(_DEPS))
OBJ = $(patsubst %.c,build/%.o,$(_DEPS))

# ERRORLOG=logs/compiler_error.log

# Executables
_EXECS = exact_diag_simulation calculate_dist_vs_gfuncsq \
calculate_imbalance output_hamiltonian sigma_exact_diag sigma_make_func \
exact_diag_batch batch_average
EXECS = $(patsubst %,build/%,$(_EXECS))

# External Dependencies
EXTDEPS = extern/unity/unity.c

# Tests
_SIMPLETESTS = test_utils test_ham_gen
_UNITYTESTS = test_gfuncsq test_utils2
SIMPLETESTS = $(patsubst %,build/tests/%,$(_SIMPLETESTS))
UNITYTESTS = $(patsubst %,build/tests/%,$(_UNITYTESTS))

default: $(EXECS)

all: $(EXECS) tests
bins: $(EXECS)
tests: $(UNITYTESTS) #$(SIMPLETESTS)

build/%.o: src/%.c
	$(CC) -c -o $@ $^ $(CFLAGS)

$(EXECS): build/%: src/%.c $(OBJ)
	$(CC) -o $@  $^ $(CFLAGS) $(LFLAGS)

# $(wildcard build/run_diag_params_c*): build/%: src/%.c $(OBJ)
# 	$(CC) -o $@  $^ $(CFLAGS) $(LFLAGS)

$(SIMPLETESTS): build/tests/%: tests/%.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LFLAGS)

$(UNITYTESTS): build/tests/%: tests/%.c $(OBJ) $(EXTDEPS)
	$(CC) -o $@ $^ $(TESTFLAGS) $(CFLAGS) $(LFLAGS)

# build/tests/test_gfuncsq: tests/test_gfuncsq.c $(OBJ) $(EXTDEPS)
# 	$(CC) -o $@ $^ $(TESTFLAGS) $(CFLAGS) $(LFLAGS)

# build/tests/%: tests/%.c $(OBJ)
# 	$(CC) -o $@ $^ $(CFLAGS) $(LFLAGS) 


clean:
	# rm -rf build/*/*.o $(EXECS)
	find build -type f -delete
	find data -size 0 -print -delete

