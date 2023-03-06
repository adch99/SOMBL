CC=gcc
IDIRS=
LDIRS=

# -----------------
# MKL Configuration
# -----------------
CMKL=
LMKL=


# On cluster use these lines
# On flock
# On leap
# MKLROOT=/apps/intel/oneapi/mkl/latest
# On local machine, though it isn't necessary
# MKLROOT=/opt/intel/mkl

ONCLUSTER=
ONLOCAL=yes
FLOCK=
LEAP=
ZEAL=

# For all clusters
ifdef ONCLUSTER
ifdef FLOCK
MKLROOT=/apps/intel_2018/mkl
endif
ifdef LEAP
MKLROOT=/apps/intel/oneapi/mkl/latest
CMKL += -I"$(MKLROOT)/include"
endif
ifdef ZEAL
endif
LMKL += -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a
LMKL += $(MKLROOT)/lib/intel64/libmkl_intel_thread.a
LMKL += $(MKLROOT)/lib/intel64/libmkl_core.a
LMKL += -Wl,--end-group
endif

# On local machine use these lines
# We are using static linking
ifdef ONLOCAL
CMKL += -DMKL_ILP64
CMKL += -m64
LMKL += -Wl,--start-group
LMKL += ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a
LMKL += ${MKLROOT}/lib/intel64/libmkl_intel_thread.a
LMKL += ${MKLROOT}/lib/intel64/libmkl_core.a
LMKL += ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_ilp64.a
LMKL += -Wl,--end-group
endif

# Common lines
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
UNITYROOT=extern/unity
TESTFLAGS=
TESTFLAGS += -DUNITY_INCLUDE_FLOAT
TESTFLAGS += -DUNITY_INCLUDE_DOUBLE
TESTFLAGS += -I$(UNITYROOT)
# TESTFLAGS += -DUNITY_INCLUDE_CONFIG_H

# --------------------------------------
# Main Compiler and Linker Configuration
# --------------------------------------
CFLAGS=
LFLAGS=
CFLAGS += -Wall
CFLAGS += -Wextra
CFLAGS += -g
CFLAGS += -fdiagnostics-color=always
CFLAGS += -ffast-math
CFLAGS += -ftrapv -fwrapv
CFLAGS += -O2
# CFLAGS += -Ofast
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
_DEPS = utils/utils.c ham_gen/ham_gen.c params/params.c io/io.c diag/diag.c gfunc/gfunc.c
DEPS = $(patsubst %,src/%,$(_DEPS))
OBJ = $(patsubst %.c,build/%.o,$(_DEPS))

# ERRORLOG=logs/compiler_error.log

# Executables
# _EXECS = exact_diag_simulation calculate_dist_vs_gfuncsq \
# calculate_imbalance output_hamiltonian sigma_exact_diag sigma_make_func \
# exact_diag_batch batch_average keldysh_window_batch check_io
_EXECS = keldysh_window_batch_nobins convert_txt_to_bin binary_reader \
keldysh_energy_batch_average_error keldysh_densities_error \
keldysh_window_batch_nobins_error
EXECS = $(patsubst %,build/%,$(_EXECS))

# External Dependencies
EXTDEPS = extern/unity/unity.c

# Tests
_SIMPLETESTS = test_utils test_ham_gen
_UNITYTESTS = test_gfuncsq test_utils2 test_io
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

$(UNITYTESTS): build/tests/%: tests/%.c $(EXTDEPS) $(OBJ) 
	$(CC) -o $@ $(TESTFLAGS) $(CFLAGS) $^ $(LFLAGS)

# build/tests/test_gfuncsq: tests/test_gfuncsq.c $(OBJ) $(EXTDEPS)
# 	$(CC) -o $@ $^ $(TESTFLAGS) $(CFLAGS) $(LFLAGS)

# build/tests/%: tests/%.c $(OBJ)
# 	$(CC) -o $@ $^ $(CFLAGS) $(LFLAGS) 


clean:
	# rm -rf build/*/*.o $(EXECS)
	find build -type f -delete
	find data -size 0 -print -delete

