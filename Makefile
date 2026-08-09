# Compiler selection: gfortran (default) or nvfortran
COMPILER ?= gfortran

ifeq ($(COMPILER),nvfortran)
  FC := nvfortran
  FFLAGS := -g -acc=gpu -fast -gpu=ccnative -Mpreprocess -Minfo=accel
else
  FC := gfortran
  FFLAGS := -O3 -march=native -g -std=f2008 -Wall -Wextra -fopenmp -cpp
endif

OBJS := m_config.o m_mser.o m_settings.o
TESTS := dmc

NDIM   ?= 3
NPART  ?= 1
POT    ?= ATOM
ATOMZ  ?= 1
FLOAT  ?= 64


CPPFLAGS = -DNDIM=$(NDIM) -DNPART=$(NPART) -DPOT_$(POT) -DATOMZ=$(ATOMZ) \
	-DFLOAT_$(FLOAT)

# Parameter string used to detect changes
PARAMS := COMPILER=$(COMPILER) NDIM=$(NDIM) NPART=$(NPART) POT=$(POT) \
	ATOMZ=$(ATOMZ) FLOAT=$(FLOAT)

.PHONY:	all clean settings force

all: 	settings $(TESTS)

settings:
	@echo "=== Build Settings ==="
	@echo "COMPILER = $(COMPILER)"
	@echo "FC       = $(FC)"
	@echo "FFLAGS   = $(FFLAGS)"
	@echo "CPPFLAGS = $(CPPFLAGS)"
	@echo "NDIM     = $(NDIM)"
	@echo "NPART    = $(NPART)"
	@echo "POT      = $(POT)"
	@echo "ATOMZ    = $(ATOMZ)"
	@echo "FLOAT    = $(FLOAT)"
	@echo "======================"

clean:
	$(RM) $(TESTS) $(OBJS) $(OBJS:.o=.mod) .params

# Rewrite .params only if the parameters changed
.params: force
	@echo '$(PARAMS)' | cmp -s - $@ || echo '$(PARAMS)' > $@

# Dependency information
$(TESTS): $(OBJS)

# Rebuild all objects when parameters change
$(OBJS): .params

# How to get .o object files from .f90 source files
%.o: %.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(CPPFLAGS)

# How to get executables from .o object files
%: %.o
	$(FC) -o $@ $^ $(FFLAGS)
