FC := gfortran
FFLAGS := -O3 -march=native -g -std=f2008 -Wall -Wextra -fopenmp -cpp
OBJS := m_config.o m_mser.o m_settings.o m_physics.o
TESTS := dmc

NDIM   ?= 3
NPART  ?= 1
POT    ?= ATOM
ATOM_Z ?= 1
REGION ?= NONE

CPPFLAGS = -DNDIM=$(NDIM) -DNPART=$(NPART) -DPOT_$(POT) -DATOM_Z=$(ATOM_Z) \
	-DREGION_$(REGION)

# Parameter string used to detect changes
PARAMS := NDIM=$(NDIM) NPART=$(NPART) POT=$(POT) ATOM_Z=$(ATOM_Z) REGION=$(REGION)

.PHONY:	all clean settings force

all: 	settings $(TESTS)

settings:
	@echo "=== Build Settings ==="
	@echo "FC       = $(FC)"
	@echo "FFLAGS   = $(FFLAGS)"
	@echo "CPPFLAGS = $(CPPFLAGS)"
	@echo "NDIM     = $(NDIM)"
	@echo "NPART    = $(NPART)"
	@echo "POT      = $(POT)"
	@echo "ATOM_Z   = $(ATOM_Z)"
	@echo "REGION   = $(REGION)"
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

# Dependencies
m_physics.o: m_settings.o
