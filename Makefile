FC := gfortran
FFLAGS := -O3 -march=native -flto -g -std=f2008 -Wall -Wextra -fopenmp
OBJS := m_config.o m_random.o m_mser.o
TESTS := dmc

.PHONY:	all clean

all: 	$(TESTS)

clean:
	$(RM) $(TESTS) $(OBJS) $(OBJS:.o=.mod)

# Dependency information
$(TESTS): $(OBJS)

# How to get .o object files from .f90 source files
%.o: %.f90
	$(FC) -c -o $@ $< $(FFLAGS)

# How to get executables from .o object files
%: %.o
	$(FC) -o $@ $^ $(FFLAGS)
