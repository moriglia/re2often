F_SRC_ROOT = src/fortran
F_SRC_LDPC = $(F_SRC_ROOT)/ldpc
F_TEST     = test/fortran
BUILD = build
O_FORTRAN = $(BUILD)/fortran
O_LDPC    = $(O_FORTRAN)/ldpc
O_MODS    = $(O_FORTRAN)/mods

OPTIMIZE ?= -O3

STDLIB = $(HOME)/.local/lib
STDMOD = $(HOME)/.local/include/fortran_stdlib/GNU-11.4.0


.PHONY: all testldpc testalpha test libraries

all: libraries test
libraries: lib/libldpc.a lib/libalpha.a lib/libsimtools.a
test: testldpc testalpha


# Missing directory
%/:
	mkdir -p $@


# LDPC Library
$(O_LDPC)/ldpc_decoder.o : $(O_LDPC)/ldpc_edge_list.o
$(O_LDPC)/%.o : $(F_SRC_LDPC)/%.f90
	@mkdir -p $(@D) $(O_MODS)
	gfortran $(OPTIMIZE) -c $(F_SRC_LDPC)/$(*F).f90 -J$(O_MODS) -I$(O_MODS) -o $@


LDPC_MODS = ldpc_edge_list ldpc_decoder
LDPC_OBJS = $(patsubst %, $(O_LDPC)/%.o, $(LDPC_MODS))
lib/libldpc.a : $(LDPC_OBJS)

TEST_LDPC = $(patsubst %, $(F_TEST)/test_ldpc_%, decoder_construction edge_list)
$(TEST_LDPC): lib/libldpc.a
testldpc: $(TEST_LDPC)



# Alphabet library
O_ALPHA     = $(O_FORTRAN)/alpha
F_SRC_ALPHA = $(F_SRC_ROOT)/alpha
$(O_ALPHA)/%.o : $(F_SRC_ALPHA)/%.f90
	@mkdir -p $(@D) $(O_MODS)
	gfortran $(OPTIMIZE) -c $(F_SRC_ALPHA)/$(*F).f90 -J$(O_MODS) -I$(O_MODS) -o $@

ALPHA_MODS = alpha_pam
ALPHA_OBJS = $(patsubst %, $(O_ALPHA)/%.o, $(ALPHA_MODS))
lib/libalpha.a: $(ALPHA_OBJS)

TEST_ALPHA = $(patsubst %, $(F_TEST)/test_alpha_%, pam)
$(TEST_ALPHA) : lib/libalpha.a
testalpha: $(TEST_ALPHA)


# Simtools library
O_SIMTOOLS     = $(O_FORTRAN)/simtools
F_SRC_SIMTOOLS = $(F_SRC_ROOT)/simtools
$(O_SIMTOOLS)/%.o : $(F_SRC_SIMTOOLS)/%.f90
	@mkdir -p $(@D) $(O_MODS)
	gfortran $(OPTIMIZE) -c $(F_SRC_SIMTOOLS)/$(*F).f90 -J$(O_MODS) -I$(O_MODS) \
		-L$(STDLIB) -I$(STDMOD) -lfortran_stdlib \
		-fcoarray=lib -lopenmpi -o $@

SIMTOOLS_MODS = sim_direct_channel
SIMTOOLS_OBJS = $(patsubst %, $(O_SIMTOOLS)/%.o, $(SIMTOOLS_MODS))
lib/libsimtools.a : $(SIMTOOLS_OBJS)




# Generic library rule
lib/lib%.a :
	ar r $@ $(O_FORTRAN)/$(*F)/*.o



# Generic executable
%: %.f90
	gfortran $< -o $@ -lldpc -lalpha -Llib -I$(O_MODS)


clean:
	rm -rf build/
	rm lib/*.a
