F_SRC_ROOT = src/fortran
F_SRC_LDPC = $(F_SRC_ROOT)/ldpc
F_TEST     = test/fortran
BUILD = build
O_FORTRAN = $(BUILD)/fortran
O_LDPC    = $(O_FORTRAN)/ldpc
O_MODS    = $(O_FORTRAN)/mods

OPTIMIZE ?= -O3

STDLIB = $(HOME)/.local/lib
# STDMOD = $(HOME)/.local/include/fortran_stdlib/GNU-11.4.0
INCLUDES = $(HOME)/.local/include


.PHONY: all testldpc testalpha test libraries py programs

all: libraries programs
libraries: lib/libldpc.a lib/libalpha.a lib/libsimtools.a lib/libnoisemapper.a
test: testldpc testalpha


# Missing directory
%/:
	mkdir -p $@


# LDPC Library
$(O_LDPC)/ldpc_decoder.o : $(O_LDPC)/ldpc_edge_list.o
$(O_LDPC)/%.o : $(F_SRC_LDPC)/%.f90
	@mkdir -p $(@D) $(O_MODS)
	gfortran -fPIC $(OPTIMIZE) -c $(F_SRC_LDPC)/$(*F).f90 -J$(O_MODS) -I$(O_MODS) -o $@


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


# Noise Mapper library
O_NM     = $(O_FORTRAN)/noisemapper
F_SRC_NM = $(F_SRC_ROOT)/noisemapper
$(O_NM)/%.o : $(F_SRC_NM)/%.f90
	@mkdir -p $(@D) $(O_MODS)
	gfortran $(OPTIMIZE) -c $(F_SRC_NM)/$(*F).f90 -J$(O_MODS) -I$(O_MODS) \
	-I$(INCLUDES) -L$(STDLIB) -lstdlib \
	-o $@

NM_MODS = noisemapper
NM_OBJS = $(patsubst %, $(O_NM)/%.o, $(NM_MODS))
lib/libnoisemapper.a: $(NM_OBJS) lib/libalpha.a


# Simtools library
O_SIMTOOLS     = $(O_FORTRAN)/simtools
F_SRC_SIMTOOLS = $(F_SRC_ROOT)/simtools
$(O_SIMTOOLS)/%.o : $(F_SRC_SIMTOOLS)/%.f90 lib/libldpc.a lib/libalpha.a lib/libnoisemapper.a
	@mkdir -p $(@D) $(O_MODS)
	caf $(OPTIMIZE) -c $(F_SRC_SIMTOOLS)/$(*F).f90 -J$(O_MODS) -I$(O_MODS) \
		-L$(STDLIB) -I$(INCLUDES) -lstdlib -lforbear \
		-fcoarray=lib -lcaf_openmpi -o $@

SIMTOOLS_MODS = sim_direct_channel sim_rrs
SIMTOOLS_OBJS = $(patsubst %, $(O_SIMTOOLS)/%.o, $(SIMTOOLS_MODS))
lib/libsimtools.a : $(SIMTOOLS_OBJS)




# Generic library rule
lib/lib%.a :
	ar r $@ $(O_FORTRAN)/$(*F)/*.o




# Programs
PROGS = reverse_soft_reconciliation direct_reconciliation try_mpi
F_SRC_PROGS = programs
O_PROGS = $(BUILD)/programs
programs: $(patsubst %, $(O_PROGS)/%, $(PROGS))

$(O_PROGS)/%: $(F_SRC_PROGS)/%.f90 lib/libsimtools.a
	@mkdir -p $(@D)
	caf $< -o $@ \
		-fcoarray=lib -lcaf_openmpi \
		-lsimtools -lFLAP -lIO-Fortran-Library -lforbear -lalpha -lldpc -lnoisemapper -lstdlib \
		-L$(STDLIB) -I$(INCLUDES) -Llib -I$(O_MODS)


# Generic executable
%: %.f90
	gfortran $< -o $@ -lldpc -lalpha -Llib -I$(O_MODS)





# Python module
CY_SRC_DIR = src/cython
PY_SO_SUFF = $(shell python3 -c "import distutils; print(distutils.sysconfig.get_config_var('EXT_SUFFIX'))")

re2often/simtools$(PY_SO_SUFF) : $(CY_SRC_DIR)/simtools.pyx lib/libsimtools.a

# Generic rule for python extension
re2often/%$(PY_SO_SUFF) :
	python3 setup.py build_ext -b . --only $(*F)


py : $(patsubst %,re2often/%$(PY_SO_SUFF), simtools)


clean:
	rm -rf build/
	rm -f $(wildcard lib/*.a)
	rm -f $(wildcard re2often/*$(PY_SO_SUFF))
	rm -rf re2often/__pycache__/

cleancycache :
	rm -rf $(wildcard build/temp.*-cpython-*/)
	rm -rf re2often/__pycache__/
