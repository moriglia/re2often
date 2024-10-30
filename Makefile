F_SRC_ROOT = src/fortran
F_SRC_LDPC = $(F_SRC_ROOT)/ldpc

BUILD = build
O_FORTRAN = $(BUILD)/fortran
O_LDPC    = $(O_FORTRAN)/ldpc
O_MODS    = $(O_FORTRAN)/mods

OPTIMIZE ?= -O3


.PHONY: all test

all: lib/libldpc.a


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


test: test/fortran/test_ldpc_decoder_construction test/fortran/test_ldpc_edge_list


# Generic library rule
lib/lib%.a :
	ar r $@ $(O_FORTRAN)/$(*F)/*.o



# Generic executable
%: %.f90
	gfortran $< -o $@ -lldpc -Llib -I$(O_MODS)


clean:
	rm -rf build/
	rm lib/*.a
