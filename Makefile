F_SRC_ROOT = src
BUILD ?= build
O_LIB = $(BUILD)/lib
O_INC = $(BUILD)/inc
O_LIBOBJ = $(BUILD)/libobj

O_PROGS = $(BUILD)/programs
S_PROGS = programs

OPTIMIZE ?= -O2
COMPILER ?= gfortran
CAF_COMP ?= caf


LIB = $(HOME)/.local/lib
INCLUDES = $(HOME)/.local/include


.PHONY: all library clean install programs

all: library
library: $(O_LIB)/libre2often.a


# Missing directory
%/:
	mkdir -p $@


# Library
$(O_INC)/%.mod $(O_LIBOBJ)/%.o : $(F_SRC_ROOT)/%.f90 $(O_LIBOBJ)/ $(O_INC)/
	$(COMPILER) $(OPTIMIZE) -c $< -J$(O_INC) -I$(O_INC) -I$(INCLUDES) -L$(LIB) -o $@


LIB_COMPONENTS = alphabet
LIB_OBJECTS    = $(patsubst %, $(O_LIBOBJ)/%.o, $(LIB_COMPONENTS))
$(O_LIB)/libre2often.a : $(LIB_OBJECTS) $(O_LIB)/
	ar r $@ $(LIB_OBJECTS)



# Programs
$(O_PROGS)/%: $(S_PROGS)/%.f90 $(O_PROGS)/ $(O_LIB)/libre2often.a
	$(CAF_COMP) $(OPTIMIZE) $< -o $@ \
		-J$(O_INC) -I$(O_INC) -I$(INCLUDES) -L$(LIB) -L$(O_LIB) \
		-lre2often -lfldpc -lfortran_stdlib -lIO-Fortran-Library -lforbear -lcaf_mpi

PROG_LIST = $(patsubst programs/%.f90, $(O_PROGS)/%, $(wildcard programs/*.f90) )
programs: $(PROG_LIST)


clean:
	rm -rf $(O_PROGS) $(O_LIB) $(O_INC) $(O_LIBOBJ)


cleanall:
	rm -rf build/


MODS = $(patsubst %, $(O_INC)/%.mod, $(LIB_COMPONENTS))
install: $(O_LIB)/libre2often.a $(MODS)
	cp $(O_LIB)/libre2often.a $(LIB)/
	cp $(O_INC)/*.mod $(INCLUDES)/
