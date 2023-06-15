PROGNAME=mcrscm

IDIR =./src/include
ODIR=./src/obj
SRCDIR=./src

CC=g++ -std=c++17
C=gcc

#CFLAGS=-I$(IDIR) `pkg-config --cflags gsl` `pkg-config --cflags RNAlib2` -O3 # for stuff
#LIBS=-lm `pkg-config --libs gsl` `pkg-config --libs RNAlib2` -fopenmp

CFLAGST=-I$(IDIR) `pkg-config --cflags gsl` -pthread -I/home/danielred/packages -I/home/danielred/packages/ViennaRNA -lboost_system -lboost_serialization -ggdb -fexceptions -Wall -pg # for testing
CFLAGS=-I$(IDIR) `pkg-config --cflags gsl` -O3 -pthread -I/home/danielred/packages -I/home/danielred/packages/ViennaRNA -lboost_system -lboost_serialization # for stuff with RNAfold 2.1.5

LIBS=-lm `pkg-config --libs gsl` -L/home/danielred/packages/ViennaRNA/lib -fno-lto -Wl,-fno-lto -lRNA -fopenmp -lgsl -lgslcblas -lpthread -lstdc++ -fopenmp # for RNAlib 2.1.5

_DEPS = randomgen.h dv_tools.h parameters.h rnarep.h annot.h cm.h rnarep_serialise.h cm_serialise.h broken.hpp 
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o dv_tools.o parameters.o rnarep.o annot.o randomgen.o cm.o broken.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

_OBJ_test = test.o annot.o parameters.o dv_tools.o bitmuveletek.o 
OBJ_test = $(patsubst %,$(ODIR)/%,$(_OBJ_test))

_OBJ_rev = reverse.o ca.o annot.o rnarep.o parameters.o randomgen.o dv_tools.o 
OBJ_rev = $(patsubst %,$(ODIR)/%,$(_OBJ_rev))

_OBJ_strgen = strgen2.o randomgen.o 
OBJ_strgen = $(patsubst %,$(ODIR)/%,$(_OBJ_strgen))

_OBJ_strdist = strdist.o randomgen.o annot.o rnarep.o parameters.o dv_tools.o 
OBJ_strdist = $(patsubst %,$(ODIR)/%,$(_OBJ_strdist))

_OBJ_randseq = randseq.o randomgen.o rnarep.o annot.o parameters.o dv_tools.o
OBJ_randseq = $(patsubst %,$(ODIR)/%,$(_OBJ_randseq))



$(ODIR)/%.o: $(SRCDIR)/%.cpp $(DEPS)
	@mkdir -p ${ODIR}
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/%.o: $(SRCDIR)/%.c $(DEPS)
	@mkdir -p ${ODIR}
	$(C) -c -o $@ $< $(CFLAGS)

$(PROGNAME): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: gdb
gdb: debug
gdb: CFLAGS=$(CFLAGST)
debug: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGST) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~  

.PHONY: run

run:
	./$(PROGNAME)  

test: $(OBJ_test)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
	#./test

rev: $(OBJ_rev)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

gen: $(OBJ_strgen)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

dist: $(OBJ_strdist)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

randseq: $(OBJ_randseq)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

