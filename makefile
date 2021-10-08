PROGNAME=simulation

IDIR =./src/include
ODIR=./src/obj
SRCDIR=./src

CC=g++
C=gcc

#CFLAGS=-I$(IDIR) `pkg-config --cflags gsl` `pkg-config --cflags RNAlib2` -ggdb -fexceptions -Wall -pg # for testing
CFLAGS=-I$(IDIR) `pkg-config --cflags gsl` `pkg-config --cflags RNAlib2` -O3 # for stuff

LIBS=-lm `pkg-config --libs gsl` `pkg-config --libs RNAlib2` -fopenmp

_DEPS = ca.h randomgen.h dv_tools.h parameters.h rnarep.h annot.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o ca.o dv_tools.o parameters.o rnarep.o annot.o randomgen.o
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
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/%.o: $(SRCDIR)/%.c $(DEPS)
	$(C) -c -o $@ $< $(CFLAGS)

$(PROGNAME): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~  

.PHONY: run

run:
	./$(PROGNAME)  

test: $(OBJ_test)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
	./test

rev: $(OBJ_rev)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

gen: $(OBJ_strgen)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

dist: $(OBJ_strdist)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

randseq: $(OBJ_randseq)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

