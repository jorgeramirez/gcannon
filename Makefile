CC = gcc
BOMPDIR = bin/openmp
BMPIDIR = bin/mpi
UDIR = util
OMPDIR = openmp
CFLAGS= -fopenmp -lm

_DEPS = matrix.h matrix.c
DEPS = $(patsubst %, $(UDIR)/%, $(_DEPS))

all: cannon gcannon cannonmpi gcannonmpi

cannon: $(OMPDIR)/cannon.c $(DEPS)
	$(CC) -o $(BOMPDIR)/$@ $^ $(CFLAGS)

gcannon: $(OMPDIR)/gcannon.c $(DEPS)
	$(CC) -o $(BOMPDIR)/$@ $^ $(CFLAGS)

cannonmpi: mpi/cannon.c $(DEPS)
	mpicc $^ -o $(BMPIDIR)/cannonmpi -lm 


gcannonmpi: mpi/gcannon.c $(DEPS)
	mpicc $^ -o $(BMPIDIR)/gcannonmpi -lm 

clean:
	rm $(BOMPDIR)/* && rm $(BMPIDIR)/*
