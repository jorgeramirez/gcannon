CC = gcc
BOMPDIR = bin/openmp
UDIR = util
OMPDIR = openmp
CFLAGS= -fopenmp -lm

_DEPS = matrix.h matrix.c
DEPS = $(patsubst %, $(UDIR)/%, $(_DEPS))

all: cannon gcannon

cannon: $(OMPDIR)/cannon.c $(DEPS)
	$(CC) -o $(BOMPDIR)/$@ $^ $(CFLAGS)

gcannon: $(OMPDIR)/gcannon.c $(DEPS)
	$(CC) -o $(BOMPDIR)/$@ $^ $(CFLAGS)

clean:
	rm $(BOMPDIR)/*
