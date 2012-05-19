CC = gcc
BOMPDIR = bin/openmp
UDIR = util
OMPDIR = openmp
CFLAGS= -fopenmp

_DEPS = matrix.h matrix.c
DEPS = $(patsubst %, $(UDIR)/%, $(_DEPS))

$(BOMPDIR)/cannon: $(OMPDIR)/cannon.c $(DEPS)
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	rm $(BOMPDIR)/*
