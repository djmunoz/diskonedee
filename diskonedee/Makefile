EXEC=diskonedee

CC=gcc

SRC_DIR=src
INCLUDE_DIR=.
GSL_INCLUDE_DIR=/usr/include/gsl
CFLAGS=-I$(INCLUDE_DIR) -I$(GSL_INCLUDE_DIR) -std=c99


OBJECT_DIR=$(SRC_DIR)

LIBS=-lm -lgsl -lgslcblas

_DEPS = disk.h
DEPS = $(patsubst %,$(INCLUDE_DIR)/%,$(_DEPS))

_OBJS = main.o global.o invert.o grid.o 
OBJS = $(patsubst %,$(OBJECT_DIR)/%,$(_OBJS))

$(OBJECT_DIR)/%.o: %.c $(DEPS)
		   $(CC) -c -o $@ $< $(CFLAGS)

$(EXEC): $(OBJS)
	$(CC) $^ $(CFLAGS) $(LIBS) -o $@ 


clean:
	rm -f $(OBJECT_DIR)/*.o $(EXEC) *~ core $(INCLUDE_DIR)/*~ 





 




















