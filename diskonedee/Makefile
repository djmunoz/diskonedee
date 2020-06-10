EXEC=diskonedee

CC=gcc

SRC_DIR=src
INCLUDE_DIR=.
CFLAGS=-I$(INCLUDE_DIR) -I$(GSL_INCLUDE_DIR) -std=c99 -g


OBJECT_DIR=$(SRC_DIR)

LIBS=-lm -lgsl -lgslcblas

_DEPS = disk.h
DEPS = $(patsubst %,$(INCLUDE_DIR)/%,$(_DEPS))

_OBJS = main.o global.o invert.o invert_gsl.o grid.o disk.o advance.o io.o ics.o parameters.o 
OBJS = $(patsubst %,$(OBJECT_DIR)/%,$(_OBJS))

$(OBJECT_DIR)/%.o: %.c $(DEPS)
		   $(CC) -c -o $@ $< $(CFLAGS)

$(EXEC): $(OBJS)
	$(CC) $^ $(CFLAGS) $(LIBS) -o $@ 


clean:
	rm -f $(OBJECT_DIR)/*.o $(EXEC) *~ core $(INCLUDE_DIR)/*~ 





 




















