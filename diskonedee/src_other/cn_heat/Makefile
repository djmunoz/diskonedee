EXEC=cnmethod

CC=gcc

INCLUDE_DIR=.
GSL_INCLUDE_DIR=/usr/include/gsl
CFLAGS=-I$(INCLUDE_DIR) -I$(GSL_INCLUDE_DIR)

OBJECT_DIR=.

LIBS=-lm

_DEPS = cnmethod.h
DEPS = $(patsubst %,$(INCLUDE_DIR)/%,$(_DEPS))

_OBJS = cnmethod.o
OBJS = $(patsubst %,$(OBJECT_DIR)/%,$(_OBJS))

$(OBJECT_DIR)/%.o: %.c $(DEPS)
		   $(CC) -c -o $@ $< $(CFLAGS)

$(EXEC): $(OBJS)
	$(CC) $^ $(CFLAGS) $(LIBS) -o $@ 


clean:
	rm -f $(OBJECT_DIR)/*.o $(EXEC) *~ core $(INCLUDE_DIR)/*~ 





 




















