# Makefile for cdd-057.

# Select ANSI C compiler
#CC = /usr/local/bin/gcc
CC = /bin/cc

# Location of include files
INCLUDEDIR = .

# Compiler options
#CFLAGS = -g -O -I$(INCLUDEDIR)
#CFLAGS = -g -pg -O -I$(INCLUDEDIR)
CFLAGS = -O -I$(INCLUDEDIR)

########## You shouldn't have to change anything after this point ##########

cddio.o: cddio.c cdd.h cdddef.h dplex.h dplexdef.h
	$(CC) $(CFLAGS) -c cddio.c

cddarith.o: cddarith.c cdd.h cdddef.h dplex.h dplexdef.h
	$(CC) $(CFLAGS) -c cddarith.c

cdd.o: cdd.c cdd.h cdddef.h dplex.h dplexdef.h
	$(CC) $(CFLAGS) -c cdd.c

dplex.o: dplex.c dplex.h dplexdef.h
	$(CC) $(CFLAGS) -c dplex.c

dplex_test.o: dplex_test.c dplex.h dplexdef.h
	$(CC) $(CFLAGS) -c dplex_test.c

setoper.o: setoper.c
	$(CC) $(CFLAGS) -c setoper.c

cdd: cdd.o cddio.o cddarith.o dplex.o setoper.o
	$(CC) $(CFLAGS) cdd.o cddio.o cddarith.o dplex.o setoper.o -o cdd 

dplex_test: dplex.o dplex_test.o setoper.o
	$(CC) $(CFLAGS) dplex.o dplex_test.o setoper.o -o dplex_test

clean:
	rm -rf core a.out cdd dplex_test *.o *~

all: cdd dplex_test

