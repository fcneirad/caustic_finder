#
#
LIBS = -lm -lgsl -lgslcblas -L/usr/lib64
CC = g++
#

all: causticfinder

causticfinder: causticfinder.c
	$(CC) -o causticfinder causticfinder.c $(LIBS)

clean: 
	rm -f *.o
