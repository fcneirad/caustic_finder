#
#
LIBS = -lm -lgsl -lgslcblas -L/usr/lib64
CC = g++
CFLAGS = -fopenmp
#

all: causticfinder

causticfinder: 
	$(CC) -o causticfinder_init causticfinder_init.c $(LIBS)
	$(CC) -o causticfinder_cmap causticfinder_cmap.c $(LIBS) $(CFLAGS)

clean: 
	rm -f causticfinder_init causticfinder_cmap
