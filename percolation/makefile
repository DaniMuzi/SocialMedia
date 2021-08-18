CC=gcc
OT=-ffast-math -O3 -g -Wall
OP=$(OT)
LOP=-lm -o

CHDR=./
CSRC=./
MSRC=./

MAIN=TS_phase_diag

TAG=$(MAIN).out

OBJ=$(MAIN).o percolation_algorithm.o my_sort.o 

$(TAG)	: $(OBJ)
	$(CC) $(OBJ) $(LOP) $(TAG)
	rm *.o

$(MAIN).o :
	$(CC) $(OP) -I$(CHDR) -c $(MSRC)$(MAIN).c

percolation_algorithm.o	:
	$(CC) $(OP) -I$(CHDR) -c $(CSRC)percolation_algorithm.c

my_sort.o	:
	$(CC) $(OP) -I$(CHDR) -c $(CSRC)my_sort.c

clear	:
	rm data/*; rm index.dat; cp work/words.lst ./

clean	:
	rm *.o
