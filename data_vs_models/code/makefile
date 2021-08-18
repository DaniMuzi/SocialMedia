CC=gcc
OT=-ffast-math -O3 -Wall #-ggdb3
# OT=-ffast-math -O3 -Wall -Wextra -Werror -Wcast-align -Wcast-qual -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-include-dirs -Wredundant-decls -Wshadow -Wundef -Wno-unused -Wno-variadic-macros -Wno-parentheses -fdiagnostics-show-option # -pedantic -Wstrict-overflow=5
OP=$(OT)
LOP=-lm -o 




CHDR=./
CSRC=./
MSRC=./

MAIN=data_vs_models
VER=


EXE=$(MAIN)$(VER)$(SIZE)
TAG=$(EXE).out

OBJ=$(MAIN).o models_algorithm.o mt19937-64.o my_sort.o stats_functions.o basic_functions.o fitting_protocol.o
$(TAG)	: $(OBJ)
	$(CC) $(OBJ) $(LOP) $(TAG)
	rm *.o

$(MAIN).o :
	$(CC) $(OP) -I$(CHDR) -c $(MSRC)$(MAIN).c

models_algorithm.o  :
	$(CC) $(OP) -I$(CHDR) -c $(CSRC)models_algorithm.c

mt19937-64.o	:
	$(CC) $(OP) -I$(CHDR) -c $(CSRC)mt19937-64.c

my_sort.o  :
	$(CC) $(OP) -I$(CHDR) -c $(CSRC)my_sort.c

stats_functions.o  :
	$(CC) $(OP) -I$(CHDR) -c $(CSRC)stats_functions.c

basic_functions.o  :
	$(CC) $(OP) -I$(CHDR) -c $(CSRC)basic_functions.c

fitting_protocol.o  :
	$(CC) $(OP) -I$(CHDR) -c $(CSRC)fitting_protocol.c

clear	:
	rm data/*; rm index.dat; cp work/words.lst ./

clean	:
	rm *.o
