#@Author: B159973
#@Date:	10/4/2019
#@Course: Performance Programming - 2020
#@University of Edinburgh

SRC=MD.c control.c util.c 
OBJ=$(SRC:.c=.o)
CC=icc 
CFLAGS= -g  -mtune=broadwell  -xCore-AVX2 -qopenmp -qopt-dynamic-align

all: MD

MD: $(OBJ)
	$(CC) $(CFLAGS)  -o $@  $(OBJ) -lm -qopt-report3


output.dat: MD input.dat
	./MD

clean:
	rm -f MD $(OBJ) 

$(OBJ) : coord.h Makefile


