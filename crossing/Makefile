#
# Target all builds everything.
#
#
CSDP_ROOT=$(HOME)/src/misc/Csdp
export CFLAGS=-m64 -march=native -mtune=native -Ofast -fopenmp -ansi -Wall -DBIT64 -DUSEOPENMP -DSETNUMTHREADS -DUSESIGTERM -DUSEGETTIME -I$(CSDP_ROOT)/include
#
# LIBS settings for 64 bit Linux/unix systems.
#
export LIBS=-static -L$(CSDP_ROOT)/lib -lsdp -llapack -lblas -lm -lgfortran -lquadmath -lm
all: example
#
# This builds the example code.
#
example: example.o 
	$(CC) $(CFLAGS) example.o $(LIBS) -o example
#
# To clean up the directory.
#
clean:
	rm -f *.o
	rm -f example
	rm -f prob.dat-s
	rm -f prob.sol



