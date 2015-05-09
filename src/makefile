include ../makefile.in

PROGRAM=fvm

SRC=$(wildcard *.f)
OBJ=$(patsubst %.f,%.o,$(SRC))

.f.o:
	$(F77) $(CFLAGS) -c $*.f

$(PROGRAM) :$(OBJ)
	$(F77) -o $(PROGRAM) $(OBJ)

clean:;    rm -f $(OBJ) core
