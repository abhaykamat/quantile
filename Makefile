

CPLUS = g++
SHELL=/bin/sh
COPT = -g0 -O2 -I.
LIBS = -lmpfr -ldl

all:test

rng.o : rng.cpp
	${CPLUS} ${COPT}  -c rng.cpp -o rng.o

test.o : test.cpp
	${CPLUS} ${COPT}  -c test.cpp -o test.o

test:test.o rng.o
	${CPLUS} ${COPT} -o test test.o rng.o ${LIBS}

clean:
	-rm -r -f *.o test

