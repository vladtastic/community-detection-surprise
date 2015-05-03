
CC=g++
CFLAGS=-g -Wall -O2

SERIAL_SRC=comm-detec-serial.cpp Routines.cpp
SERIAL_OBJ=comm-detec-serial.o

SURPRISE_SRC=comm-detec-serial-surprise.cpp Routines.cpp
SURPRISE_OBJ=comm-detec-serial-surprise.o Routines.o
SURPRISE_OUT=surprise-serial.out

OMP_SRC=comm-detec-omp.cpp
OMP_OBJ=comm-detec-omp.o

GSL_INCLUDE=/nethome/vcoxall3/fa12/mga/project/gsl/include
GSL_LINK=/nethome/vcoxall3/fa12/mga/project/gsl/lib

all:
	echo "Nothing"

serial: ${SERIAL_OBJ}
	${CC} ${SERIAL_SRC} ${CFLAGS} -o serial.out

surprise: ${SURPRISE_OBJ}
	${CC} -L${GSL_LINK} ${SURPRISE_OBJ} -lgsl -lgslcblas -lm -o ${SURPRISE_OUT}

${SURPRISE_OBJ}:
	${CC} -I${GSL_INCLUDE} ${CFLAGS} -c ${SURPRISE_SRC}

clean:
	rm -rf *.o
