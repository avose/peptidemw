# Author:  Aaron D. Vose
# License: This is software is released into the public domain by the author in March 2013.

OPTS= -O3 -Wall


all: peptide

#

fastautil.o: fastautil.c fastautil.h types.h Makefile
	gcc ${OPTS} -c fastautil.c

peptide.o: peptide.c peptide.h fastautil.h types.h Makefile
	gcc ${OPTS} -c peptide.c

peptide: peptide.o fastautil.o Makefile
	gcc ${OPTS} -o peptide peptide.o fastautil.o 

#

clean:
	rm -f *.o peptide

strip: clean
	rm -f *~
