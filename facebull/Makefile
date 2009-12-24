CC=gcc
CFLAGS=-Wall -g

all: facebull megbb pathwalk edgewalk

facebull: facebull.o

edgewalk: edgewalk.o

pathwalk: pathwalk.o

megbb: megbb.o

test: facebull
	@python runtests.py

submit:	clean facebull
	@tar -cvf facebull.tar Makefile facebull.c
	
clean:
	rm -f *.pyc *.pyo
	rm -f test/*.run
	rm -f facebull.tar
	rm -f facebull.tar.gz
	rm -f facebull.tgz
	rm -f *.o
	rm -f edgewalk pathwalk megbb facebull

clobber: clean
	rm -f *~

