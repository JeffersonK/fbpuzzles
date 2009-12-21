CC=gcc
CFLAGS=-Wall

all: facebull

facebull: facebull.o

test: facebull
	@python runtests.py

submit:	clean facebull
	@tar -cvf facebull.tar Makefile facebull.c
	
clean:
	rm -f facebull facebull.o
	rm -f *.pyc *.pyo
	rm -f test/*.run
	rm -f facebull.tar
	rm -f facebull.tar.gz
	rm -f facebull.tgz
	
clobber: clean
	rm -f *~

