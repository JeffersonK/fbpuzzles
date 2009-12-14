CC=gcc
CFLAGS=-Wall -g

all: facebull

facebull: facebull.o

test: facebull
	@python runtests.py

submit:	clean facebull
	@tar -cvf facebull.tar Makefile facebull.c
	
clean:
	rm -f facebull facebull.o
	rm -f *.pyc *.pyo
	rm -f facebull_tests/*.run
	rm -f facebull.tar
	
clobber: clean
	rm -f *~

