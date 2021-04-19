CC = gcc
CFLAGS = -Wall -std=gnu11
LDLIBS = -lm

all:Out.htm
	cat $<

Out.htm:main
	./$< > $@

main:main.o adapt.o
main.o adapt.o:adapt.h

clean:
	find -type f -executable -delete
	$(RM) *.o [Oo]ut*
