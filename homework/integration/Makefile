CC = gcc
CFLAGS = -Wall -std=gnu11 $$(gsl-config --cflags)
LDLIBS = -lm $$(gsl-config --libs)

all:Out.htm
	cat $<

Out.htm:main
	./$< > $@

main:main.o adapt.o
main.o adapt.o:adapt.h

clean:
	find -type f -executable -delete
	$(RM) *.o [Oo]ut*
