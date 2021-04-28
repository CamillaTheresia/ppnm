CFLAGS += -Wall -std=c11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

all: Out.txt

Out.txt: main
	./$< > $@

main: main.o newton.o gs.o

clean:
	rm -f main *.o [Oo]ut*
