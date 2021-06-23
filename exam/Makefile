CFLAGS = -std=gnu1x -Wall -O2
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

all: out.txt

out.txt: main
	 ./$< > $@

main: main.o cholesky.o

clean:
	$(RM) main *.o out*
