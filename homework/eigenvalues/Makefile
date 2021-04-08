CFLAGS = -Wall -O1 -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell /usr/bin/gsl-config --libs)

all: out.txt out.eigen.png

out.txt: main
	./$< > $@

out.eigen.png: out.txt Makefile
	echo '\
	set term png;\
	set title "The three lowest eigenfunctions" ;\
	set key left ;\
	set xlabel "x" ;\
	set ylabel "y" ;\
	set out "$@" ;\
	plot \
	"$<" index 2 with lines title "calculated eigenfunctions" \
	,"$<" index 3 using 1:2 with lines title "analytical eigenfunction 1" \
	,"$<" index 3 using 1:3 with lines title "analytical eigenfunction 2" \
	,"$<" index 3 using 1:4 with lines title "analytical eigenfunction 3" \
	' | pyxplot

clean:
	$(RM) out* main
