CFLAGS += -Wall -std=c11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

all: out.txt plot.txt plot.png

plot.png: plot.txt data.txt Makefile
	echo '\
	set term png ;\
	set xlabel "Energy [GeV]" ;\
	set ylabel "cross-section" ;\
	plot \
	"data.txt" using 1:2:3 with yerrorbars title "data" \
	,"plot.txt" using 1:2 with lines title "fit with m=125.972 GeV/(c*c)"\
	' | pyxplot

plot.txt: main
	./$< 1> $@

out.txt: main
	./$< 2> $@

main: main.o qnewton.o

clean:
	$(RM) main *.o *.txt plot.png
