CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

all: erf.pyxplot.png gamma.pyxplot.png

erf.pyxplot.png: data.txt Makefile tabulated.erf.data.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Error function";\
	plot \
	 "$<" using 1:2 with line title "erf from math.h"\
	,"$<" using 1:3 with line title "gsl-sf-erf"\
	,"$<" using 1:4 with line title "homemade Erf"\
	,"tabulated.erf.data.txt" using 1:2 with points title "tabulated data"\
	' | tee test.erf.gpi | pyxplot

gamma.pyxplot.png: datag.txt Makefile tabulated.gamma.data.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Gamma function";\
	plot \
	 "$<" using 1:2 with line title "gamma from math.h"\
	,"$<" using 1:3 with line title "gsl-sf-gamma"\
	,"$<" using 1:4 with line title "homemade gamma"\
	,"tabulated.gamma.data.txt" using 1:2 with points title "tabulated data"\
	' | tee test.gamma.gpi |  pyxplot

test.erf.gpi:erf.pyxplot.png

test.gamma.gpi: gamma.pyxplot.png

data.txt: main
	./$< > $@

datag.txt: maing
	./$< > $@

main: main.o erf.o

maing: maing.o gamma.o

.PHONEY:clean
clean:
	$(RM) *.o main maing *.png data.txt datag.txt test*
test:
	echo $(CFLAGS)
	echo $(LDLIBS)
