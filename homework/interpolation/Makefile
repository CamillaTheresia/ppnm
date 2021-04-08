CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default
default: linear.plot.png linear.integ.plot.png out.quadrac.txt quadratic.plot.png

linear.plot.png: out.linear.txt out.linear.xy.txt Makefile
	echo '\
	set terminal png; \
	set output "$@"; \
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Linear interpolation"; \
	plot \
	"out.linear.xy.txt" using 1:2 with line title "x and y used in interpolation"\
	,"out.linear.txt" using 1:2 with line title "S(z)"\
	,"out.linear.txt" using 1:4 with line title "gsl.interp.linear"\
	' | pyxplot

linear.integ.plot.png: out.linear.txt out.linear.xy.txt Makefile
	echo '\
	set terminal png; \
	set output "$@"; \
	set key top left;\
	set tics out;\
	set xlabel "x"; \
	set ylabel "y"; \
	set title "Linear interpolation integrated"; \
	plot \
	"out.linear.xy.txt" using 1:3 with line title "y integrated as it should be"\
	,"out.linear.txt" using 1:3 with line title "S(z) integrated"\
	,"out.linear.txt" using 1:5 with line title "gsl.integ.linear"\
	' | pyxplot

quadratic.plot.png: out.quadrac.txt Makefile
	echo '\
	set term png;\
	set title "Quadratic spline, y=x*x" ;\
	set key left ;\
	set xlabel "x" ;\
	set ylabel "y" ;\
	set out "$@" ;\
	plot \
	"$<" index 0 using 1:4 with lines title "y=x*x" \
	,"$<" index 1 using 1:4 with lines title "qspline y=x*x" \
	,"$<" index 2 using 1:4 with lines title "diff y" \
	,"$<" index 3 using 1:4 with lines title "qspline diff y" \
	,"$<" index 4 using 1:4 with lines title "integ y" \
	,"$<" index 5 using 1:4 with lines title "qspline integ y" \
	' | pyxplot

out.linear.txt: linear
	./$< > $@

out.quadrac.txt: quadrac
	./$< > $@

.PHONEY:clean
clean:
	$(RM) *txt linear *png quadrac



