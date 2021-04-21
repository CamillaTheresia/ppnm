CFLAGS = -Wall -std=gnu11 -Ofast
LDLIBS = -lm

all: Out.txt Out.error.png Out.samples.png

Out.error.png:out.data.txt Makefile
	echo '\
set term png ;\
set out "$@" ;\
set log y ;\
set title "Monte Carlo integration of x*x+y*y less than R*R?1:0" ;\
set xlabel "number of points N" ;\
set ylabel "actual error" ;\
set xtics 50000 ;\
set mxtics 5 ;\
a=1; b=-0.5; f(x)=a*x**b ;\
c=1; d=-1.0; g(x)=c*x**d ;\
fit f(x) "$<" using 1:2 via a,b ;\
fit g(x) "$<" using 1:3 via c,d ;\
plot \
 "$<" using 1:2 with points title "pseudo-random" \
,"$<" using 1:3 with points title "quasi-random" \
, f(x) with lines title "%3.1f*pow(n,{%3.1f})"%(a,b) \
, g(x) with lines title "%3.1f*pow(n,{%3.1f})"%(c,d) \
	'|pyxplot

ns=$(shell seq 1000 1250 200000)
out.data.txt:main
	>$@
	@for n in $(ns); do echo n=$$n; ./main $$n >>$@; done

Out.txt: main
	./$< > $@

main: main.o plainmc.o quasimc.c

Out.samples.png:out.data Makefile
	echo '\
set term png ;\
set out "$@" ;\
set title "Stratified sampling of x*x+y*y less than R*R?1:0" ;\
set xlabel "x" ;\
set ylabel "y" ;\
set size square ;\
plot \
 "$<" with points notitle \
	'|pyxplot

out.data:mains
	./$< 1>Log 2>$@

mains:mains.o strata.o

clean:
	$(RM) main mains *.o [Oo]ut* [Ll]og*
