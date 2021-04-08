CFLAGS = -Wall
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

all: out.txt out.u.png out.Tc.png

out.u.png: out.txt Makefile
	echo '\
	set term png;\
	set title "d2u=-u, u(0)=0, du(0)=1" ; \
	set out "$@" ;\
	plot \
	"$<" index 0 using 1:2 with lines title "u" \
	,"$<" index 0 using 1:3 with lines title "du" \
	,1 notitle \
	,-1 notitle \
	' | pyxplot

out.Tc.png: out.txt Makefile
	echo '\
	set term png;\
	set title "SIR model" ; \
	set out "$@" ;\
	plot \
	"$<" index 1 using 1:2 with lines title "S, Tc=2" \
	,"$<" index 1 using 1:3 with lines title "I, Tc=2" \
	,"$<" index 1 using 1:4 with lines title "R, Tc=2" \
	,"$<" index 2 using 1:2 with lines title "S, Tc=4" \
	,"$<" index 2 using 1:3 with lines title "I, Tc=4" \
	,"$<" index 2 using 1:4 with lines title "R, Tc=4" \
	,"$<" index 3 using 1:2 with lines title "S, Tc=8" \
	,"$<" index 3 using 1:3 with lines title "I, Tc=8" \
	,"$<" index 3 using 1:4 with lines title "R, Tc=8" \
	' | pyxplot

out.txt: main
	./$< 2> $@

clean:
	$(RM) out* main
