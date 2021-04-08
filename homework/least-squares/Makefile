CFLAGS = -Wall -O1 -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell /usr/bin/gsl-config --libs)

all: out.txt out.example.png out.halftime.png out.unc.png

out.txt: main
	./$< > $@

out.example.png: out.txt Makefile
	echo '\
	set term png;\
	set title "Check that the least-square fit works";\
	set key left ;\
	set xlabel "x" ;\
	set ylabel "y" ;\
	set out "$@" ;\
	plot \
	 "$<" index 2 using 1:2:3 with errorbars title "data" \
	,"$<" index 2 using 1:4 with lines title "fit function" \
	' | pyxplot

out.halftime.png: out.txt Makefile
	echo '\
	set term png;\
	set title "Decay of 224Ra";\
	set key right ;\
	set xlabel "x" ;\
	set ylabel "y" ;\
	set out "$@" ;\
	plot \
	 "$<" index 4 using 1:2:3 with errorbars title "data" \
	,"$<" index 4 using 1:4 with lines title "fit function" \
	' | pyxplot

out.unc.png: out.txt Makefile
	echo '\
	set term png;\
	set title "Decay of 224Ra with uncertainties on ck";\
	set key right ;\
	set xlabel "x" ;\
	set ylabel "y" ;\
	set out "$@" ;\
	plot \
	 "$<" index 7 using 1:2:3 with errorbars title "data" \
	,"$<" index 7 using 1:4 with lines title "fit function" \
	,"$<" index 7 using 1:5 with lines title "fit with plus uncertainty on ck" \
	,"$<" index 7 using 1:6 with lines title "fit with minus uncertainty on ck" \
	' | pyxplot

clean:
	$(RM) out* main
