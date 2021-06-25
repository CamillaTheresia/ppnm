CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

all: out.txt out.times.png

out.txt: main
	./$< 4 3 > $@

out.times.png: out.times.txt out.times.gsl.txt Makefile
	echo '\
	set term png ;\
	set title "QR-decomposition times as function of matrix size" ;\
	set xlabel "matrix size N" ;\
	set ylabel "time, sec" ;\
	set key left ;\
	set out "$@" ;\
	N=10 ;\
	f(x)=(x/N)**3 ;\
	fit f(x) "$<" via N ;\
	plot \
	 "$<" with points title "Measurement for own QR-decomp" \
	, f(x) with lines title "Fit gives $$(x/%3.0f)^3$$ for own QR-decomp"%(N)\
	, "out.times.gsl.txt" with points title "Measurement for gsl\_linalg\_QR\_decomp" ;\
	'| pyxplot

out.times.txt: main Makefile
	>$@
	for n in $$(seq 220 10 380);\
	do \time --format "$$n %U" --append --output $@ ./$< $$n ;\
	done

out.times.gsl.txt: partcgsl Makefile
	>$@
	for n in $$(seq 220 10 380);\
	do \time --format "$$n %U" --append --output $@ ./$< $$n ;\
	done
clean:
	$(RM) out* main *.o *.log partcgsl
