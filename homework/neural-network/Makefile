CFLAGS += -Wall -std=gnu11 $$(gsl-config --cflags) -Ofast
LDLIBS += -lm $$(gsl-config --libs)

out.png:out.txt Makefile
	echo '\
	set term png;\
	set out "$@";\
	set title "neural network curve fitting";\
	set xlabel "x";\
	set ylabel "y";\
	plot\
	 "$<" index 0 with points title "data"\
	,"$<" index 1 with lines title "fit"\
	'|pyxplot

out.txt:main
	./$< 1>$@ 2>log

main:main.o ann.o qnewton.o

clean:
	$(RM) main *.o out* log*
