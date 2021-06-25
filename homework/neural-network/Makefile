CFLAGS += -Wall -std=gnu11 $$(gsl-config --cflags) -Ofast
LDLIBS += -lm $$(gsl-config --libs)

out.png:out.txt Makefile
	echo '\
	set term png;\
	set out "$@";\
	set title "Neural network curve fitting";\
	set xlabel "x";\
	set ylabel "y";\
	plot\
	 "$<" index 0 with points title "Datapoints"\
	,"$<" index 1 using 1:2 with lines title "Neural network fit"\
	,"$<" index 1 using 1:3 with lines title "Neural network fit derivative"\
	,"$<" index 1 using 1:4 with lines title "Neural network fit integral"\
	'|pyxplot

out.txt:main
	./$< 1>$@ 2>log

main:main.o ann.o qnewton.o

clean:
	$(RM) main *.o out* log*
