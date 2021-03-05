CFLAGS = -Wall -O1 -std=gnu11
LDLIBS = -lm

ms.pdf: ms.tex fig.pdf
	pdflatex $<
	pdflatex $<

out.data.txt: main
	./$< > $@

main: main.o exponential.o

fig.pdf: fig.gpi out.data.txt
	pyxplot $<

clean:
	$(RM) exponential.o fig.pdf main.o ms.aux ms.log out.data.txt main
