CFLAGS = -O -Wall -std=gnu11
LDLIBS = -lm

all: out.cmdline.txt out.stdin.txt out.file.txt

out.cmdline.txt: cmdline input.txt
	./cmdline $$(cat input.txt) > $@

out.stdin.txt: stdin input.txt
	./stdin < input.txt 1> $@

out.file.txt: file input.txt
	./file input.txt $@

input.txt: Makefile
	cat /dev/null >$@
	for x in $$(seq 1 1 5); do echo "$$x" >> $@; done

.PHONEY:clean
clean:
	$(RM) input.txt out.cmdline.txt cmdline out.stdin.txt stdin out.file.txt file
