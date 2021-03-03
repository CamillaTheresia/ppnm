CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

.PHONEY: all
all: out.txt

out.txt: main
	./$< > out.txt

.PHONEY: clean
clean:
	$(RM) main.o main out.txt

test:
	@echo $(CFLAGS)
	@echo $(LDLIBS)
