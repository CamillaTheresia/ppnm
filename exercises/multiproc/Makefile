CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += -pthread
LDLIBS = -lm

out.txt: main
	./$< > $@

clean:
	$(RM) out* main
