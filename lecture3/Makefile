CC = gcc
CFLAGS = -std=gnu99 -O -Wall
LDLIBS = -lm

out.txt: math variables
	./math > out.txt
	./variables >> out.txt

math.o: math.c
	$(CC) $(CFLAGS) -c $< -o $@

math: math.o
	$(CC) $(LDFLAGS) $< -o $@ $(LDLIBS)

variables.o: variables.c
	$(CC) $(CFLAGS) -c $< -o $@

variables: variables.o
	$(CC) $(LDFLAGS) $< -o $@ $(LDLIBS)

clean:
	$(RM) *.o out* math
