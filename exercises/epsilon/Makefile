CC = gcc -fwrapv
CFLAGS = -std=gnu11 -O -Wall
LDLIBS = -lm

out.txt: epsilon equal
	./epsilon > out.txt
	./equal >> out.txt

epsilon: epsilon.o
	$(CC) $(LDFLAGS) $< -o $@ $(LDLIBS)

equal: equal.o
	$(CC) $(LDFLAGS) $< -o $@ $(LDLIBS)

epsilon.o: epsilon.c
	$(CC) $(CFLAGS) -c epsilon.c -o epsilon.o

equal.o: equal.c
	$(CC) $(CFLAGS) -c equal.c -o equal.o

clean:
	$(RM) *.o out* epsilon equal
