CC=gcc
CFLAGS=-O3 -Wall -Wextra -g -finline-functions -unroll-all-loop

all: test

test: test.o conversion.o mask2n.o maskq.o

test.o: src/test.c include/conversion.h include/mask2n.h include/maskq.h
	$(CC) $(CFLAGS) -Iinclude -c src/test.c

conversion.o: src/conversion.c include/conversion.h include/mask2n.h include/maskq.h
	$(CC) $(CFLAGS) -Iinclude -c src/conversion.c

mask2n.o: src/mask2n.c include/mask2n.h
	$(CC) $(CFLAGS) -Iinclude -c src/mask2n.c

maskq.o: src/maskq.c include/maskq.h
	$(CC) $(CFLAGS) -Iinclude -c src/maskq.c

clean:
	rm -rf *.o test 