CC = gcc
CFLAGS = -std=gnu99 -lgmp -Wall

all: demo

demo: ec_gmp_p_mul.o
	$(CC) $(CFLAGS) ec_gmp_p_mul.c -o demo

.PHONY: clean
clean:
	rm -rf *o demo