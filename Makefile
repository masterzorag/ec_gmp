CC=		gcc
CFLAGS=	-std=gnu99 -lgmp -Wall

all: demo

test: ec_gmp_p_mul.o
	$(CC) $(CFLAGS) demo.o -o demo

.PHONY: clean
clean:
	rm -rf *o demo