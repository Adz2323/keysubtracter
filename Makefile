.PHONY: default clean all

CC = gcc
CFLAGS = -O0 -g
INCLUDES = -I./xxhash
LIBS = -lgmp -pthread -lm

default: clean all

all: keysubtracter

keysubtracter: keysubtracter.o gmpecc.o util.o sha256.o base58.o rmd160.o bloom.o xxhash.o
	$(CC) $(CFLAGS) -o keysubtracter $^ $(LIBS)

keysubtracter.o: keysubtracter.c
	$(CC) $(CFLAGS) $(INCLUDES) -c keysubtracter.c

sha256.o: sha256.c
	$(CC) $(CFLAGS) -c sha256.c

base58.o: base58/base58.c
	$(CC) $(CFLAGS) -c base58/base58.c -o base58.o

rmd160.o: rmd160.c
	$(CC) $(CFLAGS) -c rmd160.c

gmpecc.o: gmpecc.c
	$(CC) $(CFLAGS) -c gmpecc.c

util.o: util.c
	$(CC) $(CFLAGS) -c util.c

bloom.o: bloom.c
	$(CC) $(CFLAGS) $(INCLUDES) -c bloom.c

xxhash.o: xxhash/xxhash.c
	$(CC) $(CFLAGS) -c xxhash/xxhash.c -o xxhash.o

clean:
	rm -f *.o keysubtracter
