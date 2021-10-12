CC=gcc

default:./src/*.c ./include/*.h
	$(CC) -c ./src/*.c -I./include
	ar rv libcquadpack.a *.o
	rm *.o

clean:
	rm *.o *.a
