test:
	
	gcc -g -Wall -I/home/gsl/include -c main.c functions.c read_fn.c 
	
	gcc -L/home/gsl/lib main.o functions.o read_fn.o  -lgsl -lgslcblas -lumfpack -lm -o test

clean:
	rm test main.o functions.o read_fn.o

#gcc -g main.c read_fn.c functions.c -o test -lm