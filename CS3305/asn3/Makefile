GCC = gcc

default: clean asn3

asn3: asn3.c
	$(GCC) $^ -o asn3.out -lpthread

clean:
	rm -rf *.out