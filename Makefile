parallel:
	mpicxx main_2d_par.cpp -o main_2d_par.o
seq:
	mpicxx main_2d_seq.cpp -o main_2d_seq.o
clean:
	rm -f *.o