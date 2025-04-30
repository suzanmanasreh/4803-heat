par:
	mpicxx -fopenmp main_2d_par.cpp -o main_2d_par.o
par_old:
	mpicxx main_2d_par_old.cpp -I/opt/homebrew/Cellar/boost/1.87.0_1/include -L/opt/homebrew/Cellar/boost/1.87.0_1/lib -o main_2d_par_old.o
seq:
	mpicxx main_2d_seq.cpp -o main_2d_seq.o
clean:
	rm -f *.o