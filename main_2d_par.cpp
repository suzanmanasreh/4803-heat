#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <sstream>
#include <mpi.h>
#include <omp.h>
#include "boost/multi_array.hpp"

using namespace std;

typedef boost::multi_array<double, 2> arr_2d;
typedef arr_2d::index idx;

void print_x(arr_2d x, double time, int rank);
void print_1d(vector<double> x, double time, string name, int rank);
void output_to_file(arr_2d x, int step, int x_dim, int y_dim);
void gather_and_output(arr_2d &x, MPI_Comm comm, int rank, int p, int x_cells, int y_cells, int x_domains, int y_domains, int x_dim, int y_dim, int x_start, int x_end, int y_start, int y_end, double curr_time, int k);

template <typename T>
T square(T num) {
    return num * num;
}

int main(int argc, char *argv[]) {
    int rank, p;

    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);

    double t1, t2;
    double init_inner = 10.0;
    double init_border = 0;

    double curr_time = 0;
    double dt = .01;
    double Lx = 1.0;
    double Ly = 1.0;
    int x_dim = 600;
    int y_dim = 600;
    double dx = Lx / x_dim;
    double dy = Ly / y_dim;
    int num_steps = 100;
    // heat coeff
    double alpha = .1;

    int w_out = 200;

    int x_domains = 6;
    int y_domains = 6;

    int num_domains = x_domains * y_domains;

    if (p != num_domains) {
        if (rank == 0)
            printf("Aborting: Mismatch between number of processors %d and number of domains %d\n", p, num_domains);
        return 0;
    }
    if (x_dim % x_domains != 0 || y_dim % y_domains != 0) {
        if (rank == 0)
            printf("Aborting: Dimensions not divisible by domain size\n");
        return 0;
    }

    int x_cells = x_dim / x_domains;
    int y_cells = y_dim / y_domains;

    int x_my_total = x_cells + 2;
    int y_my_total = y_cells + 2;

    int x_global = x_dim + 2;
    int y_global = y_dim + 2;

    int x_total = x_global + (2 * x_domains);
    int y_total = y_global + (2 * y_domains);

    double min_h = dx;
    if (dy < dx) min_h = dy;

    double max_dt = (min_h * min_h) / (4 * alpha);

    if (dt >= max_dt) {
        if (rank == 0) {
            printf("The parameter dt is too large. Setting it to max_dt: %f, dx: %f, dy: %f, min_h: %f, alpha: %f\n", max_dt, dx, dy, min_h, alpha);
        }
        dt = max_dt;
    }

    int n_threads = omp_get_max_threads();
    if (rank == 0) {
        printf("Using %d OpenMP threads\n", n_threads);
    }


    MPI_Comm comm_cart;
    int dims[2] = {y_domains, x_domains};
    int periods[2] = {0, 0};

    MPI_Cart_create(comm, 2, dims, periods, 0, &comm_cart);

    int new_rank;
    // rank should not change in the new communicators
    MPI_Comm_rank(comm_cart, &new_rank);
 
    // get my coordinates in the new communicator
    int my_coords[2];
    MPI_Cart_coords(comm_cart, new_rank, 2, my_coords);

    // print my location in the 2D torus.
    // printf("[MPI process %d] I am located at (%d, %d).\n", new_rank, my_coords[0], my_coords[1]);

    // might have to use loops if this causes problems later
    int x_start = 1;
    int x_end = x_start + x_cells - 1;
    int y_start = 1;
    int y_end = y_start + y_cells - 1; 

    enum DIRECTIONS {DOWN, UP, LEFT, RIGHT};
    string neighbor_names[4] = {"down", "up", "left", "right"};
    vector<int> neighbor(4);
 
    // get left and right neighbours
    MPI_Cart_shift(comm_cart, 0, 1, &neighbor[LEFT], &neighbor[RIGHT]);
 
    // get up and down neighbours
    MPI_Cart_shift(comm_cart, 1, 1, &neighbor[DOWN], &neighbor[UP]);


    if (neighbor[DOWN] == MPI_PROC_NULL) {
        x_my_total += 1;
    }
    if (neighbor[LEFT] == MPI_PROC_NULL) {
        // the y-axis is in the x-direction
        // so we 1 to the start and end of y due to left shift
        y_start += 1;
        y_end += 1;

        y_my_total += 1;
    }
    if (neighbor[UP] == MPI_PROC_NULL) {
        x_start += 1;
        x_end += 1;

        x_my_total += 1;

    }
    if (neighbor[RIGHT] == MPI_PROC_NULL) {
        y_my_total += 1;
    }

    MPI_Datatype column;
    MPI_Type_vector(x_cells, 1, y_my_total, MPI_DOUBLE, &column);
    MPI_Type_commit(&column);

    arr_2d x(boost::extents[x_my_total][y_my_total]);
    arr_2d prev(boost::extents[x_my_total][y_my_total]);

    fill(x.data(), x.data() + x.num_elements(), init_inner);
    fill(prev.data(), prev.data() + prev.num_elements(), init_inner);


    // initilize borders + ghost cells on left and right border to border values
    if (neighbor[LEFT] == MPI_PROC_NULL) {
        #pragma omp parallel for
        for (idx i = 0; i < x_my_total; i++) {
            x[i][0] = init_border;
            x[i][1] = init_border;
        }
    }
    if (neighbor[RIGHT] == MPI_PROC_NULL) {
        #pragma omp parallel for
        for (idx i = 0; i < x_my_total; i++) {
            x[i][y_my_total - 1] = init_border;
            x[i][y_my_total - 2] = init_border;
        }
    }

    // initialize borders + ghost cells on up + down border to border values
    if (neighbor[UP] == MPI_PROC_NULL) {
        #pragma omp parallel for
        for (idx j = 0; j < y_my_total; j++) {
            x[0][j] = init_border;
            x[1][j] = init_border;
        }
    }
    if (neighbor[DOWN] == MPI_PROC_NULL) {
        #pragma omp parallel for
        for (idx j = 0; j < y_my_total; j++) {
            x[x_my_total - 1][j] = init_border;
            x[x_my_total - 2][j] = init_border;
        }
    }


    double sx = (alpha * dt) / (dx * dx);
    double sy = (alpha * dt) / (dy * dy);

    t1 = MPI_Wtime();
    int k = 0;
    gather_and_output(x, comm, rank, p, x_cells, y_cells, x_domains, y_domains, x_dim, y_dim, x_start, x_end, y_start, y_end, curr_time, k);

    // time stepper
    for (k = 1; k <= num_steps; k++) {
        double diff = 0;
        double global_diff = 0;

        #pragma omp parallel for reduction(+:diff)
        for (int i = x_start; i <= x_end; i++) {
            for (int j = y_start; j <= y_end; j++) {
                prev[i][j] = x[i][j];
                x[i][j] += (sx * (x[i + 1][j] - (2*x[i][j]) + x[i - 1][j])) + (sy * (x[i][j + 1] - (2*x[i][j]) + x[i][j - 1]));
                diff += square(x[i][j] - prev[i][j]);
            }
        }
        MPI_Status status;
        int tag = 1;

        // update up/down boundaries at ghost cells on receiver end
        MPI_Sendrecv(&x[x_start][y_start], y_cells, MPI_DOUBLE, neighbor[UP], tag, &x[x_end + 1][y_start], y_cells, MPI_DOUBLE, neighbor[DOWN], tag, comm_cart, &status);

        MPI_Sendrecv(&x[x_end][y_start], y_cells, MPI_DOUBLE, neighbor[DOWN], tag, &x[x_start - 1][y_start], y_cells, MPI_DOUBLE, neighbor[UP], tag, comm_cart, &status);

        tag = 2;
        // update left/right boundaries at ghost cells on receiver end
        MPI_Sendrecv(&x[x_start][y_end], 1, column, neighbor[RIGHT], tag, &x[x_start][y_start - 1], 1, column, neighbor[LEFT], tag, comm_cart, &status);

        MPI_Sendrecv(&x[x_start][y_start], 1, column, neighbor[LEFT], tag, &x[x_start][y_end + 1], 1, column, neighbor[RIGHT], tag, comm_cart, &status);
        curr_time += dt;

        if (k % w_out == 0)
            gather_and_output(x, comm, rank, p, x_cells, y_cells, x_domains, y_domains, x_dim, y_dim, x_start, x_end, y_start, y_end, curr_time, k);

        MPI_Allreduce(&diff, &global_diff, 1, MPI_DOUBLE, MPI_SUM, comm);
        global_diff = sqrt(global_diff);

        if (global_diff < .01) {
            printf("rank %d convergence at step %d\n", rank, k);
            break;
        }
    }

    t2 = MPI_Wtime();

    double tot_time = t2 - t1;
    if (rank == 0) printf("total time: %f\n", tot_time);

    MPI_Finalize();

    return 0;
}

// print our grid of values (debug function)
void print_x(arr_2d x, double time, int rank) {
    printf("rank: %d, curr_time: %f\n", rank, time);
    for(idx i = 0; i < x.size(); i++){
		for(idx j = 0; j < x[i].size(); j++) {
            printf("%5.02f ", x[i][j]);
        }
        printf("\n");
	}
}

// print 1d vector (debug function)
void print_1d(vector<double> x, double time, string name, int rank) {
    printf("curr_time: %f\n", time);
    printf("rank %d %s: ", rank, name.c_str());
    for(int i = 0; i < x.size(); i++){
        printf("%5.02f ", x[i]);
	}
    printf("\n");
}




void gather_and_output(arr_2d &x, MPI_Comm comm, int rank, int p, int x_cells, int y_cells, int x_domains, int y_domains, int x_dim, int y_dim, int x_start, int x_end, int y_start, int y_end, double curr_time, int k) {
    
    int block_size = x_cells * y_cells;

    vector<double> x_flat(block_size, 0);
    vector<double> x_final(x_dim * y_dim, 0);

    int j = 0;
    for (int i = x_end; i >= x_start; i--, j++) {
        for (int k = 0; k < y_cells; k++) {
            x_flat[j*y_cells + k] = x[i][y_start + k];
        }
    }

    MPI_Gather(x_flat.data(), block_size, MPI_DOUBLE, x_final.data(), block_size, MPI_DOUBLE, 0, comm);
    MPI_Barrier(comm);

    if (rank == 0) {
        arr_2d x_final_2d(boost::extents[x_dim][y_dim]);

        for (idx d = 0; d < p; d++) {
            long block_x_offset = ((d % x_domains) * x_cells);
            long block_y_offset = (d / x_domains) * y_cells;

            #pragma omp parallel for collapse(2)
            for (idx i = 0; i < x_cells; i++) {
                for (idx j = 0; j < y_cells; j++) {
                    long x_idx = block_x_offset + i;
                    long y_idx = block_y_offset + j;
                    long final_idx = d*block_size + i*y_cells + j;

                    x_final_2d[x_idx][y_idx] = x_final[final_idx];
                }
            }
        }

        output_to_file(x_final_2d, k, x_dim, y_dim);
    }
}

void output_to_file(arr_2d x, int step, int x_dim, int y_dim) {
    ostringstream file_name;
    file_name << "2dpomp/ex_";
    file_name << setw(4) << setfill('0') << step;
    file_name << ".vtk";

    double x_spacing = 1.0 / (x_dim - 1);
    double y_spacing = 1.0 / (y_dim - 1);

    FILE *file = fopen(file_name.str().c_str(), "w");
    fprintf(file, "# vtk DataFile Version 3.0\n");
    fprintf(file, "Temperature data\n");
    fprintf(file, "ASCII\n");
    fprintf(file, "DATASET STRUCTURED_POINTS\n");
    fprintf(file, "DIMENSIONS %d %d 1\n", y_dim, x_dim);
    fprintf(file, "ORIGIN 0 0 0\n");
    fprintf(file, "SPACING %f %f 1\n", y_spacing, x_spacing);
    fprintf(file, "POINT_DATA %d\n", x_dim * y_dim);
    fprintf(file, "SCALARS Temp float 1\n");
    fprintf(file, "LOOKUP_TABLE default\n");
    
    for(int i = 0; i < x_dim; i++){
		for(int j = 0; j < y_dim; j++) {
            fprintf(file, "%4.02f\n", x[i][j]);
        }
	}
}