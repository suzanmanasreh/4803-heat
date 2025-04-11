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

using namespace std;

void print_x(vector<vector<double> > x, double time);
void output_to_file(vector<vector<double> > x, int step, int x_dim, int y_dim);

template <typename T>
T square(T num) {
    return num * num;
}

int main(int argc, char *argv[]) {
    double t1, t2;
    double init_inner = 10.0;
    double init_border = 0;

    double curr_time = 0;
    double dt = .01;
    double Lx = 1.0;
    double Ly = 1.0;
    int x_dim = 100;
    int y_dim = 100;
    double dx = Lx / x_dim;
    double dy = Ly / y_dim;
    int num_steps = 100;
    // heat coeff
    double alpha = .1;

    double w_out = 1;
    double writes = 0;

    int x_total = x_dim + 2;
    int y_total = y_dim + 2;

    double min_h = dx;
    if (dy < dx) min_h = dy;

    double max_dt = (min_h * min_h) / (4 * alpha);
    // printf("dx: %f, dy: %f, min_h: %f, max_dt: %f\n", dx, dy, min_h, max_dt);

    if (dt >= max_dt) {
        printf("dt too large. setting it to %f\n", max_dt);
        dt = max_dt;
    }

    vector<vector<double> > x(x_total, vector<double> (y_total, init_inner));
    vector<vector<double> > prev(x_total, vector<double> (y_total, init_inner));
    // printf("init_inner: %lf\n", init_inner);

    for (int i = 0; i < x_total; i++) {
        x[i][0] = init_border;
        x[i][y_total - 1] = init_border;
    }

    for (int j = 0; j < y_total; j++) {
        x[0][j] = init_border;
        x[x_total - 1][j] = init_border;
    }

    // print_x(x, curr_time);
    // output_to_file(x, 0, x_dim, y_dim);

    double sx = (alpha * dt) / (dx * dx);
    double sy = (alpha * dt) / (dy * dy);

    MPI_Init(&argc, &argv);
    int rank, p;
    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);

    t1 = MPI_Wtime();

    for (int k = 1; k <= num_steps; k++) {
        double diff = 0;
        for (int i = 1; i <= x_dim; i++) {
            for (int j = 1; j <= y_dim; j++) {
                prev[i][j] = x[i][j];
                x[i][j] += (sx * (x[i + 1][j] - (2*x[i][j]) + x[i - 1][j])) + (sy * (x[i][j + 1] - (2*x[i][j]) + x[i][j - 1]));
                diff += square(x[i][j] - prev[i][j]);
            }
        }
        curr_time += dt;
        // print_x(x, curr_time);
        output_to_file(x, k, x_dim, y_dim);
        if (diff < .01) {
            printf("convergence at step %d\n", k);
            break;
        }
    }

    // if (rank == 0) {
    //     print_x(x, curr_time);
    // }

    t2 = MPI_Wtime();

    double tot_time = t2 - t1;
    printf("total time: %f\n", tot_time);

    MPI_Finalize();
}

// print our grid of values
void print_x(vector<vector<double> > x, double time) {
    printf("curr_time: %f\n", time);
    for(int i = 0; i < x.size(); i++){
		for(int j = 0; j < x[i].size(); j++) {
            printf("%5.02f ", x[i][j]);
        }
        printf("\n");
	}
}

void output_to_file(vector<vector<double> > x, int step, int x_dim, int y_dim) {
    ostringstream file_name;
    file_name << "2d/ex_";
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
    for(int i = 1; i <= x_dim; i++){
		for(int j = 1; j <= y_dim; j++) {
            fprintf(file, "%4.02f\n", x[i][j]);
        }
	}
}