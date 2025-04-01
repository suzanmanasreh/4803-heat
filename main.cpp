#include <vector>
#include <iostream>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>

using namespace std;

int main() {
    string vtk_str = "# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET RECTILINEAR_GRID\nDIMENSIONS 6 2 1\nX_COORDINATES 6 float\n0 1 2 3 4 5\nY_COORDINATES 2 float\n0 1\nZ_COORDINATES 1 float\n0\nCELL_DATA 5\nSCALARS temp float 1\nLOOKUP_TABLE default\n";
    float init_inner = 10;
    float init_border = 0;

    float init_time = 0;
    float curr_time = init_time;
    float dt = .1;
    float L = 5;
    float dx = .2;
    int x_dim = (int) (L / dx);
    int num_steps = 20;
    // heat coeff
    float alpha = 1;

    vector<float> x_vals(x_dim + 2, init_inner);
    x_vals[0] = init_border;
    x_vals[x_dim + 1] = init_border;

    for (int i = 0; i < num_steps; i++) {
        char str[20]; // = format("examples/ex_%04d.vtk", i);
        sprintf(str, "examples/ex_%04d.vtk", i);
        printf("%f,", curr_time);
        ofstream file(str);
        file << vtk_str;
        for (int j = 1; j < (x_dim + 1); j++) {
            file << x_vals[j] << " ";
            printf("%f, ", x_vals[j]);
            x_vals[j] += ((alpha * dt) / (dx*dx)) * (x_vals[j - 1] - (2 * x_vals[j]) + x_vals[j + 1]);
        }
        printf("\n");
        curr_time += dt;
        // printf("7\n");
        file.close();
    }
}