#include <vector>
#include <iostream>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>

using namespace std;

int main() {
    float init_inner = 10;
    float init_border = 0;

    float init_time = 0;
    float curr_time = init_time;
    float dt = .1;
    float L = 5;
    float dx = 1;
    int x_dim = (int) (L / dx);
    int num_steps = 21;
    // heat coeff
    float alpha = 1;
    float w_out = 1;
    float writes = 0;

    vector<float> x_vals(x_dim + 2, init_inner);
    x_vals[0] = init_border;
    x_vals[x_dim + 1] = init_border;

    ofstream file("data.csv");
    file << "x,";
    for (int j = 1; j < (x_dim + 1); j++) {
        file << j << ",";
    }
    file << endl;

    for (int i = 0; i < num_steps; i++) {
        if (abs(curr_time - (w_out * writes)) < .0001) {
            printf("%f,", curr_time);
            file << curr_time << ",";
        }
        for (int j = 1; j < (x_dim + 1); j++) {
            if (abs(curr_time - (w_out * writes)) < .0001) {
                file << x_vals[j];
                printf("%f, ", x_vals[j]);
                if (j != x_dim) {
                    file << ",";
                } else {
                    writes++;
                    printf("\n");
                    file << "\n";
                }
            }
            x_vals[j] += ((alpha * dt) / (dx*dx)) * (x_vals[j - 1] - (2 * x_vals[j]) + x_vals[j + 1]);
        }
        curr_time += dt;
    }
    file.close();
}