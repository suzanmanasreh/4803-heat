#include <vector>
#include <iostream>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <cmath>

using namespace std;

void print_x(vector<vector<float>> x, float time);

int main() {
    float init_inner = 10;
    float init_border = 0;

    float curr_time = 0;
    float dt = .1;
    float Lx = 1;
    float Ly = 1;
    int x_dim = 4;
    int y_dim = 4;
    float dx = Lx / x_dim;
    float dy = Ly / y_dim;
    int num_steps = 11;
    // heat coeff
    float alpha = .1;

    float w_out = 1;
    float writes = 0;

    int x_total = x_dim + 2;
    int y_total = y_dim + 2;

    vector<vector<float>> x(x_total, vector<float> (y_total, init_inner));

    for (int i = 0; i < x_total; i++) {
        x[i][0] = init_border;
        x[i][y_total - 1] = init_border;
    }

    for (int j = 0; j < y_total; j++) {
        x[0][j] = init_border;
        x[x_total - 1][j] = init_border;
    }

    print_x(x, curr_time);

    float sx = (alpha * dt) / (dx * dx);
    float sy = (alpha * dt) / (dy * dy);

    for (int k = 0; k < num_steps; k++) {
        for (int i = 1; i <= x_dim; i++) {
            for (int j = 1; j <= y_dim; j++) {
                x[i][j] += (sx * (x[i + 1][j] - (2*x[i][j]) + x[i - 1][j])) + (sy * (x[i][j + 1] - (2*x[i][j]) + x[i][j - 1]));
            }
        }
        curr_time += dt;
        print_x(x, curr_time);
    }

    // for (int i = 0; i < num_steps; i++) {
    //     for (int j = 1; j < (x_dim + 1); j++) {
    //         // actual update
    //         x_vals[j] += ((alpha * dt) / (dx*dx)) * (x_vals[j - 1] - (2 * x_vals[j]) + x_vals[j + 1]);
    //     }
    //     curr_time += dt;
    // }
}

// print our grid of values
void print_x(vector<vector<float>> x, float time) {
    printf("curr_time: %f\n", time);
    for(int i = 0; i < x.size(); i++){
		for(int j = 0; j < x[i].size(); j++) {
            printf("%5.02f ", x[i][j]);
        }
		cout<<endl;
	}
}