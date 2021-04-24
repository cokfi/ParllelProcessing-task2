// title
#pragma warning(disable : 4996)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include <time.h>
#include <math.h>
int main(argc, argv)
int argc;
char** argv;
{
    int rank, size, slice_start, slice_end, n;
    int vx_sign, vy_sign;
    n = 1000;
    double gen_mat[6][1000], M, G, v_zero, light_y, t_step, r, x, y;//general matrix:(x,y,v_x,v_y,a_x,a_y)
    double temp_mat[6][1000];// for each processor calculations
    M = 2 * (pow(10, 30));//mass
    G = 6.674 * (pow(10, -17));//gravitational constant
    v_zero = 200;// velocity 200km/sec
    light_y = 9 * (pow(10, 12));//light year in kilometers
    t_step = 450000;
    //FILE* fp = fopen("result.txt", "w");
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    srand(time(NULL) + rank);
    slice_start = rank * n / size;//might not be equal
    //init
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 6; j++)
            temp_mat[j][i] = 0;
    }
    if (rank == size - 1 && n % size != 0)
        slice_end = n;
    else
        slice_end = (rank + 1) * n / size;

    for (int i = slice_start; i < slice_end; i++) {// generate starting possition and velocity
        gen_mat[0][i] = (double)(rand() / RAND_MAX) * 100.0 * light_y;
        //fprintf(fp, "%lf ", gen_mat[0][i]);
        gen_mat[1][i] = (double)(rand() / RAND_MAX) * 100.0 * light_y;
        vx_sign = 2 * rand() % 4 - 1; // sign generate
        vy_sign = 2 * rand() % 4 - 1; // sign generate
        gen_mat[2][i] = (double)((rand() / RAND_MAX) + 0.5) * v_zero * vx_sign;
        gen_mat[3][i] = (double)((rand() / RAND_MAX) + 0.5) * v_zero * vy_sign;
        gen_mat[4][i] = 0;
        gen_mat[5][i] = 0;
        temp_mat[0][i] = gen_mat[0][i];
        temp_mat[1][i] = gen_mat[1][i];
        temp_mat[2][i] = gen_mat[2][i];
        temp_mat[3][i] = gen_mat[3][i];
    }
    //fclose(fp);
    //if(rank ==0){
    //	for (int i = 0; i<1000;i++)
    //        	fprintf(fp,"%f ",gen_mat[0][i]);
    //}

    MPI_Allreduce(&temp_mat[0], &gen_mat[0], 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&temp_mat[1], &gen_mat[1], 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (int j = slice_start; j < slice_end; j++) {// generate starting acc
        for (int i = 0; i < n; i++) {
            x = gen_mat[0][i] - gen_mat[0][j];
            y = gen_mat[1][i] - gen_mat[1][j];
            r = sqrt(pow(x, 2) + pow(y, 2));
            if (r == 0 || i == j)
                continue;
            temp_mat[4][j] = G * M * x / pow(r, 3);
            temp_mat[5][j] = G * M * y / pow(r, 3);
        }
    }
    for (int step = 0; step < 100; step++) {
        for (int j = slice_start; j < slice_end; j++) {// calculate new position and new speed
            temp_mat[0][j] += temp_mat[2][j] * t_step; // x = x0+v*t
            temp_mat[1][j] += temp_mat[3][j] * t_step;// y = y0 +v*t
            temp_mat[2][j] += temp_mat[4][j] * t_step; // vx = vx0 +ax*t
            temp_mat[3][j] += temp_mat[5][j] * t_step;// vy = vy0 + ay*t
        }

        MPI_Allreduce(&temp_mat[0], &gen_mat[0], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&temp_mat[1], &gen_mat[1], 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        for (int j = slice_start; j < slice_end; j++) {// calculate new acceleration.
            for (int i = 0; i < n; i++) {
                x = gen_mat[0][i] - gen_mat[0][j];
                y = gen_mat[1][i] - gen_mat[1][j];
                r = sqrt(pow(x, 2) + pow(y, 2));
                if (r == 0 || i == j)
                    continue;
                temp_mat[4][j] = G * M * x / pow(r, 3); // a = f/m
                temp_mat[5][j] = G * M * y / pow(r, 3);
            }
        }
    }
    MPI_Finalize();
    //for (int i = 0; i<1000;i++)
    //	fprintf(fp,"%f ",gen_mat[0][i]);
    return(0);
}