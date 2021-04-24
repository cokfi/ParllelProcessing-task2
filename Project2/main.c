//C program that simulates N-body Problem
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
    int vx_sign, vy_sign, total_step, check;
    total_step = 1000;
    check = total_step / 2;
    n = 1000;
    double starttime, endtime;
    double M, G, v_zero, light_y, t_step, r, x, y;//general matrix:(x,y,v_x,v_y,a_x,a_y)
    double temp_mat[6][1000], gen_mat[6][1000];// for each processor calculations
    M = 2 * (pow(10, 30));//mass
    G = 6.674 * (pow(10, -17));//gravitational constant
    v_zero = 200;// velocity 200km/sec
    light_y = 9 * (pow(10, 12));//light year in kilometers
    t_step = 45000000;
    MPI_Init(&argc, &argv);
    starttime = MPI_Wtime();
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    srand(time(NULL) + rank);
    slice_start = rank * n / size;//might not be equal
    /*FILE *fpx_start;
    FILE *fpx_mid;
    FILE *fpx_final;
    FILE *fpy_start;
    FILE *fpy_mid;
    FILE *fpy_final;

    if(rank ==0){
        fpx_start = fopen("x_start.txt","w");
        fpx_mid = fopen("x_mid.txt","w");
        fpx_final = fopen("x_final.txt","w");
        fpy_start = fopen("y_start.txt","w");
        fpy_mid = fopen("y_mid.txt","w");
        fpy_final = fopen("y_final.txt","w");
    }
    */
    //init
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 6; j++) {
            temp_mat[j][i] = 0;
            gen_mat[j][i] = 0;
        }
    }
    if (rank == size - 1 && n % size != 0)
        slice_end = n;
    else
        slice_end = (rank + 1) * n / size;

    for (int i = slice_start; i < slice_end; i++) {// generate starting possition and velocity
        temp_mat[0][i] = (double)rand() / RAND_MAX * 100.0 * light_y;
        temp_mat[1][i] = (double)rand() / RAND_MAX * 100.0 * light_y;
        vx_sign = 2 * rand() % 4 - 1; // sign generate
        vy_sign = 2 * rand() % 4 - 1; // sign generate
        temp_mat[2][i] = ((double)rand() / RAND_MAX + 0.5) * v_zero * (double)vx_sign;
        temp_mat[3][i] = ((double)rand() / RAND_MAX + 0.5) * v_zero * (double)vy_sign;
        temp_mat[4][i] = 0;
        temp_mat[5][i] = 0;
    }

    MPI_Allreduce(&temp_mat[0], &gen_mat[0], n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&temp_mat[1], &gen_mat[1], n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    /*if (rank ==0){
    for (int i = 0; i< n;i++){
            fprintf(fpx_start,"%lf ",gen_mat[0][i]);
            fprintf(fpy_start,"%lf ",gen_mat[1][i]);
    }
    }
    */
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
    for (int step = 0; step < total_step; step++) {
        for (int j = slice_start; j < slice_end; j++) {// calculate new position and new speed
            temp_mat[0][j] += temp_mat[2][j] * t_step; // x = x0+v*t
            temp_mat[1][j] += temp_mat[3][j] * t_step;// y = y0 +v*t
            temp_mat[2][j] += temp_mat[4][j] * t_step; // vx = vx0 +ax*t
            temp_mat[3][j] += temp_mat[5][j] * t_step;// vy = vy0 + ay*t
        }
        MPI_Allreduce(&temp_mat[0], &gen_mat[0], n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&temp_mat[1], &gen_mat[1], n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        /*	if (step == check){
                if (rank ==0){
                    for (int i = 0; i< n;i++){
                            fprintf(fpx_mid,"%f ",gen_mat[0][i]);
                            fprintf(fpy_mid,"%f ",gen_mat[1][i]);
                    }
                }
            }*/
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
    if (rank == 0) {
        /*for (int i = 0; i<n;i++){
            fprintf(fpx_final,"%f ",gen_mat[0][i]);
            fprintf(fpy_final,"%f ",gen_mat[1][i]);
        }
            fclose(fpx_start);
            fclose(fpx_mid);
            fclose(fpx_final);
            fclose(fpy_start);
            fclose(fpy_mid);
            fclose(fpy_final);
        */
        endtime = MPI_Wtime();
        printf("That took %f seconds \n", endtime - starttime);
    }
    MPI_Finalize();
    return(0);
}

