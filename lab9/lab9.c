// gcc -o lab9 lab9.c -lgsl -lgslcblas -lm

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

int j(int l, int nx) {
    return l / (nx + 1);
}

int i(int l, int nx) {
    return l - j(l, nx) * (nx + 1);
}

int l(int i, int j, double nx) {
    return i + j * (nx + 1);
}

void dyfuzja_macierzowo(int dt_value) {
    int nx = 40;
    int ny = 40;
    int N = (nx + 1) * (ny + 1);
    double delta = 1.0;
    double dt = 1.0;
    double Ta = 40.0;
    double Tb = 0.0;
    double Tc = 30.0;
    double Td = 0.0;

    double kb = 0.1;
    double kd = 0.6;
    int it_max = 2000;
    FILE* f;

    gsl_matrix *A = gsl_matrix_calloc(N, N);
    gsl_matrix *B = gsl_matrix_calloc(N, N);
    gsl_vector *c = gsl_vector_calloc(N);
    gsl_vector *T = gsl_vector_calloc(N);
    gsl_vector *T_old = gsl_vector_calloc(N);
    gsl_vector *d = gsl_vector_calloc(N);
    gsl_vector *dT = gsl_vector_calloc(N);

    int l_tmp = 0;

    // Macierze A B i vector c
    for(int i = 0; i <= nx; i++) {
        for(int j = 0; j <= ny; j++) {
            l_tmp = l(i, j, nx);

            // Wnętrze obszaru
            if(1 <= i && i <= nx - 1 && 1 <= j && j <= ny - 1) {
                gsl_matrix_set(A, l_tmp, l_tmp - nx - 1, dt / (2 * delta * delta));
                gsl_matrix_set(A, l_tmp, l_tmp - 1, dt / (2 * delta * delta));
                gsl_matrix_set(A, l_tmp, l_tmp + 1, dt / (2 * delta * delta));
                gsl_matrix_set(A, l_tmp, l_tmp + nx + 1, dt / (2 * delta * delta));
                gsl_matrix_set(A, l_tmp, l_tmp, -2 * dt / (delta * delta) - 1);

                gsl_matrix_set(B, l_tmp, l_tmp - nx - 1, -dt / (2 * delta * delta));
                gsl_matrix_set(B, l_tmp, l_tmp - 1, -dt / (2 * delta * delta));
                gsl_matrix_set(B, l_tmp, l_tmp + 1, -dt / (2 * delta * delta));
                gsl_matrix_set(B, l_tmp, l_tmp + nx + 1, -dt / (2 * delta * delta));
                gsl_matrix_set(B, l_tmp, l_tmp, 2 * dt / (delta * delta) - 1);
            }

            // WB Dirichleta
            if(i == 0) {
                gsl_matrix_set(A, l_tmp, l_tmp, 1);
                gsl_matrix_set(B, l_tmp, l_tmp, 1);
                gsl_vector_set(c, l_tmp, 0);
            }

            if(i == nx) {
                gsl_matrix_set(A, l_tmp, l_tmp, 1);
                gsl_matrix_set(B, l_tmp, l_tmp, 1);
                gsl_vector_set(c, l_tmp, 0);
            }

            // WB von Neumanna na górnym brzegu dla chwili n + 1
            if (1 <= i && i <= nx - 1 && j == ny) {
                gsl_matrix_set(A, l_tmp, l_tmp - nx - 1, -1 / (kb * delta));
                gsl_matrix_set(A, l_tmp, l_tmp, 1 + 1 / (kb * delta));
                for(int k = 0; k < N; k++)
                    gsl_matrix_set(B, l_tmp, k, 0);
                gsl_vector_set(c, l_tmp, Tb);
            }

            // WB von Neumanna na dolnym brzegu dla chwili n + 1
            if (1 <= i && i <= nx - 1 && j == 0)
            {
                gsl_matrix_set(A, l_tmp, l_tmp + nx + 1, -1 / (kd * delta));
                gsl_matrix_set(A, l_tmp, l_tmp, 1 + 1 / (kd * delta));
                for(int k = 0; k < N; k++)
                    gsl_matrix_set(B, l_tmp, k, 0);
                gsl_vector_set(c, l_tmp, Td);
            }
        }
    }

    // Vector t
    for (int i = 0; i <= nx; i++) {
        for (int j = 0; j <= ny; j++) {
            l_tmp = l(i, j, nx);

            if (i == 0)
                gsl_vector_set(T, l_tmp, Ta);
            else if (i == nx)
                gsl_vector_set(T, l_tmp, Tc);
            else
                gsl_vector_set(T, l_tmp, 0);
        }
    }

    int s;
    gsl_permutation *p = gsl_permutation_calloc(N);
    gsl_linalg_LU_decomp(A, p, &s);


    for (int it = 1; it <= it_max; it++) {
        // d = B * T;
        gsl_blas_dgemv(CblasNoTrans, 1.0, B, T, 0, d);
        // d = d + c;
        gsl_blas_daxpy(1.0, c, d);

        for(int k = 0; k < N; k++)
            gsl_vector_set(T_old, k, gsl_vector_get(T, k));

        // A * T = d
        gsl_linalg_LU_solve(A, p, d, T);

    
        if(it == 100) {
            if(dt_value == 1) {
                f = fopen("7_100.txt", "w");
                gsl_vector_fprintf(f, T, "%f");
            }
            else {
                f = fopen("8_100.txt", "w");
                for(int k = 0; k < N; k++)
                    gsl_vector_set(dT, k, (gsl_vector_get(T, k)-gsl_vector_get(T_old, k))/dt);
                gsl_vector_fprintf(f, dT, "%f");
            }
            fclose(f);
        }
        if(it == 200) {
            if(dt_value == 1) {
                f = fopen("7_200.txt", "w");
                gsl_vector_fprintf(f, T, "%f");
            }
            else {
                f = fopen("8_200.txt", "w");
                for(int k = 0; k < N; k++)
                    gsl_vector_set(dT, k, (gsl_vector_get(T, k)-gsl_vector_get(T_old, k))/dt);
                gsl_vector_fprintf(f, dT, "%f");
            }
            fclose(f);
        }
        if(it == 500) {
            if(dt_value == 1) {
                f = fopen("7_500.txt", "w");
                gsl_vector_fprintf(f, T, "%f");
            }
            else {
                f = fopen("8_500.txt", "w");
                for(int k = 0; k < N; k++)
                    gsl_vector_set(dT, k, (gsl_vector_get(T, k)-gsl_vector_get(T_old, k))/dt);
                gsl_vector_fprintf(f, dT, "%f");
            }
            fclose(f);
        }
        if(it == 1000) {
            if(dt_value == 1) {
                f = fopen("7_1000.txt", "w");
                gsl_vector_fprintf(f, T, "%f");
            }
            else {
                f = fopen("8_1000.txt", "w");
                for(int k = 0; k < N; k++)
                    gsl_vector_set(dT, k, (gsl_vector_get(T, k)-gsl_vector_get(T_old, k))/dt);
                gsl_vector_fprintf(f, dT, "%e");
            }
            fclose(f);
        }
        if(it == 2000) {
            if(dt_value == 1) {
                f = fopen("7_2000.txt", "w");
                gsl_vector_fprintf(f, T, "%f");
            }
            else {
                f = fopen("8_2000.txt", "w");
                for(int k = 0; k < N; k++)
                    gsl_vector_set(dT, k, (gsl_vector_get(T, k)-gsl_vector_get(T_old, k))/dt);
                gsl_vector_fprintf(f, dT, "%e");
            }
            fclose(f);
        }
    }

    gsl_permutation_free(p);
    gsl_matrix_free(A);
    gsl_matrix_free(B);
    gsl_vector_free(c);
    gsl_vector_free(T);
    gsl_vector_free(T_old);
    gsl_vector_free(d);
}

int main()
{
    // dyfuzja_macierzowo(1);
    dyfuzja_macierzowo(0);

    return 0;
}
