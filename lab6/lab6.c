// gcc -o lab6 lab6.c mgmres.c -lm


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mgmres.h"


double j(int l, int nx) {
    return l/(nx+1);
}

double i(int l, int nx) {
    return l - j(l,nx)*(nx+1);
}

double ro1(double x, double y, double x_max, double y_max, double sigma) {
    return exp(-pow(x-0.25*x_max,2)/sigma/sigma-pow(y-0.5*y_max,2)/sigma/sigma);
}

double ro2(double x, double y, double x_max, double y_max, double sigma) {
    return -exp(-pow(x-0.75*x_max,2)/sigma/sigma-pow(y-0.5*y_max,2)/sigma/sigma);
}

double ro(int l, int nx, double x_max, double y_max, double sigma, double delta) {
    int tmp_i = i(l, nx);
    int tmp_j = j(l, nx);
    return ro1(tmp_i*delta,tmp_j*delta, x_max, y_max, sigma)+ro2(tmp_i*delta,tmp_j*delta, x_max, y_max, sigma);
}

double epsilon(int l, int nx, int epsilon1, int epsilon2) {
    if(i(l,nx) <= (double)nx/2)
        return epsilon1;
    else
        return epsilon2;
}



void pot_sparse(int nx, int ny, int epsilon1, int epsilon2, int is_test, int is_zero, double V1, double V2, double V3, double V4, char* file_name) {
    double delta = 0.1;
    int N = (nx+1)*(ny+1);
    double x_max = nx*delta;
    double y_max = ny*delta;
    double sigma = x_max/10;
    int k = -1;
    int nz_num;
    int edge;
    int V_edge;

    double a[5*N];
    int ja[5*N];
    int ia[N+1];
    double b[N];
    double V[N];

    for(int s = 0; s < 5*N; s++) {
        a[s] = 0;
        ja[s] = 0;
    }
    for(int s = 0; s < N+1; s++)
        ia[s] = -1;
    for(int s = 0; s < N; s++) {
        b[s] = 0;
        V[s] = 0;
    }

    for(int l = 0; l < N; l++) {
        edge = 0;
        V_edge = 0;

        if(i(l,nx) == 0) {
            edge = 1;
            V_edge = V1;
        }

        if(j(l,nx) == ny) {
            edge = 1;
            V_edge = V2;
        }

        if(i(l,nx) == nx) {
            edge = 1;
            V_edge = V3;
        }

        if(j(l,nx) == 0) {
            edge = 1;
            V_edge = V4;
        }

        if(is_zero == 1)
            b[l] = 0;
        else
            b[l]=-ro(l, nx, x_max, y_max, sigma, delta);

        if(edge == 1)
            b[l] = V_edge;


        ia[l] = -1; 

        if(l-nx-1>=0 && edge == 0) {
            k++;
            if(ia[l] < 0)
                ia[l] = k;
            a[k] = epsilon(l,nx,epsilon1,epsilon2)/delta/delta; 
            ja[k] = l - nx - 1;
        }

        if(l-1>=0 && edge == 0) {
            k++;
            if(ia[l] < 0)
                ia[l] = k;
            a[k] = epsilon(l,nx,epsilon1,epsilon2)/delta/delta; 
            ja[k] = l - 1;
        }

        k++;
        if(ia[l] < 0)
            ia[l] = k;
        if(edge == 0) 
            a[k] = -(2*epsilon(l,nx,epsilon1,epsilon2)+epsilon(l+1,nx,epsilon1,epsilon2)+epsilon(l+nx+1,nx,epsilon1,epsilon2))/delta/delta;
        else
            a[k] = 1;
        ja[k] = l;

        if(1 < N && edge == 0) {
            k++;
            a[k] = epsilon(l+1,nx,epsilon1,epsilon2)/delta/delta;
            ja[k] = l + 1;
        }

        if(1< N-nx-1 && edge == 0) {
            k++;
            a[k] = epsilon(l+nx+1,nx,epsilon1,epsilon2)/delta/delta;
            ja[k] = l + nx + 1;
        }
    }

    nz_num = k+1;
    ia[N] = nz_num;


    pmgmres_ilu_cr(N, nz_num, ia, ja, a, V, b, 500, 500, 1e-8, 1e-8);

    if(is_test == 1) {
        FILE *fptrB;
        fptrB = fopen("B.txt","w");
        for(int l = 0; l < N; l++) {
            fprintf(fptrB,"%d\t%f\t%f\t%f\n",l,i(l,nx),j(l,nx),b[l]);
        }
        fclose(fptrB);

        FILE *fptrA;
        fptrA = fopen("A.txt","w");
        for(int q = 0; q < 5*N; q++) 
        fprintf(fptrA,"%d\t%f\n",q,a[q]);
        fclose(fptrA);
    } else {
        FILE *fptr;
        fptr = fopen(file_name,"w");
        for(int l = 0; l < N; l++) {
            fprintf(fptr,"%f\t%f\t%f\n",i(l,nx),j(l,nx),V[l]);
        }
        fclose(fptr);
    }
}


int main() {
    pot_sparse(4, 4, 1, 1, 1, 1, 10, -10, 10, -10, "");
    pot_sparse(50, 50, 1, 1, 2, 1, 10, -10, 10, -10, "5a.txt");
    pot_sparse(100, 100, 1, 1, 2, 1, 10, -10, 10, -10, "5b.txt");
    pot_sparse(200, 200, 1, 1, 2, 1, 10, -10, 10, -10, "5c.txt");
    pot_sparse(100, 100, 1, 1, 2, 0, 0, 0, 0, 0, "6a.txt");
    pot_sparse(100, 100, 1, 2, 2, 0, 0, 0, 0, 0, "6b.txt");
    pot_sparse(100, 100, 1, 10, 2, 0, 0, 0, 0, 0, "6c.txt");

    return 0;
}

