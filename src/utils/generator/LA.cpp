#include "LA.h"

#include <math.h>
#include <stdio.h>

void print_vector(double *a, int n) {
    int i;
    printf("[");
    for(i = 0; i < n-1; i++) {
        printf("%lf,", a[i]);
    }
    printf("%lf]\n", a[n-1]);
}

void print_matrix(double **A, int n) {
    int i, j;
    printf("[");
    for(i = 0; i < n; i++) {
        printf("[");
        for(j = 0; j < n-1; j++) {
            printf("%lf,", A[i][j]);
        }
        printf("%lf]", A[i][j]);
        if(i != n-1) {
            printf(",");
        }
    }
    printf("]\n");
}

double vec_abs(double v[3], int n) {
    return sqrt(vec_scalar_product(v, v, n));
}

double vec_scalar_product(double v1[3], double v2[3], int n) {
    int i;
    double result = 0;
    for(i = 0; i < n; i++) {
        result += v1[i]*v2[i];
    }
    return result;
}
