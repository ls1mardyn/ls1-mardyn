#include "LU.h"

#include <stdio.h>

void LU(double **A, int n, double **L, double **R) {
    int i, j, k;
    for(k = 0; k < n; k++) {
        R[0][k] = A[0][k];
        L[k][k] = 1;
    }
    for(i = 1; i < n; i++) {
        L[i][0] = A[i][0] / A[0][0];
        for(k = 1; k < i; k++) {
            L[i][k] = A[i][k];
            for(j = 0; j < k; j++) {
                L[i][k] -= L[i][j]*R[j][k];
            }
            L[i][k] /= R[k][k];
        }
        for(k = i; k < n; k++) {
            R[i][k] = A[i][k];
            for(j = 0; j < i; j++) {
                R[i][k] -= L[i][j]*R[j][k];
            }
        }
    }
}

