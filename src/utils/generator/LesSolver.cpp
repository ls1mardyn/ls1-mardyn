#include "LesSolver.h"

#include "LA.h"
#include "LU.h"

void LesSolve(double **A, double *a, int n, double *x) {
    int i, k;
    double **L;
    double **R;
    L = new double*[n];
    R = new double*[n];
    for(i = 0; i < n; i++) {
        L[i] = new double[n];
        R[i] = new double[n];
    }
    LU(A, n, L, R);
    //print_matrix(L,n);
    //print_matrix(R,n);
    //print_vector(a,n);

    double *r = new double[n];

    /* forward substitution */
    r[0] = a[0];
    for(i = 1; i < n; i++) {
        r[i] = a[i];
        for(k = 0; k < i; k++) {
            r[i] -= L[i][k] * r[k];
        }
    }

    //print_vector(r,n);

    /* backward substitution */
    x[n-1] = r[n-1] / R[n-1][n-1];
    for(i = n-2; i >= 0; i--) {
        x[i] = r[i];
        for(k = i+1; k < n; k++) {
            x[i] -= R[i][k] * x[k];
        }
        x[i] /= R[i][i];
    }
    for(i = 0; i < n; i++) {
		delete[] L[i];
		delete[] R[i];
	}
    delete[] L;
	delete[] R;
	delete[] r;
}
