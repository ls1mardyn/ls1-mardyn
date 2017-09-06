#ifndef LU_H
#define LU_H


/** Standard LR decomposition
 *
 * Calculates the LR factorization of a square matrix LR = A.
 * For algorithm see e.g. Numerik-Algorithmen (ISBN:3-540-62669-7).
 *
 * @param[in]  A  pointer to NxN matrix
 * @param[in]  n  dimensionality N of the matrix
 * @param[out] L  lower triangular matrix with L_{ii} = 1
 * @param[out] R  upper triangular matrix
 */
void LU(double **A, int n, double **L, double **U);

#endif /* LU_H */
