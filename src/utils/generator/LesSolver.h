#ifndef LES_SOLVER_H
#define LES_SOLVER_H

#ifdef __cplusplus
extern "C" {
#endif

/** Linear equation system solver
 *
 * Solves simple n dimensional linear systems Ax = a
 * The determinant of A has to be unequal zero.
 *
 * @param[in]  A  NxN matrix
 * @param[in]  a  vector
 * @param[in]  n  dimensionality of the system
 * @param[out] x  result vector
 */
void LesSolve(double **A, double *a, int n, double *x);

#ifdef __cplusplus
} /*extern "C" */
#endif

#endif /*LES_SOLVER_H */
