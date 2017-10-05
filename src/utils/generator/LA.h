#ifndef LA_H
#define LA_H

/** Print vector
 * @param[in]  a  pointer to array
 * @param[in]  n  length of array
 */
void print_vector(double *a, int n);

/** Print square matrix
 * @param[in]  A  pointer to array
 * @param[in]  n  length of array
 */
void print_matrix(double **A, int n);

/** Calculate absolute value of vector
 * @param[in]  v  pointer to array
 * @param[in]  n  length of array
 * @return  absolute avlue
 */
double vec_abs(double v[3], int n);

/** Calculate scalar product of two vectors
 * @param[in]  v  pointer to array
 * @param[in]  n  length of array
 * @return  scalar product
 */
double vec_scalar_product(double v1[3], double v2[3], int n);

#endif /* LA_H */
