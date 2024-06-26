
/**

C++ module 'eig3' by Connelly Barnes
------------------------------------

License: public domain.

The source files in this directory have been copied from the public domain
Java matrix library JAMA.  The derived source code is in the public domain
as well.

*/

/* Eigen-decomposition for symmetric 3x3 real matrices.
   Public domain, copied from the public domain Java library JAMA. */

#ifndef _eig_h

/* Symmetric matrix A => eigenvectors in columns of V, corresponding
   eigenvalues in d. */
void eigen_decomposition(double A[3][3], double V[3][3], double d[3]);

#endif
