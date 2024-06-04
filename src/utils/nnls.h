/**
 * Algorithm NNLS: NONNEGATIVE LEAST SQUARES
 *
 * The original version of this code was developed by
 * Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
 * 1973 JUN 15, and published in the book
 * "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
 * Revised FEB 1995 to accompany reprinting of the book by SIAM.
 *
 * Given an m by n matrix, A, and an m-vector, b,  compute an n-vector, x, that solves the least squares problem
 *          A * x = b  subject to x >= 0
 * @param An m by n dimensional matrix (column-major order). On entry a contains the m by n matrix, A. On exit, A
 * contains the product matrix, Q*A , where Q is an m by m orthogonal matrix generated implicitly by this subroutine.
 *
 * @param mda The first dimensioning parameter for the ARRAY A. Should be the same as m.
 * @param m The first dimension of A.
 * @param n The second dimension of A.
 * @param b On entry b contains the m-vector, b. On exit b contains Q*b.
 * @param x On entry x need not be initialized. On exit x will contain the solution vector.
 * @param rnorm On exit rnorm contains the euclidean norm of the residual vector.
 * @param w     An n-array of working space. On exit w() will contain the dual solution vector.
 * w will satisfy w(i) = 0. for all i in set p.
 * w(i) <= 0. for all i in set z.
 * @param zz     an m-array of working space.
 * @param index      an integer working array of length at least n.
 *                 on exit the contents of this array define the sets
 *                 p and z as follows..
 *                 index(1)   thru index(nsetp) = set p.
 *                 index(iz1) thru index(iz2)   = set z.
 *                 iz1 = nsetp + 1 = npp1
 *                 iz2 = n
 * @param mode This is a success-failure flag with the following
 *             meanings.
 *             1     the solution has been computed successfully.
 *             2     the dimensions of the problem are bad.
 *                   either m .le. 0 or n .le. 0.
 *             3    iteration count exceeded.  more than 3*n iterations.
 */
int nnls_c(double *A, const int *mda, const int *m, const int *n, double *b, double *x, double *rnorm, double *w,
		   double *zz, int *index, int *mode);