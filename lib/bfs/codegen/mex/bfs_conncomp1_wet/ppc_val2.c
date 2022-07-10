/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * ppc_val2.c
 *
 * Code generation for function 'ppc_val2'
 *
 */

/* Include files */
#include "ppc_val2.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"

/* Function Definitions */
void ppc_val2(const real_T X_data[], int32_T X_size, const real_T C_data[],
              const int32_T C_size[2], const real_T D_data[], real_T x,
              real_T *y, real_T *z)
{
  real_T p;
  real_T t;
  int32_T b_i;
  int32_T i;
  int32_T i1;
  int32_T j;
  int32_T k;
  int32_T o;
  int32_T y_size_idx_1;
  int8_T y_data[15];
  int8_T szC_idx_0;
  /* PPC_VAL2  Piecewise Polynomial Evaluation, quick and twice */
  /*  */
  /*  */
  /*  [y,z] = ppc_val2(X, C, D, x) */
  /*  evaluates the piecewise polynomials whose coefficients are C and D and */
  /*  whose knots are X, at data sites x. */
  /*  */
  /*  [y,z] = ppc_val2(..., d) */
  /*  as above but evaluates the d'th derivative. */
  /*  */
  /*  */
  /*  --- Input: */
  /*  X [    K   x M], knots of the piecewise polynomials */
  /*  C [O x K-1 x N], coefficients of the first piecewise polynomial */
  /*  D [O x K-1 x N], coefficients of the second piecewise polynomial */
  /*  x [    L   x M], evaluation sites */
  /*  d [1, 1]       , order of derivative to evaluate; defaults to 0 if not */
  /*                   given, which means to evaluate with no differentiation.
   */
  /*  */
  /*  */
  /*  --- Output: */
  /*  y [L x N], the first piecewise polynomial evaluated at x */
  /*  z [L x N], the second piecewise polynomial evaluated at x */
  /*  */
  /*  */
  /*  --- Notes: */
  /*  X(:,n) must be monotonically increasing for all n. */
  /*  */
  /*  NaN's in X are treated as +Inf (and as such must come at the end of each
   */
  /*  column). */
  /*  */
  /*  M, the last dimension of X and x, can be either 1 or N. */
  /*  */
  /*  Any dimension N can actually be higher-dimensional, so long as it has N */
  /*  elements. */
  /*  */
  /*  Even if L == 1, x does need a leading singleton dimension. */
  /*  */
  /*  If L == 1, this dimension is squeeze()'d out of y and z. */
  /*  */
  /*  */
  /*  --- Acknowledgements: */
  /*  This code is adapted from MATLAB's ppval.m */
  /*  Author    : Geoff Stanley */
  /*  Email     : geoffstanley@gmail.com */
  /*  Version   : 1.0 */
  /*  History   : 24/10/2019 - initial release */
  /*  Set default to evaluate, not differentiate: */
  szC_idx_0 = (int8_T)C_size[0];
  /*  Order of the piecewise polynomial */
  /*  number of knots of the piecewise polynomials */
  /*  number of levels to interpolate */
  /*  Handle case when L should be 1. */
  /*  Add a trailing 1, to ensure 2 dimensions at least */
  *y = rtNaN;
  *z = rtNaN;
  /*  Evaluate each piecewise polynomial */
  /*  used for linear indexing */
  /*  used for linear indexing */
  /*  used for linear indexing */
  /*  used for linear indexing */
  /*  x(l,n) or x(l,1) as appropriate */
  if ((!muDoubleScalarIsNaN(x)) && (!muDoubleScalarIsNaN(X_data[0])) &&
      (!(x < X_data[0])) && (!(x > X_data[(int8_T)X_size - 1]))) {
    /*  Leftmost binary search to find i such that: */
    /*  i = 1                      if xln <= X(1), or */
    /*  i = M                      if X(M) < xln */
    /*  X(i-1,n) < xln <= X(i,n)   otherwise */
    /*  We use leftmost so that NaN's at the end of X are treated as though they
     * are Inf. */
    i = 1;
    /*  Result will be >= 1 always */
    k = (int8_T)X_size;
    /*  Result will be <= M always */
    while (i < k) {
      j = (int32_T)muDoubleScalarFloor((real_T)(i + k) / 2.0);
      if (X_data[j - 1] < x) {
        /*  [X(j,n) or X(j,1) as appropriate]  < xln */
        i = j + 1;
      } else {
        k = j;
      }
    }
    if (i == 1) {
      /*  Note: X(nX + 1) == xln   is guaranteed */
      /*  Build integer multiplying the coefficient */
      *y = C_data[(int8_T)C_size[0] - 1];
      /*  y(l,n) = C(O - d, i, n) * p; */
      *z = D_data[(int8_T)C_size[0] - 1];
    } else {
      /*  Evaluate this piece of the polynomial (see ppval.m) */
      t = x - X_data[i - 2];
      /*  Switch to local coordinates */
      /*  Overload variable i, to speed with indexing */
      /*  subtract 1 from i so that 1 <= i <= M-1, and X(i) <= xln < X(i+1), */
      /*  subtract another 1 for indexing.  */
      i = (i - 2) * (int8_T)C_size[0];
      *y = 0.0;
      *z = 0.0;
      b_i = (int8_T)C_size[0];
      for (o = 0; o < b_i; o++) {
        j = szC_idx_0 - o;
        if (j - 1 < j) {
          y_size_idx_1 = 0;
        } else {
          k = (int8_T)(j - 1) - (int8_T)j;
          y_size_idx_1 = k + 1;
          for (i1 = 0; i1 <= k; i1++) {
            y_data[i1] = (int8_T)((int8_T)j + (int8_T)i1);
          }
        }
        if (y_size_idx_1 == 0) {
          p = 1.0;
        } else {
          p = y_data[0];
          for (k = 2; k <= y_size_idx_1; k++) {
            p *= (real_T)y_data[k - 1];
          }
        }
        /*  Build integer multiplying the coefficient */
        j = o + i;
        *y = t * *y + C_data[j] * p;
        /*  y(l,n) = t * y(l,n) + C(o,i+1,n+1) * p; */
        *z = t * *z + D_data[j] * p;
      }
    }
  }
  /*  for l */
  /*  for n */
  /*  Reshape output to be like input.  Also remove leading dimensions if L ==
   */
  /*  1, but leave row vectors as row vectors. */
}

/* End of code generation (ppc_val2.c) */
