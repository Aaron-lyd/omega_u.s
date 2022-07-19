/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * ntp_slope.c
 *
 * Code generation for function 'ntp_slope'
 *
 */

/* Include files */
#include "ntp_slope.h"
#include "ntp_midpoint_to_casts.h"
#include "ntp_slope_data.h"
#include "ntp_slope_emxutil.h"
#include "ntp_slope_types.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"

/* Function Definitions */
void ntp_slope(const emxArray_real_T *Sppc, const emxArray_real_T *Tppc,
               const emxArray_real_T *Z, const emxArray_real_T *z, real_T tolz,
               emxArray_real_T *dx, emxArray_real_T *dy, emxArray_real_T *sx,
               emxArray_real_T *sy)
{
  emxArray_boolean_T *b;
  emxArray_boolean_T *r;
  emxArray_boolean_T *x;
  emxArray_int32_T *y;
  emxArray_real_T *K;
  real_T Sppc_data[152];
  real_T Tppc_data[152];
  real_T b_Sppc_data[152];
  real_T b_Tppc_data[152];
  real_T Zm_data[20];
  real_T Zn_data[20];
  real_T km;
  real_T zm;
  int32_T Sppc_size[2];
  int32_T b_Sppc_size[2];
  int32_T Zm_size;
  int32_T Zn_size;
  int32_T b_i;
  int32_T i;
  int32_T i1;
  int32_T i2;
  int32_T j;
  int32_T k;
  int32_T loop_ub;
  int32_T ni;
  int32_T npages;
  int32_T vlen;
  int32_T xpageoffset;
  int16_T sqsz[2];
  boolean_T Zmat;
  boolean_T exitg1;
  emlrtHeapReferenceStackEnterFcnR2012b(emlrtRootTLSGlobal);
  /* NTP_SLOPE  Find the slope of the Neutral Tangent Plane at all points along
   * a surface */
  /*  */
  /*  */
  /*  [sx,sy] = ntp_slope(Sppc, Tppc, Z, z, tolz, dx, dy) */
  /*  finds the slopes in the two horizontal directions, sx and sy, of the */
  /*  neutral tangent plane at all points on a surface.  Essentially, this runs
   */
  /*  ntp_midpoint_to_casts at each point on a surface in the two directions, */
  /*  so see ntp_midpoint_to_casts for further documentation. */
  /*  */
  /*  */
  /*  --- Input: */
  /*  Sppc [O, K-1, ni, nj]: coefficients for piecewise polynomial for */
  /*                         practical / Absolute Salinity in terms of Z */
  /*  Tppc [O, K-1, ni, nj]: coefficients for piecewise polynomial for */
  /*                         potential / Conservative Temperature in terms of Z
   */
  /*  Z [K, ni, nj] or [K, 1]: pressure or depth in water column */
  /*  z [ni, nj]: pressure or depth of a surface */
  /*  tolz [1, 1]: tolerance for solving the level of neutral buoyancy (same */
  /*               units as Z) */
  /*  dx [ni, nj]:  distance between grid points (i,j) and (i-1,j) */
  /*  dy [ni, nj]:  distance between  grid points (i,j) and (i,j-1) */
  /*  */
  /*  Note: physical units for Sppc, Tppc, Z, z are determined by eos.m.  See */
  /*        "Non-Boussinesq case" below. */
  /*  */
  /*  Note: Z must increase monotonically along its first dimension. */
  /*  */
  /*  Note: dx and dy can have any dimension replaced by a singleton dimension.
   */
  /*  */
  /*  */
  /*  --- Output: */
  /*  sx [ni, nj]: slope of the neutral tangent plane between (i,j) and (i-1,j)
   */
  /*  sy [ni, nj]: slope of the neutral tangent plane between (i,j) and (i,j-1)
   */
  /*  */
  /*  */
  /*  --- Non-Boussinesq case: */
  /*  in the non-Boussinesq case where Z and z are actually pressure then the */
  /*  units of sx and sy are [dbar / m].  To convert to [dbar / dbar] = [1], */
  /*  divide sx and sy by dp/dz which can be obtained from hydrostatic balance.
   */
  /*  Specifically, do the following. */
  /*                Pa2db = 1e-4; */
  /*                grav = 9.81; % [m s-2]  adjust as needed */
  /*                [s,t] = ppc_val2(Z, Sppc, Tppc, z); */
  /*                if eos(34.5, 3, 1000) > 1  % eos gives density */
  /*                  dpdz = Pa2db * grav * eos(s,t,z); */
  /*                else  % eos gives specific volume */
  /*                  dpdz = Pa2db * grav ./ eos(s,t,z); */
  /*                end */
  /*                sx = sx * 2 ./ (dpdz + circshift(dpdz, [+1, 0])); */
  /*                sy = sy * 2 ./ (dpdz + circshift(dpdz, [0, +1])); */
  /*  A better approach still would be to know pressure as a function of depth.
   */
  /*  This is not provided for here. */
  /*  */
  /*  */
  /*  --- See Also: */
  /*  ntp_midpoint_to_casts */
  /*  ntp_slope_error */
  /*  ppc_linterp, ppc_pchip */
  /*  Author(s) : Geoff Stanley */
  /*  Email     : g.stanley@unsw.edu.au */
  /*  Email     : geoffstanley@gmail.com */
  /*  --- Input checks, set parameters */
  if ((dx->size[0] == 0) || (dx->size[1] == 0)) {
    i = dx->size[0] * dx->size[1];
    dx->size[0] = 1;
    dx->size[1] = 1;
    emxEnsureCapacity_real_T(dx, i);
    dx->data[0] = 1.0;
  }
  if ((dy->size[0] == 0) || (dy->size[1] == 0)) {
    i = dy->size[0] * dy->size[1];
    dy->size[0] = 1;
    dy->size[1] = 1;
    emxEnsureCapacity_real_T(dy, i);
    dy->data[0] = 1.0;
  }
  ni = z->size[0];
  /*  --- Solve nonlinear problem for NTP depth difference between each pair of
   * adjacent casts */
  Zmat = (((Z->size[0] != 1) && (Z->size[1] != 1)) || (Z->size[2] != 1));
  k = Z->size[0];
  Zm_size = Z->size[0];
  for (i = 0; i < k; i++) {
    Zm_data[i] = Z->data[i];
  }
  k = Z->size[0];
  Zn_size = Z->size[0];
  for (i = 0; i < k; i++) {
    Zn_data[i] = Z->data[i];
  }
  emxInit_boolean_T(&b, 4, true);
  /*  K gives the number of valid bottles in each water column. */
  /*  Note that K >= 1.  Even where all bottles are invalid, i.e. land, K = 1.
   */
  k = Sppc->size[1];
  loop_ub = Sppc->size[2];
  xpageoffset = Sppc->size[3];
  i = b->size[0] * b->size[1] * b->size[2] * b->size[3];
  b->size[0] = 1;
  b->size[1] = Sppc->size[1];
  b->size[2] = Sppc->size[2];
  b->size[3] = Sppc->size[3];
  emxEnsureCapacity_boolean_T(b, i);
  for (i = 0; i < xpageoffset; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      for (i2 = 0; i2 < k; i2++) {
        b->data[(i2 + b->size[1] * i1) + b->size[1] * b->size[2] * i] =
            muDoubleScalarIsInf(
                Sppc->data[(Sppc->size[0] * i2 +
                            Sppc->size[0] * Sppc->size[1] * i1) +
                           Sppc->size[0] * Sppc->size[1] * Sppc->size[2] * i]);
      }
    }
  }
  emxInit_boolean_T(&r, 4, true);
  k = Sppc->size[1];
  loop_ub = Sppc->size[2];
  xpageoffset = Sppc->size[3];
  i = r->size[0] * r->size[1] * r->size[2] * r->size[3];
  r->size[0] = 1;
  r->size[1] = Sppc->size[1];
  r->size[2] = Sppc->size[2];
  r->size[3] = Sppc->size[3];
  emxEnsureCapacity_boolean_T(r, i);
  for (i = 0; i < xpageoffset; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      for (i2 = 0; i2 < k; i2++) {
        r->data[(i2 + r->size[1] * i1) + r->size[1] * r->size[2] * i] =
            muDoubleScalarIsNaN(
                Sppc->data[(Sppc->size[0] * i2 +
                            Sppc->size[0] * Sppc->size[1] * i1) +
                           Sppc->size[0] * Sppc->size[1] * Sppc->size[2] * i]);
      }
    }
  }
  k = b->size[1] * b->size[2] * b->size[3];
  i = b->size[0] * b->size[1] * b->size[2] * b->size[3];
  b->size[0] = 1;
  emxEnsureCapacity_boolean_T(b, i);
  for (i = 0; i < k; i++) {
    b->data[i] = ((!b->data[i]) && (!r->data[i]));
  }
  emxFree_boolean_T(&r);
  emxInit_int32_T(&y, 4, true);
  vlen = b->size[1];
  if ((b->size[1] == 0) || (b->size[2] == 0) || (b->size[3] == 0)) {
    i = y->size[0] * y->size[1] * y->size[2] * y->size[3];
    y->size[0] = 1;
    y->size[1] = 1;
    y->size[2] = (int16_T)b->size[2];
    y->size[3] = (int16_T)b->size[3];
    emxEnsureCapacity_int32_T(y, i);
    k = (int16_T)b->size[2] * (int16_T)b->size[3];
    for (i = 0; i < k; i++) {
      y->data[i] = 0;
    }
  } else {
    npages = 1;
    k = 4;
    while ((k > 2) && (b->size[k - 1] == 1)) {
      k--;
    }
    for (xpageoffset = 3; xpageoffset <= k; xpageoffset++) {
      npages *= b->size[xpageoffset - 1];
    }
    i = y->size[0] * y->size[1] * y->size[2] * y->size[3];
    y->size[0] = 1;
    y->size[1] = 1;
    y->size[2] = (int16_T)b->size[2];
    y->size[3] = (int16_T)b->size[3];
    emxEnsureCapacity_int32_T(y, i);
    for (b_i = 0; b_i < npages; b_i++) {
      xpageoffset = b_i * b->size[1];
      y->data[b_i] = b->data[xpageoffset];
      for (k = 2; k <= vlen; k++) {
        y->data[b_i] += b->data[(xpageoffset + k) - 1];
      }
    }
  }
  emxFree_boolean_T(&b);
  sqsz[0] = 1;
  sqsz[1] = 1;
  k = 4;
  while ((k > 2) && (y->size[k - 1] == 1)) {
    k--;
  }
  if (k == 2) {
    sqsz[1] = 1;
  } else {
    j = 0;
    if (y->size[2] != 1) {
      j = 1;
      sqsz[0] = (int16_T)y->size[2];
    }
    if (y->size[3] != 1) {
      sqsz[j] = (int16_T)y->size[3];
    }
  }
  emxInit_real_T(&K, 2, true);
  i = K->size[0] * K->size[1];
  K->size[0] = sqsz[0];
  K->size[1] = sqsz[1];
  emxEnsureCapacity_real_T(K, i);
  xpageoffset = sqsz[0] * sqsz[1];
  for (i = 0; i < xpageoffset; i++) {
    K->data[i] = (real_T)y->data[i] + 1.0;
  }
  emxFree_int32_T(&y);
  /*  Loop over each water column.   */
  /*  Note this calculation assumes a doubly-periodic domain */
  /*  The NTP slope is grad_n z, where z < 0 */
  /*  Negative sign added to ntp_midpoint_to_casts results because we've been z
   * > 0 people. */
  i = z->size[1];
  i1 = sx->size[0] * sx->size[1];
  sx->size[0] = (int16_T)z->size[0];
  sx->size[1] = (int16_T)z->size[1];
  emxEnsureCapacity_real_T(sx, i1);
  i1 = sy->size[0] * sy->size[1];
  sy->size[0] = (int16_T)z->size[0];
  sy->size[1] = (int16_T)z->size[1];
  emxEnsureCapacity_real_T(sy, i1);
  for (j = 0; j < i; j++) {
    xpageoffset = j - 1;
    if (z->size[1] == 0) {
      if (j - 1 == 0) {
        xpageoffset = 0;
      }
    } else if (j - 1 == 0) {
      xpageoffset = 0;
    } else {
      xpageoffset =
          (int32_T)muDoubleScalarRem(((real_T)j + 1.0) - 2.0, z->size[1]);
      if ((xpageoffset != 0) && (j - 1 < 0)) {
        xpageoffset += z->size[1];
      }
    }
    for (b_i = 0; b_i < ni; b_i++) {
      vlen = b_i - 1;
      if (ni == 0) {
        if (b_i - 1 == 0) {
          vlen = 0;
        }
      } else if (b_i - 1 == 0) {
        vlen = 0;
      } else {
        vlen = (int32_T)muDoubleScalarRem(((real_T)b_i + 1.0) - 2.0, ni);
        if ((vlen != 0) && (b_i - 1 < 0)) {
          vlen += ni;
        }
      }
      zm = z->data[b_i + z->size[0] * j];
      km = K->data[b_i + K->size[0] * j];
      if (Zmat) {
        k = Z->size[0];
        Zm_size = Z->size[0];
        for (i1 = 0; i1 < k; i1++) {
          Zm_data[i1] =
              Z->data[(i1 + Z->size[0] * b_i) + Z->size[0] * Z->size[1] * j];
        }
        k = Z->size[0];
        Zn_size = Z->size[0];
        for (i1 = 0; i1 < k; i1++) {
          Zn_data[i1] =
              Z->data[(i1 + Z->size[0] * vlen) + Z->size[0] * Z->size[1] * j];
        }
      }
      /*  --- NTP with neighbour in i dimension */
      k = Sppc->size[0];
      loop_ub = Sppc->size[1];
      Sppc_size[0] = Sppc->size[0];
      Sppc_size[1] = Sppc->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        for (i2 = 0; i2 < k; i2++) {
          Sppc_data[i2 + Sppc_size[0] * i1] =
              Sppc->data[((i2 + Sppc->size[0] * i1) +
                          Sppc->size[0] * Sppc->size[1] * b_i) +
                         Sppc->size[0] * Sppc->size[1] * Sppc->size[2] * j];
        }
      }
      k = Tppc->size[0];
      loop_ub = Tppc->size[1];
      npages = Tppc->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        for (i2 = 0; i2 < k; i2++) {
          Tppc_data[i2 + npages * i1] =
              Tppc->data[((i2 + Tppc->size[0] * i1) +
                          Tppc->size[0] * Tppc->size[1] * b_i) +
                         Tppc->size[0] * Tppc->size[1] * Tppc->size[2] * j];
        }
      }
      k = Sppc->size[0];
      loop_ub = Sppc->size[1];
      b_Sppc_size[0] = Sppc->size[0];
      b_Sppc_size[1] = Sppc->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        for (i2 = 0; i2 < k; i2++) {
          b_Sppc_data[i2 + b_Sppc_size[0] * i1] =
              Sppc->data[((i2 + Sppc->size[0] * i1) +
                          Sppc->size[0] * Sppc->size[1] * vlen) +
                         Sppc->size[0] * Sppc->size[1] * Sppc->size[2] * j];
        }
      }
      k = Tppc->size[0];
      loop_ub = Tppc->size[1];
      npages = Tppc->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        for (i2 = 0; i2 < k; i2++) {
          b_Tppc_data[i2 + npages * i1] =
              Tppc->data[((i2 + Tppc->size[0] * i1) +
                          Tppc->size[0] * Tppc->size[1] * vlen) +
                         Tppc->size[0] * Tppc->size[1] * Tppc->size[2] * j];
        }
      }
      sx->data[b_i + sx->size[0] * j] = -ntp_midpoint_to_casts(
          Sppc_data, Sppc_size, Tppc_data, Zm_data, Zm_size, km, b_Sppc_data,
          b_Sppc_size, b_Tppc_data, Zn_data, Zn_size,
          K->data[vlen + K->size[0] * j], zm, z->data[vlen + z->size[0] * j],
          tolz);
      /*  --- NTP with neighbour in j dimension */
      if (Zmat) {
        k = Z->size[0];
        Zn_size = Z->size[0];
        for (i1 = 0; i1 < k; i1++) {
          Zn_data[i1] = Z->data[(i1 + Z->size[0] * b_i) +
                                Z->size[0] * Z->size[1] * xpageoffset];
        }
      }
      k = Sppc->size[0];
      loop_ub = Sppc->size[1];
      Sppc_size[0] = Sppc->size[0];
      Sppc_size[1] = Sppc->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        for (i2 = 0; i2 < k; i2++) {
          Sppc_data[i2 + Sppc_size[0] * i1] =
              Sppc->data[((i2 + Sppc->size[0] * i1) +
                          Sppc->size[0] * Sppc->size[1] * b_i) +
                         Sppc->size[0] * Sppc->size[1] * Sppc->size[2] * j];
        }
      }
      k = Tppc->size[0];
      loop_ub = Tppc->size[1];
      npages = Tppc->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        for (i2 = 0; i2 < k; i2++) {
          Tppc_data[i2 + npages * i1] =
              Tppc->data[((i2 + Tppc->size[0] * i1) +
                          Tppc->size[0] * Tppc->size[1] * b_i) +
                         Tppc->size[0] * Tppc->size[1] * Tppc->size[2] * j];
        }
      }
      k = Sppc->size[0];
      loop_ub = Sppc->size[1];
      b_Sppc_size[0] = Sppc->size[0];
      b_Sppc_size[1] = Sppc->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        for (i2 = 0; i2 < k; i2++) {
          b_Sppc_data[i2 + b_Sppc_size[0] * i1] =
              Sppc->data[((i2 + Sppc->size[0] * i1) +
                          Sppc->size[0] * Sppc->size[1] * b_i) +
                         Sppc->size[0] * Sppc->size[1] * Sppc->size[2] *
                             xpageoffset];
        }
      }
      k = Tppc->size[0];
      loop_ub = Tppc->size[1];
      npages = Tppc->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        for (i2 = 0; i2 < k; i2++) {
          b_Tppc_data[i2 + npages * i1] =
              Tppc->data[((i2 + Tppc->size[0] * i1) +
                          Tppc->size[0] * Tppc->size[1] * b_i) +
                         Tppc->size[0] * Tppc->size[1] * Tppc->size[2] *
                             xpageoffset];
        }
      }
      sy->data[b_i + sy->size[0] * j] = -ntp_midpoint_to_casts(
          Sppc_data, Sppc_size, Tppc_data, Zm_data, Zm_size, km, b_Sppc_data,
          b_Sppc_size, b_Tppc_data, Zn_data, Zn_size,
          K->data[b_i + K->size[0] * xpageoffset], zm,
          z->data[b_i + z->size[0] * xpageoffset], tolz);
    }
    /*  i */
  }
  emxFree_real_T(&K);
  emxInit_boolean_T(&x, 1, true);
  /*  j */
  /*  Divide by horizontal distances */
  i = x->size[0];
  x->size[0] = dx->size[0] * dx->size[1];
  emxEnsureCapacity_boolean_T(x, i);
  k = dx->size[0] * dx->size[1];
  for (i = 0; i < k; i++) {
    x->data[i] = (dx->data[i] == 1.0);
  }
  Zmat = true;
  xpageoffset = 1;
  exitg1 = false;
  while ((!exitg1) && (xpageoffset <= x->size[0])) {
    if (!x->data[xpageoffset - 1]) {
      Zmat = false;
      exitg1 = true;
    } else {
      xpageoffset++;
    }
  }
  if (!Zmat) {
    i = x->size[0];
    x->size[0] = dy->size[0] * dy->size[1];
    emxEnsureCapacity_boolean_T(x, i);
    k = dy->size[0] * dy->size[1];
    for (i = 0; i < k; i++) {
      x->data[i] = (dy->data[i] == 1.0);
    }
    Zmat = true;
    xpageoffset = 1;
    exitg1 = false;
    while ((!exitg1) && (xpageoffset <= x->size[0])) {
      if (!x->data[xpageoffset - 1]) {
        Zmat = false;
        exitg1 = true;
      } else {
        xpageoffset++;
      }
    }
    if (!Zmat) {
      k = sx->size[0] * sx->size[1];
      for (i = 0; i < k; i++) {
        sx->data[i] /= dx->data[i];
      }
      k = sy->size[0] * sy->size[1];
      for (i = 0; i < k; i++) {
        sy->data[i] /= dy->data[i];
      }
    }
  }
  emxFree_boolean_T(&x);
  emlrtHeapReferenceStackLeaveFcnR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (ntp_slope.c) */
