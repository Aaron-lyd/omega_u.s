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
#include <string.h>

/* Function Definitions */
void ntp_slope(const emxArray_real_T *Sppc, const emxArray_real_T *Tppc,
               const emxArray_real_T *Z, const real_T z_data[],
               const int32_T z_size[2], real_T tolz, real_T dx_data[],
               int32_T dx_size[2], real_T dy_data[], int32_T dy_size[2],
               real_T sx_data[], int32_T sx_size[2], real_T sy_data[],
               int32_T sy_size[2])
{
  emxArray_boolean_T *b;
  emxArray_boolean_T *r;
  real_T K_data[4225];
  real_T Sppc_data[152];
  real_T Tppc_data[152];
  real_T b_Sppc_data[152];
  real_T b_Tppc_data[152];
  real_T Zm_data[20];
  real_T Zn_data[20];
  real_T km;
  real_T zm;
  int32_T y_data[4225];
  int32_T Sppc_size[2];
  int32_T b_Sppc_size[2];
  int32_T K_size_idx_0;
  int32_T Tppc_size_idx_0;
  int32_T Zm_size;
  int32_T Zn_size;
  int32_T b_i;
  int32_T b_loop_ub;
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
  int8_T y_size[4];
  int8_T sqsz[2];
  boolean_T x_data[4225];
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
  if ((dx_size[0] == 0) || (dx_size[1] == 0)) {
    dx_size[0] = 1;
    dx_size[1] = 1;
    dx_data[0] = 1.0;
  }
  if ((dy_size[0] == 0) || (dy_size[1] == 0)) {
    dy_size[0] = 1;
    dy_size[1] = 1;
    dy_data[0] = 1.0;
  }
  ni = z_size[0];
  /*  --- Solve nonlinear problem for NTP depth difference between each pair of
   * adjacent casts */
  Zmat = (((Z->size[0] != 1) && (Z->size[1] != 1)) || (Z->size[2] != 1));
  loop_ub = Z->size[0];
  Zm_size = Z->size[0];
  for (i = 0; i < loop_ub; i++) {
    Zm_data[i] = Z->data[i];
  }
  loop_ub = Z->size[0];
  Zn_size = Z->size[0];
  for (i = 0; i < loop_ub; i++) {
    Zn_data[i] = Z->data[i];
  }
  emxInit_boolean_T(&b, 4, true);
  /*  K gives the number of valid bottles in each water column. */
  /*  Note that K >= 1.  Even where all bottles are invalid, i.e. land, K = 1.
   */
  loop_ub = Sppc->size[1];
  b_loop_ub = Sppc->size[2];
  vlen = Sppc->size[3];
  i = b->size[0] * b->size[1] * b->size[2] * b->size[3];
  b->size[0] = 1;
  b->size[1] = Sppc->size[1];
  b->size[2] = Sppc->size[2];
  b->size[3] = Sppc->size[3];
  emxEnsureCapacity_boolean_T(b, i);
  for (i = 0; i < vlen; i++) {
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      for (i2 = 0; i2 < loop_ub; i2++) {
        b->data[(i2 + b->size[1] * i1) + b->size[1] * b->size[2] * i] =
            muDoubleScalarIsInf(
                Sppc->data[(Sppc->size[0] * i2 +
                            Sppc->size[0] * Sppc->size[1] * i1) +
                           Sppc->size[0] * Sppc->size[1] * Sppc->size[2] * i]);
      }
    }
  }
  emxInit_boolean_T(&r, 4, true);
  loop_ub = Sppc->size[1];
  b_loop_ub = Sppc->size[2];
  vlen = Sppc->size[3];
  i = r->size[0] * r->size[1] * r->size[2] * r->size[3];
  r->size[0] = 1;
  r->size[1] = Sppc->size[1];
  r->size[2] = Sppc->size[2];
  r->size[3] = Sppc->size[3];
  emxEnsureCapacity_boolean_T(r, i);
  for (i = 0; i < vlen; i++) {
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      for (i2 = 0; i2 < loop_ub; i2++) {
        r->data[(i2 + r->size[1] * i1) + r->size[1] * r->size[2] * i] =
            muDoubleScalarIsNaN(
                Sppc->data[(Sppc->size[0] * i2 +
                            Sppc->size[0] * Sppc->size[1] * i1) +
                           Sppc->size[0] * Sppc->size[1] * Sppc->size[2] * i]);
      }
    }
  }
  loop_ub = b->size[1] * b->size[2] * b->size[3];
  i = b->size[0] * b->size[1] * b->size[2] * b->size[3];
  b->size[0] = 1;
  emxEnsureCapacity_boolean_T(b, i);
  for (i = 0; i < loop_ub; i++) {
    b->data[i] = ((!b->data[i]) && (!r->data[i]));
  }
  emxFree_boolean_T(&r);
  vlen = b->size[1];
  if ((b->size[1] == 0) || (b->size[2] == 0) || (b->size[3] == 0)) {
    y_size[0] = 1;
    y_size[1] = 1;
    y_size[2] = (int8_T)b->size[2];
    y_size[3] = (int8_T)b->size[3];
    loop_ub = (int8_T)b->size[2] * (int8_T)b->size[3];
    if (0 <= loop_ub - 1) {
      memset(&y_data[0], 0, loop_ub * sizeof(int32_T));
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
    y_size[0] = 1;
    y_size[1] = 1;
    y_size[2] = (int8_T)b->size[2];
    y_size[3] = (int8_T)b->size[3];
    for (b_i = 0; b_i < npages; b_i++) {
      xpageoffset = b_i * b->size[1];
      y_data[b_i] = b->data[xpageoffset];
      for (k = 2; k <= vlen; k++) {
        y_data[b_i] += b->data[(xpageoffset + k) - 1];
      }
    }
  }
  emxFree_boolean_T(&b);
  sqsz[0] = 1;
  sqsz[1] = 1;
  k = 4;
  while ((k > 2) && (y_size[k - 1] == 1)) {
    k--;
  }
  if (k == 2) {
    sqsz[1] = 1;
  } else {
    j = 0;
    if (y_size[2] != 1) {
      j = 1;
      sqsz[0] = y_size[2];
    }
    if (y_size[3] != 1) {
      sqsz[j] = y_size[3];
    }
  }
  K_size_idx_0 = sqsz[0];
  vlen = sqsz[0] * sqsz[1];
  for (i = 0; i < vlen; i++) {
    K_data[i] = (real_T)y_data[i] + 1.0;
  }
  /*  Loop over each water column.   */
  /*  Note this calculation assumes a doubly-periodic domain */
  /*  The NTP slope is grad_n z, where z < 0 */
  /*  Negative sign added to ntp_midpoint_to_casts results because we've been z
   * > 0 people. */
  i = z_size[1];
  sx_size[0] = (int8_T)z_size[0];
  sx_size[1] = (int8_T)z_size[1];
  sy_size[0] = (int8_T)z_size[0];
  sy_size[1] = (int8_T)z_size[1];
  for (j = 0; j < i; j++) {
    xpageoffset = j - 1;
    if (z_size[1] == 0) {
      if (j - 1 == 0) {
        xpageoffset = 0;
      }
    } else if (j - 1 == 0) {
      xpageoffset = 0;
    } else {
      xpageoffset =
          (int32_T)muDoubleScalarRem(((real_T)j + 1.0) - 2.0, z_size[1]);
      if ((xpageoffset != 0) && (j - 1 < 0)) {
        xpageoffset += z_size[1];
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
      npages = z_size[0] * j;
      zm = z_data[b_i + npages];
      k = K_size_idx_0 * j;
      km = K_data[b_i + k];
      if (Zmat) {
        loop_ub = Z->size[0];
        Zm_size = Z->size[0];
        for (i1 = 0; i1 < loop_ub; i1++) {
          Zm_data[i1] =
              Z->data[(i1 + Z->size[0] * b_i) + Z->size[0] * Z->size[1] * j];
        }
        loop_ub = Z->size[0];
        Zn_size = Z->size[0];
        for (i1 = 0; i1 < loop_ub; i1++) {
          Zn_data[i1] =
              Z->data[(i1 + Z->size[0] * vlen) + Z->size[0] * Z->size[1] * j];
        }
      }
      /*  --- NTP with neighbour in i dimension */
      loop_ub = Sppc->size[0];
      b_loop_ub = Sppc->size[1];
      Sppc_size[0] = Sppc->size[0];
      Sppc_size[1] = Sppc->size[1];
      for (i1 = 0; i1 < b_loop_ub; i1++) {
        for (i2 = 0; i2 < loop_ub; i2++) {
          Sppc_data[i2 + Sppc_size[0] * i1] =
              Sppc->data[((i2 + Sppc->size[0] * i1) +
                          Sppc->size[0] * Sppc->size[1] * b_i) +
                         Sppc->size[0] * Sppc->size[1] * Sppc->size[2] * j];
        }
      }
      loop_ub = Tppc->size[0];
      b_loop_ub = Tppc->size[1];
      Tppc_size_idx_0 = Tppc->size[0];
      for (i1 = 0; i1 < b_loop_ub; i1++) {
        for (i2 = 0; i2 < loop_ub; i2++) {
          Tppc_data[i2 + Tppc_size_idx_0 * i1] =
              Tppc->data[((i2 + Tppc->size[0] * i1) +
                          Tppc->size[0] * Tppc->size[1] * b_i) +
                         Tppc->size[0] * Tppc->size[1] * Tppc->size[2] * j];
        }
      }
      loop_ub = Sppc->size[0];
      b_loop_ub = Sppc->size[1];
      b_Sppc_size[0] = Sppc->size[0];
      b_Sppc_size[1] = Sppc->size[1];
      for (i1 = 0; i1 < b_loop_ub; i1++) {
        for (i2 = 0; i2 < loop_ub; i2++) {
          b_Sppc_data[i2 + b_Sppc_size[0] * i1] =
              Sppc->data[((i2 + Sppc->size[0] * i1) +
                          Sppc->size[0] * Sppc->size[1] * vlen) +
                         Sppc->size[0] * Sppc->size[1] * Sppc->size[2] * j];
        }
      }
      loop_ub = Tppc->size[0];
      b_loop_ub = Tppc->size[1];
      Tppc_size_idx_0 = Tppc->size[0];
      for (i1 = 0; i1 < b_loop_ub; i1++) {
        for (i2 = 0; i2 < loop_ub; i2++) {
          b_Tppc_data[i2 + Tppc_size_idx_0 * i1] =
              Tppc->data[((i2 + Tppc->size[0] * i1) +
                          Tppc->size[0] * Tppc->size[1] * vlen) +
                         Tppc->size[0] * Tppc->size[1] * Tppc->size[2] * j];
        }
      }
      sx_data[b_i + sx_size[0] * j] = -ntp_midpoint_to_casts(
          Sppc_data, Sppc_size, Tppc_data, Zm_data, Zm_size, km, b_Sppc_data,
          b_Sppc_size, b_Tppc_data, Zn_data, Zn_size, K_data[vlen + k], zm,
          z_data[vlen + npages], tolz);
      /*  --- NTP with neighbour in j dimension */
      if (Zmat) {
        loop_ub = Z->size[0];
        Zn_size = Z->size[0];
        for (i1 = 0; i1 < loop_ub; i1++) {
          Zn_data[i1] = Z->data[(i1 + Z->size[0] * b_i) +
                                Z->size[0] * Z->size[1] * xpageoffset];
        }
      }
      loop_ub = Sppc->size[0];
      b_loop_ub = Sppc->size[1];
      Sppc_size[0] = Sppc->size[0];
      Sppc_size[1] = Sppc->size[1];
      for (i1 = 0; i1 < b_loop_ub; i1++) {
        for (i2 = 0; i2 < loop_ub; i2++) {
          Sppc_data[i2 + Sppc_size[0] * i1] =
              Sppc->data[((i2 + Sppc->size[0] * i1) +
                          Sppc->size[0] * Sppc->size[1] * b_i) +
                         Sppc->size[0] * Sppc->size[1] * Sppc->size[2] * j];
        }
      }
      loop_ub = Tppc->size[0];
      b_loop_ub = Tppc->size[1];
      Tppc_size_idx_0 = Tppc->size[0];
      for (i1 = 0; i1 < b_loop_ub; i1++) {
        for (i2 = 0; i2 < loop_ub; i2++) {
          Tppc_data[i2 + Tppc_size_idx_0 * i1] =
              Tppc->data[((i2 + Tppc->size[0] * i1) +
                          Tppc->size[0] * Tppc->size[1] * b_i) +
                         Tppc->size[0] * Tppc->size[1] * Tppc->size[2] * j];
        }
      }
      loop_ub = Sppc->size[0];
      b_loop_ub = Sppc->size[1];
      b_Sppc_size[0] = Sppc->size[0];
      b_Sppc_size[1] = Sppc->size[1];
      for (i1 = 0; i1 < b_loop_ub; i1++) {
        for (i2 = 0; i2 < loop_ub; i2++) {
          b_Sppc_data[i2 + b_Sppc_size[0] * i1] =
              Sppc->data[((i2 + Sppc->size[0] * i1) +
                          Sppc->size[0] * Sppc->size[1] * b_i) +
                         Sppc->size[0] * Sppc->size[1] * Sppc->size[2] *
                             xpageoffset];
        }
      }
      loop_ub = Tppc->size[0];
      b_loop_ub = Tppc->size[1];
      Tppc_size_idx_0 = Tppc->size[0];
      for (i1 = 0; i1 < b_loop_ub; i1++) {
        for (i2 = 0; i2 < loop_ub; i2++) {
          b_Tppc_data[i2 + Tppc_size_idx_0 * i1] =
              Tppc->data[((i2 + Tppc->size[0] * i1) +
                          Tppc->size[0] * Tppc->size[1] * b_i) +
                         Tppc->size[0] * Tppc->size[1] * Tppc->size[2] *
                             xpageoffset];
        }
      }
      sy_data[b_i + sy_size[0] * j] = -ntp_midpoint_to_casts(
          Sppc_data, Sppc_size, Tppc_data, Zm_data, Zm_size, km, b_Sppc_data,
          b_Sppc_size, b_Tppc_data, Zn_data, Zn_size,
          K_data[b_i + K_size_idx_0 * xpageoffset], zm,
          z_data[b_i + z_size[0] * xpageoffset], tolz);
    }
    /*  i */
  }
  /*  j */
  /*  Divide by horizontal distances */
  xpageoffset = dx_size[0] * dx_size[1];
  for (i = 0; i < xpageoffset; i++) {
    x_data[i] = (dx_data[i] == 1.0);
  }
  Zmat = true;
  vlen = 1;
  exitg1 = false;
  while ((!exitg1) && (vlen <= xpageoffset)) {
    if (!x_data[vlen - 1]) {
      Zmat = false;
      exitg1 = true;
    } else {
      vlen++;
    }
  }
  if (!Zmat) {
    xpageoffset = dy_size[0] * dy_size[1];
    for (i = 0; i < xpageoffset; i++) {
      x_data[i] = (dy_data[i] == 1.0);
    }
    Zmat = true;
    vlen = 1;
    exitg1 = false;
    while ((!exitg1) && (vlen <= xpageoffset)) {
      if (!x_data[vlen - 1]) {
        Zmat = false;
        exitg1 = true;
      } else {
        vlen++;
      }
    }
    if (!Zmat) {
      loop_ub = (int8_T)z_size[0] * (int8_T)z_size[1];
      for (i = 0; i < loop_ub; i++) {
        sx_data[i] /= dx_data[i];
      }
      loop_ub = (int8_T)z_size[0] * (int8_T)z_size[1];
      for (i = 0; i < loop_ub; i++) {
        sy_data[i] /= dy_data[i];
      }
    }
  }
  emlrtHeapReferenceStackLeaveFcnR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (ntp_slope.c) */
