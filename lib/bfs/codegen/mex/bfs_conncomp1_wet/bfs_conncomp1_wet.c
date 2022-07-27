/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * bfs_conncomp1_wet.c
 *
 * Code generation for function 'bfs_conncomp1_wet'
 *
 */

/* Include files */
#include "bfs_conncomp1_wet.h"
#include "bfs_conncomp1_wet_data.h"
#include "bfs_conncomp1_wet_emxutil.h"
#include "bfs_conncomp1_wet_types.h"
#include "ppc_val2.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"

/* Function Definitions */
void bfs_conncomp1_wet(const emxArray_real_T *Sppc, const emxArray_real_T *Tppc,
                       const emxArray_real_T *P, emxArray_real_T *s,
                       emxArray_real_T *t, emxArray_real_T *p,
                       const emxArray_real_T *ML, real_T TOL_P,
                       const emxArray_real_T *A, const emxArray_real_T *BotK,
                       real_T r, emxArray_real_T *qu, real_T *freshly_wet,
                       real_T *qt)
{
  emxArray_boolean_T *G;
  emxArray_boolean_T *dry;
  real_T Sppc_data[392];
  real_T Tppc_data[392];
  real_T b_Sppc_data[392];
  real_T b_Tppc_data[392];
  real_T Pn_data[50];
  real_T b_d;
  real_T b_r;
  real_T b_s;
  real_T b_s1o2;
  real_T b_t;
  real_T c;
  real_T dxm;
  real_T dxp;
  real_T e;
  real_T fa;
  real_T fb;
  real_T fc;
  real_T lb;
  real_T m;
  real_T n;
  real_T pB;
  real_T s1o2;
  real_T sB;
  real_T tB;
  real_T ub;
  real_T z_tmp;
  int32_T Sppc_size[2];
  int32_T b_Sppc_size[2];
  int32_T Pmat;
  int32_T Pn_size;
  int32_T Sppc_idx_1;
  int32_T Tppc_idx_0;
  int32_T Tppc_idx_1;
  int32_T b_loop_ub;
  int32_T d;
  int32_T exitg1;
  int32_T i;
  int32_T i1;
  int32_T i2;
  int32_T loop_ub;
  int32_T nij;
  int32_T nij_idx_0;
  int32_T qh;
  boolean_T fapos;
  boolean_T fbpos;
  boolean_T guard1 = false;
  emlrtHeapReferenceStackEnterFcnR2012b(emlrtRootTLSGlobal);
  /* BFS_CONNCOMP1_WET  Find one connected component using Breadth First Search,
   */
  /*                    and test neutral tangent plane connections from the
   * perimeter */
  /*  */
  /*  */
  /*  [s, t, p, freshly_wet] = bfs_conncomp1_wet(Sppc, Tppc, P, s, t, p, ML,
   * TOL_P, A, BotK, r, qu) */
  /*  works as bfs_conncomp1.m, with G in that function given by isfinite(p). */
  /*  That is, nodes in the surface are walked from the root node r.  Where an
   */
  /*  invalid node is reached which is not part of the surface but nonetheless
   */
  /*  ocean (BotK > 1 here), a neutral tangent plane calculation is performed */
  /*  from the surface to the cast at this node.  If successful, and if the NTP
   */
  /*  link reaches this neighbouring cast below the mixed layer, then this node
   */
  /*  is added to the surface, the salinity and temperature are interpolated */
  /*  using the piecewise polynomials Sppc and Tppc with knots at P, and the */
  /*  BFS continues.  NTP connections are accurate to a pressure or depth of */
  /*  TOL_P. The input qu is optional, simply to save memory by working */
  /*  in-place.  On output, qu(1:qt) are linear indices to the nodes on the */
  /*  surface (whether or not they are freshly wet), in the order that they */
  /*  were added to the queue.  Note qu(1) == r. */
  /*  */
  /*  */
  /*  --- Input: */
  /*  Sppc [O,K-1,ni,nj]: coefficients for piecewise polynomial for practical */
  /*                    / Absolute Salinity in terms of P */
  /*  Tppc [O,K-1,ni,nj]: coefficients for piecewise polynomial for potential */
  /*                    / Conservative Temperature in terms of P */
  /*  P [K,ni,nj]: knots for the pressure or depth of the casts */
  /*  s [ni,nj]: practical / Absolute salinity on the surface */
  /*  t [ni,nj]: potential / Conservative temperature on the surface */
  /*  p [ni,nj]: pressure or depth on the surface */
  /*  ML [ni,nj]: pressure or depth of the mixed layer. */
  /*  TOL_P [1,1]: tolerance in p for finding neutral connections */
  /*  A [D, ni*nj]: adjacency, where D is the most neighbours possible */
  /*  qu [N,1]: vector to work in-place (optional) */
  /*  BotK [ni,nj]: number of valid data points on each cast */
  /*  r [1, 1]: linear index to the reference cast */
  /*  qu [N,1]: vector to work in-place (optional) */
  /*  */
  /*  */
  /*  --- Output: */
  /*  s [ni,nj]: updated practical/Absolute salinity on the surface */
  /*  t [ni,nj]: updated potential/Conservative temperature on the surface */
  /*  p [ni,nj]: updated pressure [dbar] or depth [m] of the surface */
  /*  freshly_wet [1,1]: number of casts freshly wet */
  /*  qu [N,1]: the nodes visited by the BFS's in order from 1 to qt */
  /*  qt [1,1]: the queue tail index */
  /*  */
  /*  */
  /*  See also BFS_CONNCOMP1, GRID_ADJACENCY */
  /*  */
  /*  */
  /*  Author(s) : Geoff Stanley */
  /*  Email     : g.stanley@unsw.edu.au */
  /*  Email     : geoffstanley@gmail.com */
  nij = p->size[0] * p->size[1];
  Pmat = (((P->size[0] != 1) && (P->size[1] != 1)) || (P->size[2] != 1));
  /*  0 if P is a vector, 1 otherwise */
  *freshly_wet = 0.0;
  /*  count number of casts that wetting adds to the surface */
  if ((qu->size[0] == 0) || (qu->size[1] == 0)) {
    i = qu->size[0] * qu->size[1];
    qu->size[0] = nij;
    qu->size[1] = 1;
    emxEnsureCapacity_real_T(qu, i);
    for (i = 0; i < nij; i++) {
      qu->data[i] = 0.0;
    }
    /*  pre-allocate queue storing linear indices to pixels */
  }
  emxInit_boolean_T(&G, 2, true);
  /*  maximal degree */
  i = G->size[0] * G->size[1];
  G->size[0] = p->size[0];
  G->size[1] = p->size[1];
  emxEnsureCapacity_boolean_T(G, i);
  loop_ub = p->size[0] * p->size[1];
  for (i = 0; i < loop_ub; i++) {
    G->data[i] = muDoubleScalarIsInf(p->data[i]);
  }
  emxInit_boolean_T(&dry, 2, true);
  i = dry->size[0] * dry->size[1];
  dry->size[0] = p->size[0];
  dry->size[1] = p->size[1];
  emxEnsureCapacity_boolean_T(dry, i);
  loop_ub = p->size[0] * p->size[1];
  for (i = 0; i < loop_ub; i++) {
    dry->data[i] = muDoubleScalarIsNaN(p->data[i]);
  }
  loop_ub = G->size[0] * G->size[1];
  for (i = 0; i < loop_ub; i++) {
    G->data[i] = ((!G->data[i]) && (!dry->data[i]));
  }
  /*  Good nodes */
  i = dry->size[0] * dry->size[1];
  dry->size[0] = BotK->size[0];
  dry->size[1] = BotK->size[1];
  emxEnsureCapacity_boolean_T(dry, i);
  loop_ub = BotK->size[0] * BotK->size[1];
  for (i = 0; i < loop_ub; i++) {
    dry->data[i] = ((BotK->data[i] > 1.0) && (!G->data[i]));
  }
  /*  Try wetting only these locations: ocean and not currently in the surface
   */
  /*  Queue Tail */
  qh = 0;
  /*  Queue Head */
  *qt = 1.0;
  /*  Add r to queue */
  qu->data[0] = r;
  G->data[(int32_T)r - 1] = false;
  /*  mark r as discovered */
  loop_ub = P->size[0];
  Pn_size = P->size[0];
  for (i = 0; i < loop_ub; i++) {
    Pn_data[i] = P->data[i];
  }
  while (*qt > qh) {
    qh++;
    /*  advance head of the queue */
    m = qu->data[qh - 1];
    /*  me node; pop from head of queue */
    i = A->size[0];
    for (d = 0; d < i; d++) {
      nij_idx_0 = A->size[0];
      n = A->data[d + nij_idx_0 * ((int32_T)m - 1)];
      /*  neighbour node */
      if (n <= nij) {
        /*  check n is not a neighbour across a non-periodic boundary */
        if (G->data[(int32_T)n - 1]) {
          /*  n is good, and undiscovered */
          (*qt)++;
          /*  Add n to queue */
          qu->data[(int32_T)*qt - 1] = n;
          G->data[(int32_T)n - 1] = false;
          /*  mark n as discovered */
        } else if (dry->data[(int32_T)n - 1]) {
          /*  n is "dry".  Try wetting. */
          if (Pmat != 0) {
            loop_ub = P->size[0];
            nij_idx_0 = P->size[0];
            Pn_size = P->size[0];
            for (i1 = 0; i1 < loop_ub; i1++) {
              Pn_data[i1] = P->data[i1 + nij_idx_0 * ((int32_T)n - 1)];
            }
          }
          loop_ub = Sppc->size[0];
          nij_idx_0 = Sppc->size[0];
          b_loop_ub = Sppc->size[1];
          Sppc_idx_1 = Sppc->size[1];
          Sppc_size[0] = Sppc->size[0];
          Sppc_size[1] = Sppc->size[1];
          for (i1 = 0; i1 < b_loop_ub; i1++) {
            for (i2 = 0; i2 < loop_ub; i2++) {
              Sppc_data[i2 + Sppc_size[0] * i1] =
                  Sppc->data[(i2 + nij_idx_0 * i1) +
                             nij_idx_0 * Sppc_idx_1 * ((int32_T)n - 1)];
            }
          }
          loop_ub = Tppc->size[0];
          Tppc_idx_0 = Tppc->size[0];
          b_loop_ub = Tppc->size[1];
          Tppc_idx_1 = Tppc->size[1];
          nij_idx_0 = Tppc->size[0];
          for (i1 = 0; i1 < b_loop_ub; i1++) {
            for (i2 = 0; i2 < loop_ub; i2++) {
              Tppc_data[i2 + nij_idx_0 * i1] =
                  Tppc->data[(i2 + Tppc_idx_0 * i1) +
                             Tppc_idx_0 * Tppc_idx_1 * ((int32_T)n - 1)];
            }
          }
          sB = s->data[(int32_T)m - 1];
          tB = t->data[(int32_T)m - 1];
          pB = p->data[(int32_T)m - 1];
          /* NTP_BOTTLE_TO_CAST  Find the Neutral Tangent Plane from a bottle to
           * a cast. */
          /*  */
          /*  [p, s, t] = ntp_bottle_to_cast(Sppc, Tppc, P, k, sB, tB, pB, tolp)
           */
          /*  finds a point on the cast with properties (s, t, p) that is
           * neutrally  */
          /*  related to the bottle of (sB, tB, pB), meaning that */
          /*    eos(s, t, p_avg) = eos(sB, tB, p_avg) */
          /*  is approximately solved (there is an exact solution within tolp of
           * p), */
          /*  where eos is the equation of state given by eos.m in MATLAB's
           * path, */
          /*  and   p_avg = (pB + p) / 2 is the average of the bottle's pressure
           */
          /*                             and the pressure on the cast, */
          /*  and   [s,t] are the salinity and temperature on the cast at
           * pressure p */
          /*              given by ppc_val2(P, Sppc, Tppc, p). */
          /*  The cast's hydrographic properties are given by piecewise
           * polynomial */
          /*  interpolants for salinity and temperature as functions of
           * pressure, */
          /*  given by coefficient arrays Sppc and Tppc with knots at P(1:k). If
           * no */
          /*  such (s, t, p) is found, each of (s, t, p) is NaN. */
          /*  */
          /*  For a non-Boussinesq ocean, P, pB, and p are pressure. */
          /*  For a Boussinesq ocean, P, pB, and p are depth. */
          /*  */
          /*  */
          /*  --- Input: */
          /*  Sppc [O, K-1]: coefficients for piecewise polynomial for practical
           * /  */
          /*      Absolute Salinity in terms of Z on the cast.  Here, O is the
           */
          /*      polynomial order and K is the number of data points (including
           * NaN's) */
          /*      on the cast. */
          /*  Tppc [O, K-1]: as above but for potential / Conservative
           * Temperature. */
          /*  P [K, 1]: pressure or depth data points on the cast. */
          /*  k [1, 1]: number of valid (non-NaN) data points on the cast. */
          /*           Specifically, Sppc(:,1:k-1) and Tppc(:,1:k-1) must be
           * non-NaN.  */
          /*  sB [1 , 1]: practical / Absolute salinity of current bottle */
          /*  tB [1 , 1]: potential / Conservative temperature of current bottle
           */
          /*  pB [1 , 1]: pressure or depth of current bottle */
          /*  tolp [1, 1]: tolerance for solving the level of neutral buoyancy
           * (same */
          /*               units as P and pB) */
          /*  */
          /*  Note: physical units for Sppc, Tppc, P, sB, tB, pB, p, s, t are */
          /*  determined by eos.m. */
          /*  */
          /*  Note: P must increase monotonically along its first dimension. */
          /*  */
          /*  */
          /*  --- Output: */
          /*  p [1, 1]: pressure or depth in the cast  */
          /*            that is neutrally related to the bottle */
          /*  s [1, 1]: practical / Absolute salinity in the cast */
          /*            that is neutrally related to the bottle */
          /*  t [1, 1]: potential / Conservative temperature in the cast */
          /*            that is neutrally related to the bottle */
          /*  */
          /*  */
          /*  --- See Also: */
          /*  ntp_midpoint_to_casts */
          /*  ppc_linterp, ppc_pchip */
          /*  Author(s) : Geoff Stanley */
          /*  Email     : g.stanley@unsw.edu.au */
          /*  Email     : geoffstanley@gmail.com */
          if (BotK->data[(int32_T)n - 1] > 1.0) {
            /*  Search for a sign-change, expanding outward from an initial
             * guess */
            fa = Pn_data[0];
            c = Pn_data[(int32_T)BotK->data[(int32_T)n - 1] - 1];
            /* FZERO_GUESS_TO_BOUNDS  Search for a sign change bounding a zero
             * of a */
            /*                        univariate function, expanding
             * geometrically */
            /*                        outward from an initial guess. */
            /*  */
            /*  */
            /*  [a, b] = fzero_guess_to_bounds(f, x) */
            /*  finds a < b such that f(a) and f(b) have different sign*,
             * meaning a */
            /*  solution exists within the interval [a,b].  The bounds a,b are
             * expanded */
            /*  outward in geometric progression from an initial guess for the
             * root of f */
            /*  at x. If f evaluates to NaN at any point during the search, then
             * a = nan */
            /*  and b = nan are immediately returned.  If the function is
             * genuinely */
            /*  single-signed, or even if it is not but its values of opposite
             * sign are */
            /*  skipped over, it is possible to enter an infinite loop.  Calling
             * the */
            /*  function in this form is therefore not recommended unless you
             * know the */
            /*  function will not result in such an infinite loop. */
            /*  */
            /*  [a, b] = fzero_guess_to_bounds(f, x, A, B) */
            /*  as above, but limits [a,b] to remain inside the subset [A, B].
             * If x is */
            /*  outside of [A, B], it is immediately moved into this range. If
             * no */
            /*  sign-change is found within [A, B], then a = nan and b = nan are
             */
            /*  returned.  Note, as above, it is possible that a sign-change is
             * skipped */
            /*  over as f is only evaluated at finitely many x values. */
            /*  */
            /*  [a,b] = fzero_guess_to_bounds(f, x, A, B, ...) */
            /*  passes all additional arguments to the function f.  */
            /*  */
            /*  * Note: for computational speed, herein the "sign" of 0 is
             * considered the */
            /*  same as the sign of a negative number. */
            /*  */
            /*  This function is compatible with MATLAB's code generation. */
            /*  */
            /*  */
            /*  --- Input: */
            /*    f       : handle to a function that accepts a real scalar as
             * its first */
            /*              input and returns a real scalar */
            /*    x       : initial scalar guess for a root of f */
            /*    A       : scalar lower bound */
            /*    B       : scalar upper bound */
            /*    varargin: All additional inputs are passed directly to f */
            /*  */
            /*  */
            /*  --- Output: */
            /*    a : lower bound for interval containing a root, scalar */
            /*    b : upper bound for interval containing a root, scalar */
            /*  */
            /*  */
            /*  --- Acknowledgements: */
            /*  Expansion from initial guess inspired by MATLAB's fzero.m. */
            /*  */
            /*  */
            /*  Author    : Geoff Stanley */
            /*  Email     : geoffstanley@gmail.com */
            /*  Version   : 1.0 */
            /*  History   : 01/07/2020 - initial release */
            /*            : 22/07/2020 - fix infinite loop in bounded case,
             * arising from machine precision rounding */
            /*  Geometrically expand from the guess x, until a sign change is
             * found */
            /*  Handle bad inputs */
            fapos = muDoubleScalarIsNaN(Pn_data[0]);
            if (fapos || muDoubleScalarIsNaN(c) || muDoubleScalarIsNaN(pB)) {
              lb = rtNaN;
              ub = rtNaN;
            } else {
              /*  Evaluate difference between (a) eos at location on the cast
               * (S, T, P) */
              /*  where the pressure or depth is p, and (b) eos of the bottle
               * (sB, tB, pB); */
              /*  here, eos is always evaluated at the average pressure or
               * depth, (p + */
              /*  pB)/2. */
              ppc_val2(Pn_data, Pn_size, Sppc_data, Sppc_size, Tppc_data, fa,
                       &b_s, &b_t);
              /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
              /*  */
              /*  */
              /*  rho = eoscg_densjmd95_bsq(s,t,z)                             [
               * kg / m^3 ] */
              /*  computes the Boussinesq JMD95 in-situ density, given practical
               * salinity */
              /*  s, potential temperature t, and depth z [m, positive and
               * increasing */
              /*  down].  The depth z is converted into hydrostatic pressure for
               * a given */
              /*  gravitational acceleration and Boussinesq reference density,
               * which are */
              /*  hard-coded into this function (edit variables grav and rhob as
               * needed). */
              /*  */
              /*  */
              /*  This function is derived from densjmd95.m, documented below.
               * Input checks */
              /*  and expansion of variables have been removed; instead
               * automatic expansion */
              /*  (requiring MATLAB 2016b or later) is used. The calculation has
               * also been */
              /*  streamlined by pre-allocating arrays for coeffcients, and
               * modifying the */
              /*  coefficients such that pressure need not be internally
               * converted from */
              /*  dbar to bar. This function is compatible with MATLAB's
               * codegen. */
              /*  */
              /*  The units of s and t are as in densjmd95, documented below. */
              /*   */
              /*  The inputs' sizes are more general: they may be arrays of any
               * dimension */
              /*  and any size, so long as their sizes match each other
               * excepting that any */
              /*  input can have any dimension as a singleton. For example, s
               * and t can be */
              /*  3D arrays of size [nz, nx, ny], while p can be a vector of
               * size [nz,1]. */
              /*  */
              /*  */
              /*  Author(s)       : Geoff Stanley */
              /*  Email           : g.stanley@unsw.edu.au */
              /*  Email           : geoffstanley@gmail.com */
              /*  Version         : 1.0 */
              /*  */
              /*  */
              /*  DENSJMD95    Density of sea water */
              /* =========================================================================
               */
              /*  */
              /*  USAGE:  dens = densjmd95(S,Theta,P) */
              /*  */
              /*  DESCRIPTION: */
              /*     Density of Sea Water using Jackett and McDougall 1995 (JAOT
               * 12) */
              /*     polynomial (modified UNESCO polynomial). */
              /*  */
              /*  INPUT:  (all must have same dimensions) */
              /*    S     = salinity    [psu      (PSS-78)] */
              /*    Theta = potential temperature [degree C (IPTS-68)] */
              /*    P     = pressure    [dbar] */
              /*        (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) ) */
              /*  */
              /*  OUTPUT: */
              /*    dens = density  [kg/m^3] */
              /*  */
              /*  AUTHOR:  Martin Losch 2002-08-09  (mlosch@mit.edu) */
              /*  Jackett and McDougall, 1995, JAOT 12(4), pp. 381-388 */
              /*  created by mlosch on 2002-08-09 */
              /*  $Header: /u/gcmpack/MITgcm/utils/matlab/densjmd95.m,v 1.2
               * 2007/02/17 23:49:43 jmc Exp $ */
              /*  $Name:  $ */
              /*  */
              /*  Convert depth [m] to the hydrostatic pressure [dbar] implied
               * by the */
              /*  following two hard-coded parameters (edit as needed): */
              /*  gravitational acceleration [m /s^2] */
              /*  Boussinesq reference density [kg / m^3] */
              /*  Pascal to dbar conversion [dbar / Pa] */
              /*  depth to pressure conversion [dbar / m] */
              z_tmp = (pB + fa) / 2.0 * 1.015335;
              /*  Henceforth z is actually pressure [dbar] */
              /*  coefficients nonlinear equation of state in pressure
               * coordinates for */
              /*  1. density of fresh water at p = 0 */
              /*  2. density of sea water at p = 0 */
              /*  coefficients in pressure coordinates for */
              /*  3. secant bulk modulus K of fresh water at p = 0 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  4. secant bulk modulus K of sea water at p = 0 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  5. secant bulk modulus K of sea water at p */
              /*  == original / 10 */
              /*  == original / 10 */
              /*  == original / 10 */
              /*  == original / 10 */
              /*  == original / 10 */
              /*  == original / 10 */
              s1o2 = muDoubleScalarSqrt(sB);
              /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
              /*  */
              /*  */
              /*  rho = eoscg_densjmd95_bsq(s,t,z)                             [
               * kg / m^3 ] */
              /*  computes the Boussinesq JMD95 in-situ density, given practical
               * salinity */
              /*  s, potential temperature t, and depth z [m, positive and
               * increasing */
              /*  down].  The depth z is converted into hydrostatic pressure for
               * a given */
              /*  gravitational acceleration and Boussinesq reference density,
               * which are */
              /*  hard-coded into this function (edit variables grav and rhob as
               * needed). */
              /*  */
              /*  */
              /*  This function is derived from densjmd95.m, documented below.
               * Input checks */
              /*  and expansion of variables have been removed; instead
               * automatic expansion */
              /*  (requiring MATLAB 2016b or later) is used. The calculation has
               * also been */
              /*  streamlined by pre-allocating arrays for coeffcients, and
               * modifying the */
              /*  coefficients such that pressure need not be internally
               * converted from */
              /*  dbar to bar. This function is compatible with MATLAB's
               * codegen. */
              /*  */
              /*  The units of s and t are as in densjmd95, documented below. */
              /*   */
              /*  The inputs' sizes are more general: they may be arrays of any
               * dimension */
              /*  and any size, so long as their sizes match each other
               * excepting that any */
              /*  input can have any dimension as a singleton. For example, s
               * and t can be */
              /*  3D arrays of size [nz, nx, ny], while p can be a vector of
               * size [nz,1]. */
              /*  */
              /*  */
              /*  Author(s)       : Geoff Stanley */
              /*  Email           : g.stanley@unsw.edu.au */
              /*  Email           : geoffstanley@gmail.com */
              /*  Version         : 1.0 */
              /*  */
              /*  */
              /*  DENSJMD95    Density of sea water */
              /* =========================================================================
               */
              /*  */
              /*  USAGE:  dens = densjmd95(S,Theta,P) */
              /*  */
              /*  DESCRIPTION: */
              /*     Density of Sea Water using Jackett and McDougall 1995 (JAOT
               * 12) */
              /*     polynomial (modified UNESCO polynomial). */
              /*  */
              /*  INPUT:  (all must have same dimensions) */
              /*    S     = salinity    [psu      (PSS-78)] */
              /*    Theta = potential temperature [degree C (IPTS-68)] */
              /*    P     = pressure    [dbar] */
              /*        (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) ) */
              /*  */
              /*  OUTPUT: */
              /*    dens = density  [kg/m^3] */
              /*  */
              /*  AUTHOR:  Martin Losch 2002-08-09  (mlosch@mit.edu) */
              /*  Jackett and McDougall, 1995, JAOT 12(4), pp. 381-388 */
              /*  created by mlosch on 2002-08-09 */
              /*  $Header: /u/gcmpack/MITgcm/utils/matlab/densjmd95.m,v 1.2
               * 2007/02/17 23:49:43 jmc Exp $ */
              /*  $Name:  $ */
              /*  */
              /*  Convert depth [m] to the hydrostatic pressure [dbar] implied
               * by the */
              /*  following two hard-coded parameters (edit as needed): */
              /*  gravitational acceleration [m /s^2] */
              /*  Boussinesq reference density [kg / m^3] */
              /*  Pascal to dbar conversion [dbar / Pa] */
              /*  depth to pressure conversion [dbar / m] */
              /*  Henceforth z is actually pressure [dbar] */
              /*  coefficients nonlinear equation of state in pressure
               * coordinates for */
              /*  1. density of fresh water at p = 0 */
              /*  2. density of sea water at p = 0 */
              /*  coefficients in pressure coordinates for */
              /*  3. secant bulk modulus K of fresh water at p = 0 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  4. secant bulk modulus K of sea water at p = 0 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  5. secant bulk modulus K of sea water at p */
              /*  == original / 10 */
              /*  == original / 10 */
              /*  == original / 10 */
              /*  == original / 10 */
              /*  == original / 10 */
              /*  == original / 10 */
              b_s1o2 = muDoubleScalarSqrt(b_s);
              fc = (tB * (tB * 1.39468E-8 + -1.202016E-6) + 2.102898E-5) +
                   sB * (tB * (tB * 6.207323E-11 + 6.128773E-9) + -2.040237E-7);
              e = (tB * (tB * (tB * 1.956415E-6 + -0.0002984642) + 0.02212276) +
                   3.186519) +
                  sB *
                      ((tB * (tB * 2.059331E-7 + -0.0001847318) + 0.006704388) +
                       s1o2 * 0.0001480266);
              b_d = (tB * (tB * (tB * (tB * -0.0004190253 + 0.09648704) +
                                 -17.06103) +
                           1444.304) +
                     196593.3) +
                    sB * ((tB * (tB * (tB * -0.0005084188 + 0.06283263) +
                                 -3.101089) +
                           528.4855) +
                          s1o2 * (tB * (tB * -0.004619924 + 0.09085835) +
                                  3.88664));
              s1o2 = (tB * (tB * (tB * (tB * (tB * 6.536332E-9 + -1.120083E-6) +
                                        0.0001001685) +
                                  -0.00909529) +
                            0.06793952) +
                      999.842594) +
                     sB * (((tB * (tB * (tB * (tB * 5.3875E-9 + -8.2467E-7) +
                                         7.6438E-5) +
                                   -0.0040899) +
                             0.824493) +
                            s1o2 * (tB * (tB * -1.6546E-6 + 0.00010227) +
                                    -0.00572466)) +
                           sB * 0.00048314);
              if (s1o2 / (1.0 - z_tmp / (b_d + z_tmp * (e + z_tmp * fc))) -
                      ((b_t * (b_t * (b_t * (b_t * (b_t * 6.536332E-9 +
                                                    -1.120083E-6) +
                                             0.0001001685) +
                                      -0.00909529) +
                               0.06793952) +
                        999.842594) +
                       b_s * (((b_t * (b_t * (b_t * (b_t * 5.3875E-9 +
                                                     -8.2467E-7) +
                                              7.6438E-5) +
                                       -0.0040899) +
                                0.824493) +
                               b_s1o2 * (b_t * (b_t * -1.6546E-6 + 0.00010227) +
                                         -0.00572466)) +
                              b_s * 0.00048314)) /
                          (1.0 -
                           z_tmp /
                               (((b_t * (b_t * (b_t * (b_t * -0.0004190253 +
                                                       0.09648704) +
                                                -17.06103) +
                                         1444.304) +
                                  196593.3) +
                                 b_s * ((b_t * (b_t * (b_t * -0.0005084188 +
                                                       0.06283263) +
                                                -3.101089) +
                                         528.4855) +
                                        b_s1o2 * (b_t * (b_t * -0.004619924 +
                                                         0.09085835) +
                                                  3.88664))) +
                                z_tmp *
                                    (((b_t * (b_t * (b_t * 1.956415E-6 +
                                                     -0.0002984642) +
                                              0.02212276) +
                                       3.186519) +
                                      b_s * ((b_t * (b_t * 2.059331E-7 +
                                                     -0.0001847318) +
                                              0.006704388) +
                                             b_s1o2 * 0.0001480266)) +
                                     z_tmp * ((b_t * (b_t * 1.39468E-8 +
                                                      -1.202016E-6) +
                                               2.102898E-5) +
                                              b_s * (b_t * (b_t * 6.207323E-11 +
                                                            6.128773E-9) +
                                                     -2.040237E-7))))) ==
                  0.0) {
                lb = fa;
                ub = fa;
              } else {
                /*  Evaluate difference between (a) eos at location on the cast
                 * (S, T, P) */
                /*  where the pressure or depth is p, and (b) eos of the bottle
                 * (sB, tB, pB); */
                /*  here, eos is always evaluated at the average pressure or
                 * depth, (p + */
                /*  pB)/2. */
                ppc_val2(Pn_data, Pn_size, Sppc_data, Sppc_size, Tppc_data, c,
                         &b_s, &b_t);
                /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
                /*  */
                /*  */
                /*  rho = eoscg_densjmd95_bsq(s,t,z) [ kg / m^3 ] */
                /*  computes the Boussinesq JMD95 in-situ density, given
                 * practical salinity */
                /*  s, potential temperature t, and depth z [m, positive and
                 * increasing */
                /*  down].  The depth z is converted into hydrostatic pressure
                 * for a given */
                /*  gravitational acceleration and Boussinesq reference density,
                 * which are */
                /*  hard-coded into this function (edit variables grav and rhob
                 * as needed). */
                /*  */
                /*  */
                /*  This function is derived from densjmd95.m, documented below.
                 * Input checks */
                /*  and expansion of variables have been removed; instead
                 * automatic expansion */
                /*  (requiring MATLAB 2016b or later) is used. The calculation
                 * has also been */
                /*  streamlined by pre-allocating arrays for coeffcients, and
                 * modifying the */
                /*  coefficients such that pressure need not be internally
                 * converted from */
                /*  dbar to bar. This function is compatible with MATLAB's
                 * codegen. */
                /*  */
                /*  The units of s and t are as in densjmd95, documented below.
                 */
                /*   */
                /*  The inputs' sizes are more general: they may be arrays of
                 * any dimension */
                /*  and any size, so long as their sizes match each other
                 * excepting that any */
                /*  input can have any dimension as a singleton. For example, s
                 * and t can be */
                /*  3D arrays of size [nz, nx, ny], while p can be a vector of
                 * size [nz,1]. */
                /*  */
                /*  */
                /*  Author(s)       : Geoff Stanley */
                /*  Email           : g.stanley@unsw.edu.au */
                /*  Email           : geoffstanley@gmail.com */
                /*  Version         : 1.0 */
                /*  */
                /*  */
                /*  DENSJMD95    Density of sea water */
                /* =========================================================================
                 */
                /*  */
                /*  USAGE:  dens = densjmd95(S,Theta,P) */
                /*  */
                /*  DESCRIPTION: */
                /*     Density of Sea Water using Jackett and McDougall 1995
                 * (JAOT 12) */
                /*     polynomial (modified UNESCO polynomial). */
                /*  */
                /*  INPUT:  (all must have same dimensions) */
                /*    S     = salinity    [psu      (PSS-78)] */
                /*    Theta = potential temperature [degree C (IPTS-68)] */
                /*    P     = pressure    [dbar] */
                /*        (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) ) */
                /*  */
                /*  OUTPUT: */
                /*    dens = density  [kg/m^3] */
                /*  */
                /*  AUTHOR:  Martin Losch 2002-08-09  (mlosch@mit.edu) */
                /*  Jackett and McDougall, 1995, JAOT 12(4), pp. 381-388 */
                /*  created by mlosch on 2002-08-09 */
                /*  $Header: /u/gcmpack/MITgcm/utils/matlab/densjmd95.m,v 1.2
                 * 2007/02/17 23:49:43 jmc Exp $ */
                /*  $Name:  $ */
                /*  */
                /*  Convert depth [m] to the hydrostatic pressure [dbar] implied
                 * by the */
                /*  following two hard-coded parameters (edit as needed): */
                /*  gravitational acceleration [m /s^2] */
                /*  Boussinesq reference density [kg / m^3] */
                /*  Pascal to dbar conversion [dbar / Pa] */
                /*  depth to pressure conversion [dbar / m] */
                b_r = (pB + c) / 2.0 * 1.015335;
                /*  Henceforth z is actually pressure [dbar] */
                /*  coefficients nonlinear equation of state in pressure
                 * coordinates for */
                /*  1. density of fresh water at p = 0 */
                /*  2. density of sea water at p = 0 */
                /*  coefficients in pressure coordinates for */
                /*  3. secant bulk modulus K of fresh water at p = 0 */
                /*  == original * 10 */
                /*  == original * 10 */
                /*  == original * 10 */
                /*  == original * 10 */
                /*  == original * 10 */
                /*  4. secant bulk modulus K of sea water at p = 0 */
                /*  == original * 10 */
                /*  == original * 10 */
                /*  == original * 10 */
                /*  == original * 10 */
                /*  == original * 10 */
                /*  == original * 10 */
                /*  == original * 10 */
                /*  5. secant bulk modulus K of sea water at p */
                /*  == original / 10 */
                /*  == original / 10 */
                /*  == original / 10 */
                /*  == original / 10 */
                /*  == original / 10 */
                /*  == original / 10 */
                /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
                /*  */
                /*  */
                /*  rho = eoscg_densjmd95_bsq(s,t,z) [ kg / m^3 ] */
                /*  computes the Boussinesq JMD95 in-situ density, given
                 * practical salinity */
                /*  s, potential temperature t, and depth z [m, positive and
                 * increasing */
                /*  down].  The depth z is converted into hydrostatic pressure
                 * for a given */
                /*  gravitational acceleration and Boussinesq reference density,
                 * which are */
                /*  hard-coded into this function (edit variables grav and rhob
                 * as needed). */
                /*  */
                /*  */
                /*  This function is derived from densjmd95.m, documented below.
                 * Input checks */
                /*  and expansion of variables have been removed; instead
                 * automatic expansion */
                /*  (requiring MATLAB 2016b or later) is used. The calculation
                 * has also been */
                /*  streamlined by pre-allocating arrays for coeffcients, and
                 * modifying the */
                /*  coefficients such that pressure need not be internally
                 * converted from */
                /*  dbar to bar. This function is compatible with MATLAB's
                 * codegen. */
                /*  */
                /*  The units of s and t are as in densjmd95, documented below.
                 */
                /*   */
                /*  The inputs' sizes are more general: they may be arrays of
                 * any dimension */
                /*  and any size, so long as their sizes match each other
                 * excepting that any */
                /*  input can have any dimension as a singleton. For example, s
                 * and t can be */
                /*  3D arrays of size [nz, nx, ny], while p can be a vector of
                 * size [nz,1]. */
                /*  */
                /*  */
                /*  Author(s)       : Geoff Stanley */
                /*  Email           : g.stanley@unsw.edu.au */
                /*  Email           : geoffstanley@gmail.com */
                /*  Version         : 1.0 */
                /*  */
                /*  */
                /*  DENSJMD95    Density of sea water */
                /* =========================================================================
                 */
                /*  */
                /*  USAGE:  dens = densjmd95(S,Theta,P) */
                /*  */
                /*  DESCRIPTION: */
                /*     Density of Sea Water using Jackett and McDougall 1995
                 * (JAOT 12) */
                /*     polynomial (modified UNESCO polynomial). */
                /*  */
                /*  INPUT:  (all must have same dimensions) */
                /*    S     = salinity    [psu      (PSS-78)] */
                /*    Theta = potential temperature [degree C (IPTS-68)] */
                /*    P     = pressure    [dbar] */
                /*        (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) ) */
                /*  */
                /*  OUTPUT: */
                /*    dens = density  [kg/m^3] */
                /*  */
                /*  AUTHOR:  Martin Losch 2002-08-09  (mlosch@mit.edu) */
                /*  Jackett and McDougall, 1995, JAOT 12(4), pp. 381-388 */
                /*  created by mlosch on 2002-08-09 */
                /*  $Header: /u/gcmpack/MITgcm/utils/matlab/densjmd95.m,v 1.2
                 * 2007/02/17 23:49:43 jmc Exp $ */
                /*  $Name:  $ */
                /*  */
                /*  Convert depth [m] to the hydrostatic pressure [dbar] implied
                 * by the */
                /*  following two hard-coded parameters (edit as needed): */
                /*  gravitational acceleration [m /s^2] */
                /*  Boussinesq reference density [kg / m^3] */
                /*  Pascal to dbar conversion [dbar / Pa] */
                /*  depth to pressure conversion [dbar / m] */
                /*  Henceforth z is actually pressure [dbar] */
                /*  coefficients nonlinear equation of state in pressure
                 * coordinates for */
                /*  1. density of fresh water at p = 0 */
                /*  2. density of sea water at p = 0 */
                /*  coefficients in pressure coordinates for */
                /*  3. secant bulk modulus K of fresh water at p = 0 */
                /*  == original * 10 */
                /*  == original * 10 */
                /*  == original * 10 */
                /*  == original * 10 */
                /*  == original * 10 */
                /*  4. secant bulk modulus K of sea water at p = 0 */
                /*  == original * 10 */
                /*  == original * 10 */
                /*  == original * 10 */
                /*  == original * 10 */
                /*  == original * 10 */
                /*  == original * 10 */
                /*  == original * 10 */
                /*  5. secant bulk modulus K of sea water at p */
                /*  == original / 10 */
                /*  == original / 10 */
                /*  == original / 10 */
                /*  == original / 10 */
                /*  == original / 10 */
                /*  == original / 10 */
                b_s1o2 = muDoubleScalarSqrt(b_s);
                if (s1o2 / (1.0 - b_r / (b_d + b_r * (e + b_r * fc))) -
                        ((b_t * (b_t * (b_t * (b_t * (b_t * 6.536332E-9 +
                                                      -1.120083E-6) +
                                               0.0001001685) +
                                        -0.00909529) +
                                 0.06793952) +
                          999.842594) +
                         b_s *
                             (((b_t * (b_t * (b_t * (b_t * 5.3875E-9 +
                                                     -8.2467E-7) +
                                              7.6438E-5) +
                                       -0.0040899) +
                                0.824493) +
                               b_s1o2 * (b_t * (b_t * -1.6546E-6 + 0.00010227) +
                                         -0.00572466)) +
                              b_s * 0.00048314)) /
                            (1.0 -
                             b_r /
                                 (((b_t * (b_t * (b_t * (b_t * -0.0004190253 +
                                                         0.09648704) +
                                                  -17.06103) +
                                           1444.304) +
                                    196593.3) +
                                   b_s * ((b_t * (b_t * (b_t * -0.0005084188 +
                                                         0.06283263) +
                                                  -3.101089) +
                                           528.4855) +
                                          b_s1o2 * (b_t * (b_t * -0.004619924 +
                                                           0.09085835) +
                                                    3.88664))) +
                                  b_r *
                                      (((b_t * (b_t * (b_t * 1.956415E-6 +
                                                       -0.0002984642) +
                                                0.02212276) +
                                         3.186519) +
                                        b_s * ((b_t * (b_t * 2.059331E-7 +
                                                       -0.0001847318) +
                                                0.006704388) +
                                               b_s1o2 * 0.0001480266)) +
                                       b_r * ((b_t * (b_t * 1.39468E-8 +
                                                      -1.202016E-6) +
                                               2.102898E-5) +
                                              b_s * (b_t * (b_t * 6.207323E-11 +
                                                            6.128773E-9) +
                                                     -2.040237E-7))))) ==
                    0.0) {
                  lb = c;
                  ub = c;
                } else {
                  fb = muDoubleScalarMin(muDoubleScalarMax(pB, fa), c);
                  /*  bounds are given */
                  dxp = (c - fb) / 50.0;
                  dxm = (fb - fa) / 50.0;
                  /*  Set a = x, except when x is so close to A that machine
                   * roundoff makes dxm identically 0, */
                  /*  which would lead to an infinite loop below.  In this case,
                   * set a = A. */
                  if (dxm == 0.0) {
                    lb = fa;
                  } else {
                    lb = fb;
                  }
                  /*  Evaluate difference between (a) eos at location on the
                   * cast (S, T, P) */
                  /*  where the pressure or depth is p, and (b) eos of the
                   * bottle (sB, tB, pB); */
                  /*  here, eos is always evaluated at the average pressure or
                   * depth, (p + */
                  /*  pB)/2. */
                  ppc_val2(Pn_data, Pn_size, Sppc_data, Sppc_size, Tppc_data,
                           lb, &b_s, &b_t);
                  /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density.
                   */
                  /*  */
                  /*  */
                  /*  rho = eoscg_densjmd95_bsq(s,t,z) [ kg / m^3 ] */
                  /*  computes the Boussinesq JMD95 in-situ density, given
                   * practical salinity */
                  /*  s, potential temperature t, and depth z [m, positive and
                   * increasing */
                  /*  down].  The depth z is converted into hydrostatic pressure
                   * for a given */
                  /*  gravitational acceleration and Boussinesq reference
                   * density, which are */
                  /*  hard-coded into this function (edit variables grav and
                   * rhob as needed). */
                  /*  */
                  /*  */
                  /*  This function is derived from densjmd95.m, documented
                   * below. Input checks */
                  /*  and expansion of variables have been removed; instead
                   * automatic expansion */
                  /*  (requiring MATLAB 2016b or later) is used. The calculation
                   * has also been */
                  /*  streamlined by pre-allocating arrays for coeffcients, and
                   * modifying the */
                  /*  coefficients such that pressure need not be internally
                   * converted from */
                  /*  dbar to bar. This function is compatible with MATLAB's
                   * codegen. */
                  /*  */
                  /*  The units of s and t are as in densjmd95, documented
                   * below. */
                  /*   */
                  /*  The inputs' sizes are more general: they may be arrays of
                   * any dimension */
                  /*  and any size, so long as their sizes match each other
                   * excepting that any */
                  /*  input can have any dimension as a singleton. For example,
                   * s and t can be */
                  /*  3D arrays of size [nz, nx, ny], while p can be a vector of
                   * size [nz,1]. */
                  /*  */
                  /*  */
                  /*  Author(s)       : Geoff Stanley */
                  /*  Email           : g.stanley@unsw.edu.au */
                  /*  Email           : geoffstanley@gmail.com */
                  /*  Version         : 1.0 */
                  /*  */
                  /*  */
                  /*  DENSJMD95    Density of sea water */
                  /* =========================================================================
                   */
                  /*  */
                  /*  USAGE:  dens = densjmd95(S,Theta,P) */
                  /*  */
                  /*  DESCRIPTION: */
                  /*     Density of Sea Water using Jackett and McDougall 1995
                   * (JAOT 12) */
                  /*     polynomial (modified UNESCO polynomial). */
                  /*  */
                  /*  INPUT:  (all must have same dimensions) */
                  /*    S     = salinity    [psu      (PSS-78)] */
                  /*    Theta = potential temperature [degree C (IPTS-68)] */
                  /*    P     = pressure    [dbar] */
                  /*        (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) )
                   */
                  /*  */
                  /*  OUTPUT: */
                  /*    dens = density  [kg/m^3] */
                  /*  */
                  /*  AUTHOR:  Martin Losch 2002-08-09  (mlosch@mit.edu) */
                  /*  Jackett and McDougall, 1995, JAOT 12(4), pp. 381-388 */
                  /*  created by mlosch on 2002-08-09 */
                  /*  $Header: /u/gcmpack/MITgcm/utils/matlab/densjmd95.m,v 1.2
                   * 2007/02/17 23:49:43 jmc Exp $ */
                  /*  $Name:  $ */
                  /*  */
                  /*  Convert depth [m] to the hydrostatic pressure [dbar]
                   * implied by the */
                  /*  following two hard-coded parameters (edit as needed): */
                  /*  gravitational acceleration [m /s^2] */
                  /*  Boussinesq reference density [kg / m^3] */
                  /*  Pascal to dbar conversion [dbar / Pa] */
                  /*  depth to pressure conversion [dbar / m] */
                  z_tmp = (pB + lb) / 2.0 * 1.015335;
                  /*  Henceforth z is actually pressure [dbar] */
                  /*  coefficients nonlinear equation of state in pressure
                   * coordinates for */
                  /*  1. density of fresh water at p = 0 */
                  /*  2. density of sea water at p = 0 */
                  /*  coefficients in pressure coordinates for */
                  /*  3. secant bulk modulus K of fresh water at p = 0 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  4. secant bulk modulus K of sea water at p = 0 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  5. secant bulk modulus K of sea water at p */
                  /*  == original / 10 */
                  /*  == original / 10 */
                  /*  == original / 10 */
                  /*  == original / 10 */
                  /*  == original / 10 */
                  /*  == original / 10 */
                  /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density.
                   */
                  /*  */
                  /*  */
                  /*  rho = eoscg_densjmd95_bsq(s,t,z) [ kg / m^3 ] */
                  /*  computes the Boussinesq JMD95 in-situ density, given
                   * practical salinity */
                  /*  s, potential temperature t, and depth z [m, positive and
                   * increasing */
                  /*  down].  The depth z is converted into hydrostatic pressure
                   * for a given */
                  /*  gravitational acceleration and Boussinesq reference
                   * density, which are */
                  /*  hard-coded into this function (edit variables grav and
                   * rhob as needed). */
                  /*  */
                  /*  */
                  /*  This function is derived from densjmd95.m, documented
                   * below. Input checks */
                  /*  and expansion of variables have been removed; instead
                   * automatic expansion */
                  /*  (requiring MATLAB 2016b or later) is used. The calculation
                   * has also been */
                  /*  streamlined by pre-allocating arrays for coeffcients, and
                   * modifying the */
                  /*  coefficients such that pressure need not be internally
                   * converted from */
                  /*  dbar to bar. This function is compatible with MATLAB's
                   * codegen. */
                  /*  */
                  /*  The units of s and t are as in densjmd95, documented
                   * below. */
                  /*   */
                  /*  The inputs' sizes are more general: they may be arrays of
                   * any dimension */
                  /*  and any size, so long as their sizes match each other
                   * excepting that any */
                  /*  input can have any dimension as a singleton. For example,
                   * s and t can be */
                  /*  3D arrays of size [nz, nx, ny], while p can be a vector of
                   * size [nz,1]. */
                  /*  */
                  /*  */
                  /*  Author(s)       : Geoff Stanley */
                  /*  Email           : g.stanley@unsw.edu.au */
                  /*  Email           : geoffstanley@gmail.com */
                  /*  Version         : 1.0 */
                  /*  */
                  /*  */
                  /*  DENSJMD95    Density of sea water */
                  /* =========================================================================
                   */
                  /*  */
                  /*  USAGE:  dens = densjmd95(S,Theta,P) */
                  /*  */
                  /*  DESCRIPTION: */
                  /*     Density of Sea Water using Jackett and McDougall 1995
                   * (JAOT 12) */
                  /*     polynomial (modified UNESCO polynomial). */
                  /*  */
                  /*  INPUT:  (all must have same dimensions) */
                  /*    S     = salinity    [psu      (PSS-78)] */
                  /*    Theta = potential temperature [degree C (IPTS-68)] */
                  /*    P     = pressure    [dbar] */
                  /*        (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) )
                   */
                  /*  */
                  /*  OUTPUT: */
                  /*    dens = density  [kg/m^3] */
                  /*  */
                  /*  AUTHOR:  Martin Losch 2002-08-09  (mlosch@mit.edu) */
                  /*  Jackett and McDougall, 1995, JAOT 12(4), pp. 381-388 */
                  /*  created by mlosch on 2002-08-09 */
                  /*  $Header: /u/gcmpack/MITgcm/utils/matlab/densjmd95.m,v 1.2
                   * 2007/02/17 23:49:43 jmc Exp $ */
                  /*  $Name:  $ */
                  /*  */
                  /*  Convert depth [m] to the hydrostatic pressure [dbar]
                   * implied by the */
                  /*  following two hard-coded parameters (edit as needed): */
                  /*  gravitational acceleration [m /s^2] */
                  /*  Boussinesq reference density [kg / m^3] */
                  /*  Pascal to dbar conversion [dbar / Pa] */
                  /*  depth to pressure conversion [dbar / m] */
                  /*  Henceforth z is actually pressure [dbar] */
                  /*  coefficients nonlinear equation of state in pressure
                   * coordinates for */
                  /*  1. density of fresh water at p = 0 */
                  /*  2. density of sea water at p = 0 */
                  /*  coefficients in pressure coordinates for */
                  /*  3. secant bulk modulus K of fresh water at p = 0 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  4. secant bulk modulus K of sea water at p = 0 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  5. secant bulk modulus K of sea water at p */
                  /*  == original / 10 */
                  /*  == original / 10 */
                  /*  == original / 10 */
                  /*  == original / 10 */
                  /*  == original / 10 */
                  /*  == original / 10 */
                  b_s1o2 = muDoubleScalarSqrt(b_s);
                  fapos =
                      (s1o2 / (1.0 - z_tmp / (b_d + z_tmp * (e + z_tmp * fc))) -
                           ((b_t * (b_t * (b_t * (b_t * (b_t * 6.536332E-9 +
                                                         -1.120083E-6) +
                                                  0.0001001685) +
                                           -0.00909529) +
                                    0.06793952) +
                             999.842594) +
                            b_s * (((b_t * (b_t * (b_t * (b_t * 5.3875E-9 +
                                                          -8.2467E-7) +
                                                   7.6438E-5) +
                                            -0.0040899) +
                                     0.824493) +
                                    b_s1o2 *
                                        (b_t * (b_t * -1.6546E-6 + 0.00010227) +
                                         -0.00572466)) +
                                   b_s * 0.00048314)) /
                               (1.0 -
                                z_tmp /
                                    (((b_t *
                                           (b_t * (b_t * (b_t * -0.0004190253 +
                                                          0.09648704) +
                                                   -17.06103) +
                                            1444.304) +
                                       196593.3) +
                                      b_s *
                                          ((b_t * (b_t * (b_t * -0.0005084188 +
                                                          0.06283263) +
                                                   -3.101089) +
                                            528.4855) +
                                           b_s1o2 * (b_t * (b_t * -0.004619924 +
                                                            0.09085835) +
                                                     3.88664))) +
                                     z_tmp *
                                         (((b_t * (b_t * (b_t * 1.956415E-6 +
                                                          -0.0002984642) +
                                                   0.02212276) +
                                            3.186519) +
                                           b_s * ((b_t * (b_t * 2.059331E-7 +
                                                          -0.0001847318) +
                                                   0.006704388) +
                                                  b_s1o2 * 0.0001480266)) +
                                          z_tmp *
                                              ((b_t * (b_t * 1.39468E-8 +
                                                       -1.202016E-6) +
                                                2.102898E-5) +
                                               b_s *
                                                   (b_t * (b_t * 6.207323E-11 +
                                                           6.128773E-9) +
                                                    -2.040237E-7))))) >
                       0.0);
                  /*  Similarly, set b = x, except for machine precision
                   * problems. */
                  if (dxp == 0.0) {
                    ub = c;
                    /*  Evaluate difference between (a) eos at location on the
                     * cast (S, T, P) */
                    /*  where the pressure or depth is p, and (b) eos of the
                     * bottle (sB, tB, pB); */
                    /*  here, eos is always evaluated at the average pressure or
                     * depth, (p + */
                    /*  pB)/2. */
                    ppc_val2(Pn_data, Pn_size, Sppc_data, Sppc_size, Tppc_data,
                             c, &b_s, &b_t);
                    /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ
                     * density. */
                    /*  */
                    /*  */
                    /*  rho = eoscg_densjmd95_bsq(s,t,z) [ kg / m^3 ] */
                    /*  computes the Boussinesq JMD95 in-situ density, given
                     * practical salinity */
                    /*  s, potential temperature t, and depth z [m, positive and
                     * increasing */
                    /*  down].  The depth z is converted into hydrostatic
                     * pressure for a given */
                    /*  gravitational acceleration and Boussinesq reference
                     * density, which are */
                    /*  hard-coded into this function (edit variables grav and
                     * rhob as needed). */
                    /*  */
                    /*  */
                    /*  This function is derived from densjmd95.m, documented
                     * below. Input checks */
                    /*  and expansion of variables have been removed; instead
                     * automatic expansion */
                    /*  (requiring MATLAB 2016b or later) is used. The
                     * calculation has also been */
                    /*  streamlined by pre-allocating arrays for coeffcients,
                     * and modifying the */
                    /*  coefficients such that pressure need not be internally
                     * converted from */
                    /*  dbar to bar. This function is compatible with MATLAB's
                     * codegen. */
                    /*  */
                    /*  The units of s and t are as in densjmd95, documented
                     * below. */
                    /*   */
                    /*  The inputs' sizes are more general: they may be arrays
                     * of any dimension */
                    /*  and any size, so long as their sizes match each other
                     * excepting that any */
                    /*  input can have any dimension as a singleton. For
                     * example, s and t can be */
                    /*  3D arrays of size [nz, nx, ny], while p can be a vector
                     * of size [nz,1]. */
                    /*  */
                    /*  */
                    /*  Author(s)       : Geoff Stanley */
                    /*  Email           : g.stanley@unsw.edu.au */
                    /*  Email           : geoffstanley@gmail.com */
                    /*  Version         : 1.0 */
                    /*  */
                    /*  */
                    /*  DENSJMD95    Density of sea water */
                    /* =========================================================================
                     */
                    /*  */
                    /*  USAGE:  dens = densjmd95(S,Theta,P) */
                    /*  */
                    /*  DESCRIPTION: */
                    /*     Density of Sea Water using Jackett and McDougall 1995
                     * (JAOT 12) */
                    /*     polynomial (modified UNESCO polynomial). */
                    /*  */
                    /*  INPUT:  (all must have same dimensions) */
                    /*    S     = salinity    [psu      (PSS-78)] */
                    /*    Theta = potential temperature [degree C (IPTS-68)] */
                    /*    P     = pressure    [dbar] */
                    /*        (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) )
                     */
                    /*  */
                    /*  OUTPUT: */
                    /*    dens = density  [kg/m^3] */
                    /*  */
                    /*  AUTHOR:  Martin Losch 2002-08-09  (mlosch@mit.edu) */
                    /*  Jackett and McDougall, 1995, JAOT 12(4), pp. 381-388 */
                    /*  created by mlosch on 2002-08-09 */
                    /*  $Header:
                     * /u/gcmpack/MITgcm/utils/matlab/densjmd95.m,v 1.2
                     * 2007/02/17 23:49:43 jmc Exp $ */
                    /*  $Name:  $ */
                    /*  */
                    /*  Convert depth [m] to the hydrostatic pressure [dbar]
                     * implied by the */
                    /*  following two hard-coded parameters (edit as needed): */
                    /*  gravitational acceleration [m /s^2] */
                    /*  Boussinesq reference density [kg / m^3] */
                    /*  Pascal to dbar conversion [dbar / Pa] */
                    /*  depth to pressure conversion [dbar / m] */
                    /*  Henceforth z is actually pressure [dbar] */
                    /*  coefficients nonlinear equation of state in pressure
                     * coordinates for */
                    /*  1. density of fresh water at p = 0 */
                    /*  2. density of sea water at p = 0 */
                    /*  coefficients in pressure coordinates for */
                    /*  3. secant bulk modulus K of fresh water at p = 0 */
                    /*  == original * 10 */
                    /*  == original * 10 */
                    /*  == original * 10 */
                    /*  == original * 10 */
                    /*  == original * 10 */
                    /*  4. secant bulk modulus K of sea water at p = 0 */
                    /*  == original * 10 */
                    /*  == original * 10 */
                    /*  == original * 10 */
                    /*  == original * 10 */
                    /*  == original * 10 */
                    /*  == original * 10 */
                    /*  == original * 10 */
                    /*  5. secant bulk modulus K of sea water at p */
                    /*  == original / 10 */
                    /*  == original / 10 */
                    /*  == original / 10 */
                    /*  == original / 10 */
                    /*  == original / 10 */
                    /*  == original / 10 */
                    /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ
                     * density. */
                    /*  */
                    /*  */
                    /*  rho = eoscg_densjmd95_bsq(s,t,z) [ kg / m^3 ] */
                    /*  computes the Boussinesq JMD95 in-situ density, given
                     * practical salinity */
                    /*  s, potential temperature t, and depth z [m, positive and
                     * increasing */
                    /*  down].  The depth z is converted into hydrostatic
                     * pressure for a given */
                    /*  gravitational acceleration and Boussinesq reference
                     * density, which are */
                    /*  hard-coded into this function (edit variables grav and
                     * rhob as needed). */
                    /*  */
                    /*  */
                    /*  This function is derived from densjmd95.m, documented
                     * below. Input checks */
                    /*  and expansion of variables have been removed; instead
                     * automatic expansion */
                    /*  (requiring MATLAB 2016b or later) is used. The
                     * calculation has also been */
                    /*  streamlined by pre-allocating arrays for coeffcients,
                     * and modifying the */
                    /*  coefficients such that pressure need not be internally
                     * converted from */
                    /*  dbar to bar. This function is compatible with MATLAB's
                     * codegen. */
                    /*  */
                    /*  The units of s and t are as in densjmd95, documented
                     * below. */
                    /*   */
                    /*  The inputs' sizes are more general: they may be arrays
                     * of any dimension */
                    /*  and any size, so long as their sizes match each other
                     * excepting that any */
                    /*  input can have any dimension as a singleton. For
                     * example, s and t can be */
                    /*  3D arrays of size [nz, nx, ny], while p can be a vector
                     * of size [nz,1]. */
                    /*  */
                    /*  */
                    /*  Author(s)       : Geoff Stanley */
                    /*  Email           : g.stanley@unsw.edu.au */
                    /*  Email           : geoffstanley@gmail.com */
                    /*  Version         : 1.0 */
                    /*  */
                    /*  */
                    /*  DENSJMD95    Density of sea water */
                    /* =========================================================================
                     */
                    /*  */
                    /*  USAGE:  dens = densjmd95(S,Theta,P) */
                    /*  */
                    /*  DESCRIPTION: */
                    /*     Density of Sea Water using Jackett and McDougall 1995
                     * (JAOT 12) */
                    /*     polynomial (modified UNESCO polynomial). */
                    /*  */
                    /*  INPUT:  (all must have same dimensions) */
                    /*    S     = salinity    [psu      (PSS-78)] */
                    /*    Theta = potential temperature [degree C (IPTS-68)] */
                    /*    P     = pressure    [dbar] */
                    /*        (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) )
                     */
                    /*  */
                    /*  OUTPUT: */
                    /*    dens = density  [kg/m^3] */
                    /*  */
                    /*  AUTHOR:  Martin Losch 2002-08-09  (mlosch@mit.edu) */
                    /*  Jackett and McDougall, 1995, JAOT 12(4), pp. 381-388 */
                    /*  created by mlosch on 2002-08-09 */
                    /*  $Header:
                     * /u/gcmpack/MITgcm/utils/matlab/densjmd95.m,v 1.2
                     * 2007/02/17 23:49:43 jmc Exp $ */
                    /*  $Name:  $ */
                    /*  */
                    /*  Convert depth [m] to the hydrostatic pressure [dbar]
                     * implied by the */
                    /*  following two hard-coded parameters (edit as needed): */
                    /*  gravitational acceleration [m /s^2] */
                    /*  Boussinesq reference density [kg / m^3] */
                    /*  Pascal to dbar conversion [dbar / Pa] */
                    /*  depth to pressure conversion [dbar / m] */
                    /*  Henceforth z is actually pressure [dbar] */
                    /*  coefficients nonlinear equation of state in pressure
                     * coordinates for */
                    /*  1. density of fresh water at p = 0 */
                    /*  2. density of sea water at p = 0 */
                    /*  coefficients in pressure coordinates for */
                    /*  3. secant bulk modulus K of fresh water at p = 0 */
                    /*  == original * 10 */
                    /*  == original * 10 */
                    /*  == original * 10 */
                    /*  == original * 10 */
                    /*  == original * 10 */
                    /*  4. secant bulk modulus K of sea water at p = 0 */
                    /*  == original * 10 */
                    /*  == original * 10 */
                    /*  == original * 10 */
                    /*  == original * 10 */
                    /*  == original * 10 */
                    /*  == original * 10 */
                    /*  == original * 10 */
                    /*  5. secant bulk modulus K of sea water at p */
                    /*  == original / 10 */
                    /*  == original / 10 */
                    /*  == original / 10 */
                    /*  == original / 10 */
                    /*  == original / 10 */
                    /*  == original / 10 */
                    b_s1o2 = muDoubleScalarSqrt(b_s);
                    fbpos =
                        (s1o2 / (1.0 - b_r / (b_d + b_r * (e + b_r * fc))) -
                             ((b_t * (b_t * (b_t * (b_t * (b_t * 6.536332E-9 +
                                                           -1.120083E-6) +
                                                    0.0001001685) +
                                             -0.00909529) +
                                      0.06793952) +
                               999.842594) +
                              b_s * (((b_t * (b_t * (b_t * (b_t * 5.3875E-9 +
                                                            -8.2467E-7) +
                                                     7.6438E-5) +
                                              -0.0040899) +
                                       0.824493) +
                                      b_s1o2 * (b_t * (b_t * -1.6546E-6 +
                                                       0.00010227) +
                                                -0.00572466)) +
                                     b_s * 0.00048314)) /
                                 (1.0 -
                                  b_r /
                                      (((b_t *
                                             (b_t *
                                                  (b_t * (b_t * -0.0004190253 +
                                                          0.09648704) +
                                                   -17.06103) +
                                              1444.304) +
                                         196593.3) +
                                        b_s *
                                            ((b_t *
                                                  (b_t * (b_t * -0.0005084188 +
                                                          0.06283263) +
                                                   -3.101089) +
                                              528.4855) +
                                             b_s1o2 *
                                                 (b_t * (b_t * -0.004619924 +
                                                         0.09085835) +
                                                  3.88664))) +
                                       b_r *
                                           (((b_t * (b_t * (b_t * 1.956415E-6 +
                                                            -0.0002984642) +
                                                     0.02212276) +
                                              3.186519) +
                                             b_s * ((b_t * (b_t * 2.059331E-7 +
                                                            -0.0001847318) +
                                                     0.006704388) +
                                                    b_s1o2 * 0.0001480266)) +
                                            b_r *
                                                ((b_t * (b_t * 1.39468E-8 +
                                                         -1.202016E-6) +
                                                  2.102898E-5) +
                                                 b_s *
                                                     (b_t *
                                                          (b_t * 6.207323E-11 +
                                                           6.128773E-9) +
                                                      -2.040237E-7))))) >
                         0.0);
                  } else {
                    ub = fb;
                    if (dxm == 0.0) {
                      fbpos = fapos;
                      /*  since a = b = x */
                    } else {
                      /*  Evaluate difference between (a) eos at location on the
                       * cast (S, T, P) */
                      /*  where the pressure or depth is p, and (b) eos of the
                       * bottle (sB, tB, pB); */
                      /*  here, eos is always evaluated at the average pressure
                       * or depth, (p + */
                      /*  pB)/2. */
                      ppc_val2(Pn_data, Pn_size, Sppc_data, Sppc_size,
                               Tppc_data, fb, &b_s, &b_t);
                      /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ
                       * density. */
                      /*  */
                      /*  */
                      /*  rho = eoscg_densjmd95_bsq(s,t,z) [ kg / m^3 ] */
                      /*  computes the Boussinesq JMD95 in-situ density, given
                       * practical salinity */
                      /*  s, potential temperature t, and depth z [m, positive
                       * and increasing */
                      /*  down].  The depth z is converted into hydrostatic
                       * pressure for a given */
                      /*  gravitational acceleration and Boussinesq reference
                       * density, which are */
                      /*  hard-coded into this function (edit variables grav and
                       * rhob as needed). */
                      /*  */
                      /*  */
                      /*  This function is derived from densjmd95.m, documented
                       * below. Input checks */
                      /*  and expansion of variables have been removed; instead
                       * automatic expansion */
                      /*  (requiring MATLAB 2016b or later) is used. The
                       * calculation has also been */
                      /*  streamlined by pre-allocating arrays for coeffcients,
                       * and modifying the */
                      /*  coefficients such that pressure need not be internally
                       * converted from */
                      /*  dbar to bar. This function is compatible with MATLAB's
                       * codegen. */
                      /*  */
                      /*  The units of s and t are as in densjmd95, documented
                       * below. */
                      /*   */
                      /*  The inputs' sizes are more general: they may be arrays
                       * of any dimension */
                      /*  and any size, so long as their sizes match each other
                       * excepting that any */
                      /*  input can have any dimension as a singleton. For
                       * example, s and t can be */
                      /*  3D arrays of size [nz, nx, ny], while p can be a
                       * vector of size [nz,1]. */
                      /*  */
                      /*  */
                      /*  Author(s)       : Geoff Stanley */
                      /*  Email           : g.stanley@unsw.edu.au */
                      /*  Email           : geoffstanley@gmail.com */
                      /*  Version         : 1.0 */
                      /*  */
                      /*  */
                      /*  DENSJMD95    Density of sea water */
                      /* =========================================================================
                       */
                      /*  */
                      /*  USAGE:  dens = densjmd95(S,Theta,P) */
                      /*  */
                      /*  DESCRIPTION: */
                      /*     Density of Sea Water using Jackett and McDougall
                       * 1995 (JAOT 12) */
                      /*     polynomial (modified UNESCO polynomial). */
                      /*  */
                      /*  INPUT:  (all must have same dimensions) */
                      /*    S     = salinity    [psu      (PSS-78)] */
                      /*    Theta = potential temperature [degree C (IPTS-68)]
                       */
                      /*    P     = pressure    [dbar] */
                      /*        (P may have dims 1x1, mx1, 1xn or mxn for S(mxn)
                       * ) */
                      /*  */
                      /*  OUTPUT: */
                      /*    dens = density  [kg/m^3] */
                      /*  */
                      /*  AUTHOR:  Martin Losch 2002-08-09  (mlosch@mit.edu) */
                      /*  Jackett and McDougall, 1995, JAOT 12(4), pp. 381-388
                       */
                      /*  created by mlosch on 2002-08-09 */
                      /*  $Header:
                       * /u/gcmpack/MITgcm/utils/matlab/densjmd95.m,v 1.2
                       * 2007/02/17 23:49:43 jmc Exp $ */
                      /*  $Name:  $ */
                      /*  */
                      /*  Convert depth [m] to the hydrostatic pressure [dbar]
                       * implied by the */
                      /*  following two hard-coded parameters (edit as needed):
                       */
                      /*  gravitational acceleration [m /s^2] */
                      /*  Boussinesq reference density [kg / m^3] */
                      /*  Pascal to dbar conversion [dbar / Pa] */
                      /*  depth to pressure conversion [dbar / m] */
                      z_tmp = (pB + fb) / 2.0 * 1.015335;
                      /*  Henceforth z is actually pressure [dbar] */
                      /*  coefficients nonlinear equation of state in pressure
                       * coordinates for */
                      /*  1. density of fresh water at p = 0 */
                      /*  2. density of sea water at p = 0 */
                      /*  coefficients in pressure coordinates for */
                      /*  3. secant bulk modulus K of fresh water at p = 0 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  4. secant bulk modulus K of sea water at p = 0 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  5. secant bulk modulus K of sea water at p */
                      /*  == original / 10 */
                      /*  == original / 10 */
                      /*  == original / 10 */
                      /*  == original / 10 */
                      /*  == original / 10 */
                      /*  == original / 10 */
                      /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ
                       * density. */
                      /*  */
                      /*  */
                      /*  rho = eoscg_densjmd95_bsq(s,t,z) [ kg / m^3 ] */
                      /*  computes the Boussinesq JMD95 in-situ density, given
                       * practical salinity */
                      /*  s, potential temperature t, and depth z [m, positive
                       * and increasing */
                      /*  down].  The depth z is converted into hydrostatic
                       * pressure for a given */
                      /*  gravitational acceleration and Boussinesq reference
                       * density, which are */
                      /*  hard-coded into this function (edit variables grav and
                       * rhob as needed). */
                      /*  */
                      /*  */
                      /*  This function is derived from densjmd95.m, documented
                       * below. Input checks */
                      /*  and expansion of variables have been removed; instead
                       * automatic expansion */
                      /*  (requiring MATLAB 2016b or later) is used. The
                       * calculation has also been */
                      /*  streamlined by pre-allocating arrays for coeffcients,
                       * and modifying the */
                      /*  coefficients such that pressure need not be internally
                       * converted from */
                      /*  dbar to bar. This function is compatible with MATLAB's
                       * codegen. */
                      /*  */
                      /*  The units of s and t are as in densjmd95, documented
                       * below. */
                      /*   */
                      /*  The inputs' sizes are more general: they may be arrays
                       * of any dimension */
                      /*  and any size, so long as their sizes match each other
                       * excepting that any */
                      /*  input can have any dimension as a singleton. For
                       * example, s and t can be */
                      /*  3D arrays of size [nz, nx, ny], while p can be a
                       * vector of size [nz,1]. */
                      /*  */
                      /*  */
                      /*  Author(s)       : Geoff Stanley */
                      /*  Email           : g.stanley@unsw.edu.au */
                      /*  Email           : geoffstanley@gmail.com */
                      /*  Version         : 1.0 */
                      /*  */
                      /*  */
                      /*  DENSJMD95    Density of sea water */
                      /* =========================================================================
                       */
                      /*  */
                      /*  USAGE:  dens = densjmd95(S,Theta,P) */
                      /*  */
                      /*  DESCRIPTION: */
                      /*     Density of Sea Water using Jackett and McDougall
                       * 1995 (JAOT 12) */
                      /*     polynomial (modified UNESCO polynomial). */
                      /*  */
                      /*  INPUT:  (all must have same dimensions) */
                      /*    S     = salinity    [psu      (PSS-78)] */
                      /*    Theta = potential temperature [degree C (IPTS-68)]
                       */
                      /*    P     = pressure    [dbar] */
                      /*        (P may have dims 1x1, mx1, 1xn or mxn for S(mxn)
                       * ) */
                      /*  */
                      /*  OUTPUT: */
                      /*    dens = density  [kg/m^3] */
                      /*  */
                      /*  AUTHOR:  Martin Losch 2002-08-09  (mlosch@mit.edu) */
                      /*  Jackett and McDougall, 1995, JAOT 12(4), pp. 381-388
                       */
                      /*  created by mlosch on 2002-08-09 */
                      /*  $Header:
                       * /u/gcmpack/MITgcm/utils/matlab/densjmd95.m,v 1.2
                       * 2007/02/17 23:49:43 jmc Exp $ */
                      /*  $Name:  $ */
                      /*  */
                      /*  Convert depth [m] to the hydrostatic pressure [dbar]
                       * implied by the */
                      /*  following two hard-coded parameters (edit as needed):
                       */
                      /*  gravitational acceleration [m /s^2] */
                      /*  Boussinesq reference density [kg / m^3] */
                      /*  Pascal to dbar conversion [dbar / Pa] */
                      /*  depth to pressure conversion [dbar / m] */
                      /*  Henceforth z is actually pressure [dbar] */
                      /*  coefficients nonlinear equation of state in pressure
                       * coordinates for */
                      /*  1. density of fresh water at p = 0 */
                      /*  2. density of sea water at p = 0 */
                      /*  coefficients in pressure coordinates for */
                      /*  3. secant bulk modulus K of fresh water at p = 0 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  4. secant bulk modulus K of sea water at p = 0 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  5. secant bulk modulus K of sea water at p */
                      /*  == original / 10 */
                      /*  == original / 10 */
                      /*  == original / 10 */
                      /*  == original / 10 */
                      /*  == original / 10 */
                      /*  == original / 10 */
                      b_s1o2 = muDoubleScalarSqrt(b_s);
                      fbpos =
                          (s1o2 / (1.0 -
                                   z_tmp / (b_d + z_tmp * (e + z_tmp * fc))) -
                               ((b_t * (b_t * (b_t * (b_t * (b_t * 6.536332E-9 +
                                                             -1.120083E-6) +
                                                      0.0001001685) +
                                               -0.00909529) +
                                        0.06793952) +
                                 999.842594) +
                                b_s * (((b_t * (b_t * (b_t * (b_t * 5.3875E-9 +
                                                              -8.2467E-7) +
                                                       7.6438E-5) +
                                                -0.0040899) +
                                         0.824493) +
                                        b_s1o2 * (b_t * (b_t * -1.6546E-6 +
                                                         0.00010227) +
                                                  -0.00572466)) +
                                       b_s * 0.00048314)) /
                                   (1.0 -
                                    z_tmp /
                                        (((b_t *
                                               (b_t *
                                                    (b_t *
                                                         (b_t * -0.0004190253 +
                                                          0.09648704) +
                                                     -17.06103) +
                                                1444.304) +
                                           196593.3) +
                                          b_s *
                                              ((b_t *
                                                    (b_t *
                                                         (b_t * -0.0005084188 +
                                                          0.06283263) +
                                                     -3.101089) +
                                                528.4855) +
                                               b_s1o2 *
                                                   (b_t * (b_t * -0.004619924 +
                                                           0.09085835) +
                                                    3.88664))) +
                                         z_tmp *
                                             (((b_t *
                                                    (b_t * (b_t * 1.956415E-6 +
                                                            -0.0002984642) +
                                                     0.02212276) +
                                                3.186519) +
                                               b_s *
                                                   ((b_t * (b_t * 2.059331E-7 +
                                                            -0.0001847318) +
                                                     0.006704388) +
                                                    b_s1o2 * 0.0001480266)) +
                                              z_tmp *
                                                  ((b_t * (b_t * 1.39468E-8 +
                                                           -1.202016E-6) +
                                                    2.102898E-5) +
                                                   b_s *
                                                       (b_t *
                                                            (b_t *
                                                                 6.207323E-11 +
                                                             6.128773E-9) +
                                                        -2.040237E-7))))) >
                           0.0);
                    }
                  }
                  do {
                    exitg1 = 0;
                    guard1 = false;
                    if (lb > fa) {
                      /*  Move a left, and test for a sign change */
                      dxm *= 1.4142135623730949;
                      lb = muDoubleScalarMax(fb - dxm, fa);
                      /*  Evaluate difference between (a) eos at location on the
                       * cast (S, T, P) */
                      /*  where the pressure or depth is p, and (b) eos of the
                       * bottle (sB, tB, pB); */
                      /*  here, eos is always evaluated at the average pressure
                       * or depth, (p + */
                      /*  pB)/2. */
                      ppc_val2(Pn_data, Pn_size, Sppc_data, Sppc_size,
                               Tppc_data, lb, &b_s, &b_t);
                      /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ
                       * density. */
                      /*  */
                      /*  */
                      /*  rho = eoscg_densjmd95_bsq(s,t,z) [ kg / m^3 ] */
                      /*  computes the Boussinesq JMD95 in-situ density, given
                       * practical salinity */
                      /*  s, potential temperature t, and depth z [m, positive
                       * and increasing */
                      /*  down].  The depth z is converted into hydrostatic
                       * pressure for a given */
                      /*  gravitational acceleration and Boussinesq reference
                       * density, which are */
                      /*  hard-coded into this function (edit variables grav and
                       * rhob as needed). */
                      /*  */
                      /*  */
                      /*  This function is derived from densjmd95.m, documented
                       * below. Input checks */
                      /*  and expansion of variables have been removed; instead
                       * automatic expansion */
                      /*  (requiring MATLAB 2016b or later) is used. The
                       * calculation has also been */
                      /*  streamlined by pre-allocating arrays for coeffcients,
                       * and modifying the */
                      /*  coefficients such that pressure need not be internally
                       * converted from */
                      /*  dbar to bar. This function is compatible with MATLAB's
                       * codegen. */
                      /*  */
                      /*  The units of s and t are as in densjmd95, documented
                       * below. */
                      /*   */
                      /*  The inputs' sizes are more general: they may be arrays
                       * of any dimension */
                      /*  and any size, so long as their sizes match each other
                       * excepting that any */
                      /*  input can have any dimension as a singleton. For
                       * example, s and t can be */
                      /*  3D arrays of size [nz, nx, ny], while p can be a
                       * vector of size [nz,1]. */
                      /*  */
                      /*  */
                      /*  Author(s)       : Geoff Stanley */
                      /*  Email           : g.stanley@unsw.edu.au */
                      /*  Email           : geoffstanley@gmail.com */
                      /*  Version         : 1.0 */
                      /*  */
                      /*  */
                      /*  DENSJMD95    Density of sea water */
                      /* =========================================================================
                       */
                      /*  */
                      /*  USAGE:  dens = densjmd95(S,Theta,P) */
                      /*  */
                      /*  DESCRIPTION: */
                      /*     Density of Sea Water using Jackett and McDougall
                       * 1995 (JAOT 12) */
                      /*     polynomial (modified UNESCO polynomial). */
                      /*  */
                      /*  INPUT:  (all must have same dimensions) */
                      /*    S     = salinity    [psu      (PSS-78)] */
                      /*    Theta = potential temperature [degree C (IPTS-68)]
                       */
                      /*    P     = pressure    [dbar] */
                      /*        (P may have dims 1x1, mx1, 1xn or mxn for S(mxn)
                       * ) */
                      /*  */
                      /*  OUTPUT: */
                      /*    dens = density  [kg/m^3] */
                      /*  */
                      /*  AUTHOR:  Martin Losch 2002-08-09  (mlosch@mit.edu) */
                      /*  Jackett and McDougall, 1995, JAOT 12(4), pp. 381-388
                       */
                      /*  created by mlosch on 2002-08-09 */
                      /*  $Header:
                       * /u/gcmpack/MITgcm/utils/matlab/densjmd95.m,v 1.2
                       * 2007/02/17 23:49:43 jmc Exp $ */
                      /*  $Name:  $ */
                      /*  */
                      /*  Convert depth [m] to the hydrostatic pressure [dbar]
                       * implied by the */
                      /*  following two hard-coded parameters (edit as needed):
                       */
                      /*  gravitational acceleration [m /s^2] */
                      /*  Boussinesq reference density [kg / m^3] */
                      /*  Pascal to dbar conversion [dbar / Pa] */
                      /*  depth to pressure conversion [dbar / m] */
                      z_tmp = (pB + lb) / 2.0 * 1.015335;
                      /*  Henceforth z is actually pressure [dbar] */
                      /*  coefficients nonlinear equation of state in pressure
                       * coordinates for */
                      /*  1. density of fresh water at p = 0 */
                      /*  2. density of sea water at p = 0 */
                      /*  coefficients in pressure coordinates for */
                      /*  3. secant bulk modulus K of fresh water at p = 0 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  4. secant bulk modulus K of sea water at p = 0 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  5. secant bulk modulus K of sea water at p */
                      /*  == original / 10 */
                      /*  == original / 10 */
                      /*  == original / 10 */
                      /*  == original / 10 */
                      /*  == original / 10 */
                      /*  == original / 10 */
                      /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ
                       * density. */
                      /*  */
                      /*  */
                      /*  rho = eoscg_densjmd95_bsq(s,t,z) [ kg / m^3 ] */
                      /*  computes the Boussinesq JMD95 in-situ density, given
                       * practical salinity */
                      /*  s, potential temperature t, and depth z [m, positive
                       * and increasing */
                      /*  down].  The depth z is converted into hydrostatic
                       * pressure for a given */
                      /*  gravitational acceleration and Boussinesq reference
                       * density, which are */
                      /*  hard-coded into this function (edit variables grav and
                       * rhob as needed). */
                      /*  */
                      /*  */
                      /*  This function is derived from densjmd95.m, documented
                       * below. Input checks */
                      /*  and expansion of variables have been removed; instead
                       * automatic expansion */
                      /*  (requiring MATLAB 2016b or later) is used. The
                       * calculation has also been */
                      /*  streamlined by pre-allocating arrays for coeffcients,
                       * and modifying the */
                      /*  coefficients such that pressure need not be internally
                       * converted from */
                      /*  dbar to bar. This function is compatible with MATLAB's
                       * codegen. */
                      /*  */
                      /*  The units of s and t are as in densjmd95, documented
                       * below. */
                      /*   */
                      /*  The inputs' sizes are more general: they may be arrays
                       * of any dimension */
                      /*  and any size, so long as their sizes match each other
                       * excepting that any */
                      /*  input can have any dimension as a singleton. For
                       * example, s and t can be */
                      /*  3D arrays of size [nz, nx, ny], while p can be a
                       * vector of size [nz,1]. */
                      /*  */
                      /*  */
                      /*  Author(s)       : Geoff Stanley */
                      /*  Email           : g.stanley@unsw.edu.au */
                      /*  Email           : geoffstanley@gmail.com */
                      /*  Version         : 1.0 */
                      /*  */
                      /*  */
                      /*  DENSJMD95    Density of sea water */
                      /* =========================================================================
                       */
                      /*  */
                      /*  USAGE:  dens = densjmd95(S,Theta,P) */
                      /*  */
                      /*  DESCRIPTION: */
                      /*     Density of Sea Water using Jackett and McDougall
                       * 1995 (JAOT 12) */
                      /*     polynomial (modified UNESCO polynomial). */
                      /*  */
                      /*  INPUT:  (all must have same dimensions) */
                      /*    S     = salinity    [psu      (PSS-78)] */
                      /*    Theta = potential temperature [degree C (IPTS-68)]
                       */
                      /*    P     = pressure    [dbar] */
                      /*        (P may have dims 1x1, mx1, 1xn or mxn for S(mxn)
                       * ) */
                      /*  */
                      /*  OUTPUT: */
                      /*    dens = density  [kg/m^3] */
                      /*  */
                      /*  AUTHOR:  Martin Losch 2002-08-09  (mlosch@mit.edu) */
                      /*  Jackett and McDougall, 1995, JAOT 12(4), pp. 381-388
                       */
                      /*  created by mlosch on 2002-08-09 */
                      /*  $Header:
                       * /u/gcmpack/MITgcm/utils/matlab/densjmd95.m,v 1.2
                       * 2007/02/17 23:49:43 jmc Exp $ */
                      /*  $Name:  $ */
                      /*  */
                      /*  Convert depth [m] to the hydrostatic pressure [dbar]
                       * implied by the */
                      /*  following two hard-coded parameters (edit as needed):
                       */
                      /*  gravitational acceleration [m /s^2] */
                      /*  Boussinesq reference density [kg / m^3] */
                      /*  Pascal to dbar conversion [dbar / Pa] */
                      /*  depth to pressure conversion [dbar / m] */
                      /*  Henceforth z is actually pressure [dbar] */
                      /*  coefficients nonlinear equation of state in pressure
                       * coordinates for */
                      /*  1. density of fresh water at p = 0 */
                      /*  2. density of sea water at p = 0 */
                      /*  coefficients in pressure coordinates for */
                      /*  3. secant bulk modulus K of fresh water at p = 0 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  4. secant bulk modulus K of sea water at p = 0 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  == original * 10 */
                      /*  5. secant bulk modulus K of sea water at p */
                      /*  == original / 10 */
                      /*  == original / 10 */
                      /*  == original / 10 */
                      /*  == original / 10 */
                      /*  == original / 10 */
                      /*  == original / 10 */
                      b_s1o2 = muDoubleScalarSqrt(b_s);
                      fapos =
                          (s1o2 / (1.0 -
                                   z_tmp / (b_d + z_tmp * (e + z_tmp * fc))) -
                               ((b_t * (b_t * (b_t * (b_t * (b_t * 6.536332E-9 +
                                                             -1.120083E-6) +
                                                      0.0001001685) +
                                               -0.00909529) +
                                        0.06793952) +
                                 999.842594) +
                                b_s * (((b_t * (b_t * (b_t * (b_t * 5.3875E-9 +
                                                              -8.2467E-7) +
                                                       7.6438E-5) +
                                                -0.0040899) +
                                         0.824493) +
                                        b_s1o2 * (b_t * (b_t * -1.6546E-6 +
                                                         0.00010227) +
                                                  -0.00572466)) +
                                       b_s * 0.00048314)) /
                                   (1.0 -
                                    z_tmp /
                                        (((b_t *
                                               (b_t *
                                                    (b_t *
                                                         (b_t * -0.0004190253 +
                                                          0.09648704) +
                                                     -17.06103) +
                                                1444.304) +
                                           196593.3) +
                                          b_s *
                                              ((b_t *
                                                    (b_t *
                                                         (b_t * -0.0005084188 +
                                                          0.06283263) +
                                                     -3.101089) +
                                                528.4855) +
                                               b_s1o2 *
                                                   (b_t * (b_t * -0.004619924 +
                                                           0.09085835) +
                                                    3.88664))) +
                                         z_tmp *
                                             (((b_t *
                                                    (b_t * (b_t * 1.956415E-6 +
                                                            -0.0002984642) +
                                                     0.02212276) +
                                                3.186519) +
                                               b_s *
                                                   ((b_t * (b_t * 2.059331E-7 +
                                                            -0.0001847318) +
                                                     0.006704388) +
                                                    b_s1o2 * 0.0001480266)) +
                                              z_tmp *
                                                  ((b_t * (b_t * 1.39468E-8 +
                                                           -1.202016E-6) +
                                                    2.102898E-5) +
                                                   b_s *
                                                       (b_t *
                                                            (b_t *
                                                                 6.207323E-11 +
                                                             6.128773E-9) +
                                                        -2.040237E-7))))) >
                           0.0);
                      if (fapos != fbpos) {
                        /*  fa and fb have different signs */
                        exitg1 = 1;
                      } else {
                        guard1 = true;
                      }
                    } else if (ub == c) {
                      /*  also a == A, so cannot expand anymore */
                      if (fapos == fbpos) {
                        /*  no sign change found */
                        lb = rtNaN;
                        ub = rtNaN;
                      } else {
                        /*  one last test for sign change */
                      }
                      exitg1 = 1;
                    } else {
                      guard1 = true;
                    }
                    if (guard1) {
                      if (ub < c) {
                        /*  Move b right, and test for a sign change */
                        dxp *= 1.4142135623730949;
                        ub = muDoubleScalarMin(fb + dxp, c);
                        /*  Evaluate difference between (a) eos at location on
                         * the cast (S, T, P) */
                        /*  where the pressure or depth is p, and (b) eos of the
                         * bottle (sB, tB, pB); */
                        /*  here, eos is always evaluated at the average
                         * pressure or depth, (p + */
                        /*  pB)/2. */
                        ppc_val2(Pn_data, Pn_size, Sppc_data, Sppc_size,
                                 Tppc_data, ub, &b_s, &b_t);
                        /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ
                         * density. */
                        /*  */
                        /*  */
                        /*  rho = eoscg_densjmd95_bsq(s,t,z) [ kg / m^3 ] */
                        /*  computes the Boussinesq JMD95 in-situ density, given
                         * practical salinity */
                        /*  s, potential temperature t, and depth z [m, positive
                         * and increasing */
                        /*  down].  The depth z is converted into hydrostatic
                         * pressure for a given */
                        /*  gravitational acceleration and Boussinesq reference
                         * density, which are */
                        /*  hard-coded into this function (edit variables grav
                         * and rhob as needed). */
                        /*  */
                        /*  */
                        /*  This function is derived from densjmd95.m,
                         * documented below. Input checks */
                        /*  and expansion of variables have been removed;
                         * instead automatic expansion */
                        /*  (requiring MATLAB 2016b or later) is used. The
                         * calculation has also been */
                        /*  streamlined by pre-allocating arrays for
                         * coeffcients, and modifying the */
                        /*  coefficients such that pressure need not be
                         * internally converted from */
                        /*  dbar to bar. This function is compatible with
                         * MATLAB's codegen. */
                        /*  */
                        /*  The units of s and t are as in densjmd95, documented
                         * below. */
                        /*   */
                        /*  The inputs' sizes are more general: they may be
                         * arrays of any dimension */
                        /*  and any size, so long as their sizes match each
                         * other excepting that any */
                        /*  input can have any dimension as a singleton. For
                         * example, s and t can be */
                        /*  3D arrays of size [nz, nx, ny], while p can be a
                         * vector of size [nz,1]. */
                        /*  */
                        /*  */
                        /*  Author(s)       : Geoff Stanley */
                        /*  Email           : g.stanley@unsw.edu.au */
                        /*  Email           : geoffstanley@gmail.com */
                        /*  Version         : 1.0 */
                        /*  */
                        /*  */
                        /*  DENSJMD95    Density of sea water */
                        /* =========================================================================
                         */
                        /*  */
                        /*  USAGE:  dens = densjmd95(S,Theta,P) */
                        /*  */
                        /*  DESCRIPTION: */
                        /*     Density of Sea Water using Jackett and McDougall
                         * 1995 (JAOT 12) */
                        /*     polynomial (modified UNESCO polynomial). */
                        /*  */
                        /*  INPUT:  (all must have same dimensions) */
                        /*    S     = salinity    [psu      (PSS-78)] */
                        /*    Theta = potential temperature [degree C (IPTS-68)]
                         */
                        /*    P     = pressure    [dbar] */
                        /*        (P may have dims 1x1, mx1, 1xn or mxn for
                         * S(mxn) ) */
                        /*  */
                        /*  OUTPUT: */
                        /*    dens = density  [kg/m^3] */
                        /*  */
                        /*  AUTHOR:  Martin Losch 2002-08-09  (mlosch@mit.edu)
                         */
                        /*  Jackett and McDougall, 1995, JAOT 12(4), pp. 381-388
                         */
                        /*  created by mlosch on 2002-08-09 */
                        /*  $Header:
                         * /u/gcmpack/MITgcm/utils/matlab/densjmd95.m,v 1.2
                         * 2007/02/17 23:49:43 jmc Exp $ */
                        /*  $Name:  $ */
                        /*  */
                        /*  Convert depth [m] to the hydrostatic pressure [dbar]
                         * implied by the */
                        /*  following two hard-coded parameters (edit as
                         * needed): */
                        /*  gravitational acceleration [m /s^2] */
                        /*  Boussinesq reference density [kg / m^3] */
                        /*  Pascal to dbar conversion [dbar / Pa] */
                        /*  depth to pressure conversion [dbar / m] */
                        z_tmp = (pB + ub) / 2.0 * 1.015335;
                        /*  Henceforth z is actually pressure [dbar] */
                        /*  coefficients nonlinear equation of state in pressure
                         * coordinates for */
                        /*  1. density of fresh water at p = 0 */
                        /*  2. density of sea water at p = 0 */
                        /*  coefficients in pressure coordinates for */
                        /*  3. secant bulk modulus K of fresh water at p = 0 */
                        /*  == original * 10 */
                        /*  == original * 10 */
                        /*  == original * 10 */
                        /*  == original * 10 */
                        /*  == original * 10 */
                        /*  4. secant bulk modulus K of sea water at p = 0 */
                        /*  == original * 10 */
                        /*  == original * 10 */
                        /*  == original * 10 */
                        /*  == original * 10 */
                        /*  == original * 10 */
                        /*  == original * 10 */
                        /*  == original * 10 */
                        /*  5. secant bulk modulus K of sea water at p */
                        /*  == original / 10 */
                        /*  == original / 10 */
                        /*  == original / 10 */
                        /*  == original / 10 */
                        /*  == original / 10 */
                        /*  == original / 10 */
                        /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ
                         * density. */
                        /*  */
                        /*  */
                        /*  rho = eoscg_densjmd95_bsq(s,t,z) [ kg / m^3 ] */
                        /*  computes the Boussinesq JMD95 in-situ density, given
                         * practical salinity */
                        /*  s, potential temperature t, and depth z [m, positive
                         * and increasing */
                        /*  down].  The depth z is converted into hydrostatic
                         * pressure for a given */
                        /*  gravitational acceleration and Boussinesq reference
                         * density, which are */
                        /*  hard-coded into this function (edit variables grav
                         * and rhob as needed). */
                        /*  */
                        /*  */
                        /*  This function is derived from densjmd95.m,
                         * documented below. Input checks */
                        /*  and expansion of variables have been removed;
                         * instead automatic expansion */
                        /*  (requiring MATLAB 2016b or later) is used. The
                         * calculation has also been */
                        /*  streamlined by pre-allocating arrays for
                         * coeffcients, and modifying the */
                        /*  coefficients such that pressure need not be
                         * internally converted from */
                        /*  dbar to bar. This function is compatible with
                         * MATLAB's codegen. */
                        /*  */
                        /*  The units of s and t are as in densjmd95, documented
                         * below. */
                        /*   */
                        /*  The inputs' sizes are more general: they may be
                         * arrays of any dimension */
                        /*  and any size, so long as their sizes match each
                         * other excepting that any */
                        /*  input can have any dimension as a singleton. For
                         * example, s and t can be */
                        /*  3D arrays of size [nz, nx, ny], while p can be a
                         * vector of size [nz,1]. */
                        /*  */
                        /*  */
                        /*  Author(s)       : Geoff Stanley */
                        /*  Email           : g.stanley@unsw.edu.au */
                        /*  Email           : geoffstanley@gmail.com */
                        /*  Version         : 1.0 */
                        /*  */
                        /*  */
                        /*  DENSJMD95    Density of sea water */
                        /* =========================================================================
                         */
                        /*  */
                        /*  USAGE:  dens = densjmd95(S,Theta,P) */
                        /*  */
                        /*  DESCRIPTION: */
                        /*     Density of Sea Water using Jackett and McDougall
                         * 1995 (JAOT 12) */
                        /*     polynomial (modified UNESCO polynomial). */
                        /*  */
                        /*  INPUT:  (all must have same dimensions) */
                        /*    S     = salinity    [psu      (PSS-78)] */
                        /*    Theta = potential temperature [degree C (IPTS-68)]
                         */
                        /*    P     = pressure    [dbar] */
                        /*        (P may have dims 1x1, mx1, 1xn or mxn for
                         * S(mxn) ) */
                        /*  */
                        /*  OUTPUT: */
                        /*    dens = density  [kg/m^3] */
                        /*  */
                        /*  AUTHOR:  Martin Losch 2002-08-09  (mlosch@mit.edu)
                         */
                        /*  Jackett and McDougall, 1995, JAOT 12(4), pp. 381-388
                         */
                        /*  created by mlosch on 2002-08-09 */
                        /*  $Header:
                         * /u/gcmpack/MITgcm/utils/matlab/densjmd95.m,v 1.2
                         * 2007/02/17 23:49:43 jmc Exp $ */
                        /*  $Name:  $ */
                        /*  */
                        /*  Convert depth [m] to the hydrostatic pressure [dbar]
                         * implied by the */
                        /*  following two hard-coded parameters (edit as
                         * needed): */
                        /*  gravitational acceleration [m /s^2] */
                        /*  Boussinesq reference density [kg / m^3] */
                        /*  Pascal to dbar conversion [dbar / Pa] */
                        /*  depth to pressure conversion [dbar / m] */
                        /*  Henceforth z is actually pressure [dbar] */
                        /*  coefficients nonlinear equation of state in pressure
                         * coordinates for */
                        /*  1. density of fresh water at p = 0 */
                        /*  2. density of sea water at p = 0 */
                        /*  coefficients in pressure coordinates for */
                        /*  3. secant bulk modulus K of fresh water at p = 0 */
                        /*  == original * 10 */
                        /*  == original * 10 */
                        /*  == original * 10 */
                        /*  == original * 10 */
                        /*  == original * 10 */
                        /*  4. secant bulk modulus K of sea water at p = 0 */
                        /*  == original * 10 */
                        /*  == original * 10 */
                        /*  == original * 10 */
                        /*  == original * 10 */
                        /*  == original * 10 */
                        /*  == original * 10 */
                        /*  == original * 10 */
                        /*  5. secant bulk modulus K of sea water at p */
                        /*  == original / 10 */
                        /*  == original / 10 */
                        /*  == original / 10 */
                        /*  == original / 10 */
                        /*  == original / 10 */
                        /*  == original / 10 */
                        b_s1o2 = muDoubleScalarSqrt(b_s);
                        fbpos =
                            (s1o2 / (1.0 -
                                     z_tmp / (b_d + z_tmp * (e + z_tmp * fc))) -
                                 ((b_t *
                                       (b_t * (b_t * (b_t * (b_t * 6.536332E-9 +
                                                             -1.120083E-6) +
                                                      0.0001001685) +
                                               -0.00909529) +
                                        0.06793952) +
                                   999.842594) +
                                  b_s *
                                      (((b_t * (b_t * (b_t * (b_t * 5.3875E-9 +
                                                              -8.2467E-7) +
                                                       7.6438E-5) +
                                                -0.0040899) +
                                         0.824493) +
                                        b_s1o2 * (b_t * (b_t * -1.6546E-6 +
                                                         0.00010227) +
                                                  -0.00572466)) +
                                       b_s * 0.00048314)) /
                                     (1.0 -
                                      z_tmp /
                                          (((b_t *
                                                 (b_t *
                                                      (b_t *
                                                           (b_t *
                                                                -0.0004190253 +
                                                            0.09648704) +
                                                       -17.06103) +
                                                  1444.304) +
                                             196593.3) +
                                            b_s *
                                                ((b_t *
                                                      (b_t *
                                                           (b_t *
                                                                -0.0005084188 +
                                                            0.06283263) +
                                                       -3.101089) +
                                                  528.4855) +
                                                 b_s1o2 *
                                                     (b_t *
                                                          (b_t * -0.004619924 +
                                                           0.09085835) +
                                                      3.88664))) +
                                           z_tmp *
                                               (((b_t *
                                                      (b_t *
                                                           (b_t * 1.956415E-6 +
                                                            -0.0002984642) +
                                                       0.02212276) +
                                                  3.186519) +
                                                 b_s *
                                                     ((b_t *
                                                           (b_t * 2.059331E-7 +
                                                            -0.0001847318) +
                                                       0.006704388) +
                                                      b_s1o2 * 0.0001480266)) +
                                                z_tmp *
                                                    ((b_t * (b_t * 1.39468E-8 +
                                                             -1.202016E-6) +
                                                      2.102898E-5) +
                                                     b_s *
                                                         (b_t *
                                                              (b_t *
                                                                   6.207323E-11 +
                                                               6.128773E-9) +
                                                          -2.040237E-7))))) >
                             0.0);
                        if (fapos != fbpos) {
                          /*  fa and fb have different signs */
                          exitg1 = 1;
                        }
                      } else if (lb == fa) {
                        /*  also b == B, so cannot expand anymore */
                        if (fapos == fbpos) {
                          /*  no sign change found */
                          lb = rtNaN;
                          ub = rtNaN;
                        } else {
                          /*  one last test for sign change */
                        }
                        exitg1 = 1;
                      }
                    }
                  } while (exitg1 == 0);
                }
              }
            }
            if (!muDoubleScalarIsNaN(lb)) {
              /*  A sign change was discovered, so a root exists in the
               * interval. */
              /*  Solve the nonlinear root-finding problem using Brent's method
               */
              /* FZERO_BRENT  Find a root of a univariate function within a
               * given interval */
              /*              using Brent's method */
              /*  */
              /*  x = fzero_brent(f,a,b,t)  */
              /*  finds x in the interval [a,b] satisfying |x-y| <= t/2 where
               * f(y) = 0. */
              /*  f(a) and f(b) must have opposite signs. */
              /*  */
              /*  ... = fzero_brent(f,lb,ub,t,x,...) passes additional inputs
               * directly */
              /*  to f. */
              /*  */
              /*  This function is compatible with MATLAB's code generation --
               * so long f is */
              /*  similarly compatible. Many root-finding problems can be solved
               * with */
              /*  fzero_brent by writing another function which calls
               * fzero_brent inside a */
              /*  for loop, and this can be made fast by using code generation
               * on that */
              /*  wrapper function. Note that f will only be called with a
               * scalar input as */
              /*  its first argument; codegen knows this, and might strip out
               * unnecessary */
              /*  code from the function definition underlying f. */
              /*  --- Example: */
              /*  Simple bisection between bounds of -0.5 and 1.5 would fail to
               * find the */
              /*  root. Starting with a guess of .85 and expanding outwards
               * finds a root. */
              /*  This example is shown graphically on the MATLAB File Exchange.
               */
              /*  */
              /*  c = poly([0 1]); % a polynomial with roots at 0 and 1 */
              /*  f = @(x) polyval(c,x); % encapsulate extra parameters */
              /*  [a,b] = fzero_guess_to_bounds(f, .85, -.5, 1.5); */
              /*  root = fzero_brent(f, a, b, .05); */
              /*  */
              /*   Discussion: */
              /*  */
              /*     The interval [A,B] must be a change of sign interval for F.
               */
              /*     That is, F(A) and F(B) must be of opposite signs.  Then */
              /*     assuming that F is continuous implies the existence of at
               * least */
              /*     one value C between A and B for which F(C) = 0. */
              /*  */
              /*     The location of the zero is determined to within an
               * accuracy */
              /*     of 6 * EPS * abs ( C ) + 2 * T, where EPS is the machine
               * epsilon. */
              /*  */
              /*     Thanks to Thomas Secretin for pointing out a transcription
               * error in the */
              /*     setting of the value of P, 11 February 2013. */
              /*  */
              /*     Additional parameters given by varargin are passed directly
               * to F. */
              /*  */
              /*   Licensing: */
              /*  */
              /*     This code is distributed under the GNU LGPL license. */
              /*  */
              /*   Modified: */
              /*  */
              /*     14 September 2019 */
              /*     30 June 2020 - Geoff Stanley */
              /*  */
              /*   Author: */
              /*  */
              /*     Original FORTRAN77 version by Richard Brent. */
              /*     MATLAB version by John Burkardt. */
              /*     Minor changes (passing extra arguments) by Geoff Stanley.
               */
              /*  */
              /*   Reference: */
              /*  */
              /*     Richard Brent, */
              /*     Algorithms for Minimization Without Derivatives, */
              /*     Dover, 2002, */
              /*     ISBN: 0-486-41998-3, */
              /*     LC: QA402.5.B74. */
              /*  */
              /*   Parameters: */
              /*  */
              /*     Input, real A, B, the endpoints of the change of sign
               * interval. */
              /*  */
              /*     Input, real T, a positive error tolerance. */
              /*  */
              /*     Input, real value = F ( x ), the name of a user-supplied */
              /*     function which evaluates the function whose zero is being
               * sought. */
              /*  */
              /*     Output, real VALUE, the estimated value of a zero of */
              /*     the function F. */
              /*  */
              /*  Evaluate difference between (a) eos at location on the cast
               * (S, T, P) */
              /*  where the pressure or depth is p, and (b) eos of the bottle
               * (sB, tB, pB); */
              /*  here, eos is always evaluated at the average pressure or
               * depth, (p + */
              /*  pB)/2. */
              loop_ub = Sppc->size[0];
              nij_idx_0 = Sppc->size[0];
              b_loop_ub = Sppc->size[1];
              Sppc_idx_1 = Sppc->size[1];
              Tppc_idx_0 = Tppc->size[0];
              Tppc_idx_1 = Tppc->size[1];
              b_Sppc_size[0] = Sppc->size[0];
              b_Sppc_size[1] = Sppc->size[1];
              for (i1 = 0; i1 < b_loop_ub; i1++) {
                for (i2 = 0; i2 < loop_ub; i2++) {
                  b_Sppc_data[i2 + b_Sppc_size[0] * i1] =
                      Sppc->data[(i2 + nij_idx_0 * i1) +
                                 nij_idx_0 * Sppc_idx_1 * ((int32_T)n - 1)];
                }
              }
              loop_ub = Tppc->size[0];
              b_loop_ub = Tppc->size[1];
              nij_idx_0 = Tppc->size[0];
              for (i1 = 0; i1 < b_loop_ub; i1++) {
                for (i2 = 0; i2 < loop_ub; i2++) {
                  b_Tppc_data[i2 + nij_idx_0 * i1] =
                      Tppc->data[(i2 + Tppc_idx_0 * i1) +
                                 Tppc_idx_0 * Tppc_idx_1 * ((int32_T)n - 1)];
                }
              }
              ppc_val2(Pn_data, Pn_size, b_Sppc_data, b_Sppc_size, b_Tppc_data,
                       lb, &b_s, &b_t);
              /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
              /*  */
              /*  */
              /*  rho = eoscg_densjmd95_bsq(s,t,z)                             [
               * kg / m^3 ] */
              /*  computes the Boussinesq JMD95 in-situ density, given practical
               * salinity */
              /*  s, potential temperature t, and depth z [m, positive and
               * increasing */
              /*  down].  The depth z is converted into hydrostatic pressure for
               * a given */
              /*  gravitational acceleration and Boussinesq reference density,
               * which are */
              /*  hard-coded into this function (edit variables grav and rhob as
               * needed). */
              /*  */
              /*  */
              /*  This function is derived from densjmd95.m, documented below.
               * Input checks */
              /*  and expansion of variables have been removed; instead
               * automatic expansion */
              /*  (requiring MATLAB 2016b or later) is used. The calculation has
               * also been */
              /*  streamlined by pre-allocating arrays for coeffcients, and
               * modifying the */
              /*  coefficients such that pressure need not be internally
               * converted from */
              /*  dbar to bar. This function is compatible with MATLAB's
               * codegen. */
              /*  */
              /*  The units of s and t are as in densjmd95, documented below. */
              /*   */
              /*  The inputs' sizes are more general: they may be arrays of any
               * dimension */
              /*  and any size, so long as their sizes match each other
               * excepting that any */
              /*  input can have any dimension as a singleton. For example, s
               * and t can be */
              /*  3D arrays of size [nz, nx, ny], while p can be a vector of
               * size [nz,1]. */
              /*  */
              /*  */
              /*  Author(s)       : Geoff Stanley */
              /*  Email           : g.stanley@unsw.edu.au */
              /*  Email           : geoffstanley@gmail.com */
              /*  Version         : 1.0 */
              /*  */
              /*  */
              /*  DENSJMD95    Density of sea water */
              /* =========================================================================
               */
              /*  */
              /*  USAGE:  dens = densjmd95(S,Theta,P) */
              /*  */
              /*  DESCRIPTION: */
              /*     Density of Sea Water using Jackett and McDougall 1995 (JAOT
               * 12) */
              /*     polynomial (modified UNESCO polynomial). */
              /*  */
              /*  INPUT:  (all must have same dimensions) */
              /*    S     = salinity    [psu      (PSS-78)] */
              /*    Theta = potential temperature [degree C (IPTS-68)] */
              /*    P     = pressure    [dbar] */
              /*        (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) ) */
              /*  */
              /*  OUTPUT: */
              /*    dens = density  [kg/m^3] */
              /*  */
              /*  AUTHOR:  Martin Losch 2002-08-09  (mlosch@mit.edu) */
              /*  Jackett and McDougall, 1995, JAOT 12(4), pp. 381-388 */
              /*  created by mlosch on 2002-08-09 */
              /*  $Header: /u/gcmpack/MITgcm/utils/matlab/densjmd95.m,v 1.2
               * 2007/02/17 23:49:43 jmc Exp $ */
              /*  $Name:  $ */
              /*  */
              /*  Convert depth [m] to the hydrostatic pressure [dbar] implied
               * by the */
              /*  following two hard-coded parameters (edit as needed): */
              /*  gravitational acceleration [m /s^2] */
              /*  Boussinesq reference density [kg / m^3] */
              /*  Pascal to dbar conversion [dbar / Pa] */
              /*  depth to pressure conversion [dbar / m] */
              z_tmp = (p->data[(int32_T)m - 1] + lb) / 2.0 * 1.015335;
              /*  Henceforth z is actually pressure [dbar] */
              /*  coefficients nonlinear equation of state in pressure
               * coordinates for */
              /*  1. density of fresh water at p = 0 */
              /*  2. density of sea water at p = 0 */
              /*  coefficients in pressure coordinates for */
              /*  3. secant bulk modulus K of fresh water at p = 0 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  4. secant bulk modulus K of sea water at p = 0 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  5. secant bulk modulus K of sea water at p */
              /*  == original / 10 */
              /*  == original / 10 */
              /*  == original / 10 */
              /*  == original / 10 */
              /*  == original / 10 */
              /*  == original / 10 */
              b_s1o2 = muDoubleScalarSqrt(s->data[(int32_T)m - 1]);
              /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
              /*  */
              /*  */
              /*  rho = eoscg_densjmd95_bsq(s,t,z)                             [
               * kg / m^3 ] */
              /*  computes the Boussinesq JMD95 in-situ density, given practical
               * salinity */
              /*  s, potential temperature t, and depth z [m, positive and
               * increasing */
              /*  down].  The depth z is converted into hydrostatic pressure for
               * a given */
              /*  gravitational acceleration and Boussinesq reference density,
               * which are */
              /*  hard-coded into this function (edit variables grav and rhob as
               * needed). */
              /*  */
              /*  */
              /*  This function is derived from densjmd95.m, documented below.
               * Input checks */
              /*  and expansion of variables have been removed; instead
               * automatic expansion */
              /*  (requiring MATLAB 2016b or later) is used. The calculation has
               * also been */
              /*  streamlined by pre-allocating arrays for coeffcients, and
               * modifying the */
              /*  coefficients such that pressure need not be internally
               * converted from */
              /*  dbar to bar. This function is compatible with MATLAB's
               * codegen. */
              /*  */
              /*  The units of s and t are as in densjmd95, documented below. */
              /*   */
              /*  The inputs' sizes are more general: they may be arrays of any
               * dimension */
              /*  and any size, so long as their sizes match each other
               * excepting that any */
              /*  input can have any dimension as a singleton. For example, s
               * and t can be */
              /*  3D arrays of size [nz, nx, ny], while p can be a vector of
               * size [nz,1]. */
              /*  */
              /*  */
              /*  Author(s)       : Geoff Stanley */
              /*  Email           : g.stanley@unsw.edu.au */
              /*  Email           : geoffstanley@gmail.com */
              /*  Version         : 1.0 */
              /*  */
              /*  */
              /*  DENSJMD95    Density of sea water */
              /* =========================================================================
               */
              /*  */
              /*  USAGE:  dens = densjmd95(S,Theta,P) */
              /*  */
              /*  DESCRIPTION: */
              /*     Density of Sea Water using Jackett and McDougall 1995 (JAOT
               * 12) */
              /*     polynomial (modified UNESCO polynomial). */
              /*  */
              /*  INPUT:  (all must have same dimensions) */
              /*    S     = salinity    [psu      (PSS-78)] */
              /*    Theta = potential temperature [degree C (IPTS-68)] */
              /*    P     = pressure    [dbar] */
              /*        (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) ) */
              /*  */
              /*  OUTPUT: */
              /*    dens = density  [kg/m^3] */
              /*  */
              /*  AUTHOR:  Martin Losch 2002-08-09  (mlosch@mit.edu) */
              /*  Jackett and McDougall, 1995, JAOT 12(4), pp. 381-388 */
              /*  created by mlosch on 2002-08-09 */
              /*  $Header: /u/gcmpack/MITgcm/utils/matlab/densjmd95.m,v 1.2
               * 2007/02/17 23:49:43 jmc Exp $ */
              /*  $Name:  $ */
              /*  */
              /*  Convert depth [m] to the hydrostatic pressure [dbar] implied
               * by the */
              /*  following two hard-coded parameters (edit as needed): */
              /*  gravitational acceleration [m /s^2] */
              /*  Boussinesq reference density [kg / m^3] */
              /*  Pascal to dbar conversion [dbar / Pa] */
              /*  depth to pressure conversion [dbar / m] */
              /*  Henceforth z is actually pressure [dbar] */
              /*  coefficients nonlinear equation of state in pressure
               * coordinates for */
              /*  1. density of fresh water at p = 0 */
              /*  2. density of sea water at p = 0 */
              /*  coefficients in pressure coordinates for */
              /*  3. secant bulk modulus K of fresh water at p = 0 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  4. secant bulk modulus K of sea water at p = 0 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  5. secant bulk modulus K of sea water at p */
              /*  == original / 10 */
              /*  == original / 10 */
              /*  == original / 10 */
              /*  == original / 10 */
              /*  == original / 10 */
              /*  == original / 10 */
              s1o2 = muDoubleScalarSqrt(b_s);
              fa =
                  ((t->data[(int32_T)m - 1] *
                        (t->data[(int32_T)m - 1] *
                             (t->data[(int32_T)m - 1] *
                                  (t->data[(int32_T)m - 1] *
                                       (t->data[(int32_T)m - 1] * 6.536332E-9 +
                                        -1.120083E-6) +
                                   0.0001001685) +
                              -0.00909529) +
                         0.06793952) +
                    999.842594) +
                   s->data[(int32_T)m - 1] *
                       (((t->data[(int32_T)m - 1] *
                              (t->data[(int32_T)m - 1] *
                                   (t->data[(int32_T)m - 1] *
                                        (t->data[(int32_T)m - 1] * 5.3875E-9 +
                                         -8.2467E-7) +
                                    7.6438E-5) +
                               -0.0040899) +
                          0.824493) +
                         b_s1o2 * (t->data[(int32_T)m - 1] *
                                       (t->data[(int32_T)m - 1] * -1.6546E-6 +
                                        0.00010227) +
                                   -0.00572466)) +
                        s->data[(int32_T)m - 1] * 0.00048314)) /
                      (1.0 -
                       z_tmp /
                           (((t->data[(int32_T)m - 1] *
                                  (t->data[(int32_T)m - 1] *
                                       (t->data[(int32_T)m - 1] *
                                            (t->data[(int32_T)m - 1] *
                                                 -0.0004190253 +
                                             0.09648704) +
                                        -17.06103) +
                                   1444.304) +
                              196593.3) +
                             s->data[(int32_T)m - 1] *
                                 ((t->data[(int32_T)m - 1] *
                                       (t->data[(int32_T)m - 1] *
                                            (t->data[(int32_T)m - 1] *
                                                 -0.0005084188 +
                                             0.06283263) +
                                        -3.101089) +
                                   528.4855) +
                                  b_s1o2 * (t->data[(int32_T)m - 1] *
                                                (t->data[(int32_T)m - 1] *
                                                     -0.004619924 +
                                                 0.09085835) +
                                            3.88664))) +
                            z_tmp *
                                (((t->data[(int32_T)m - 1] *
                                       (t->data[(int32_T)m - 1] *
                                            (t->data[(int32_T)m - 1] *
                                                 1.956415E-6 +
                                             -0.0002984642) +
                                        0.02212276) +
                                   3.186519) +
                                  s->data[(int32_T)m - 1] *
                                      ((t->data[(int32_T)m - 1] *
                                            (t->data[(int32_T)m - 1] *
                                                 2.059331E-7 +
                                             -0.0001847318) +
                                        0.006704388) +
                                       b_s1o2 * 0.0001480266)) +
                                 z_tmp * ((t->data[(int32_T)m - 1] *
                                               (t->data[(int32_T)m - 1] *
                                                    1.39468E-8 +
                                                -1.202016E-6) +
                                           2.102898E-5) +
                                          s->data[(int32_T)m - 1] *
                                              (t->data[(int32_T)m - 1] *
                                                   (t->data[(int32_T)m - 1] *
                                                        6.207323E-11 +
                                                    6.128773E-9) +
                                               -2.040237E-7))))) -
                  ((b_t * (b_t * (b_t * (b_t * (b_t * 6.536332E-9 +
                                                -1.120083E-6) +
                                         0.0001001685) +
                                  -0.00909529) +
                           0.06793952) +
                    999.842594) +
                   b_s * (((b_t * (b_t * (b_t * (b_t * 5.3875E-9 + -8.2467E-7) +
                                          7.6438E-5) +
                                   -0.0040899) +
                            0.824493) +
                           s1o2 * (b_t * (b_t * -1.6546E-6 + 0.00010227) +
                                   -0.00572466)) +
                          b_s * 0.00048314)) /
                      (1.0 -
                       z_tmp /
                           (((b_t * (b_t * (b_t * (b_t * -0.0004190253 +
                                                   0.09648704) +
                                            -17.06103) +
                                     1444.304) +
                              196593.3) +
                             b_s * ((b_t * (b_t * (b_t * -0.0005084188 +
                                                   0.06283263) +
                                            -3.101089) +
                                     528.4855) +
                                    s1o2 * (b_t * (b_t * -0.004619924 +
                                                   0.09085835) +
                                            3.88664))) +
                            z_tmp * (((b_t * (b_t * (b_t * 1.956415E-6 +
                                                     -0.0002984642) +
                                              0.02212276) +
                                       3.186519) +
                                      b_s * ((b_t * (b_t * 2.059331E-7 +
                                                     -0.0001847318) +
                                              0.006704388) +
                                             s1o2 * 0.0001480266)) +
                                     z_tmp * ((b_t * (b_t * 1.39468E-8 +
                                                      -1.202016E-6) +
                                               2.102898E-5) +
                                              b_s * (b_t * (b_t * 6.207323E-11 +
                                                            6.128773E-9) +
                                                     -2.040237E-7)))));
              /*  Evaluate difference between (a) eos at location on the cast
               * (S, T, P) */
              /*  where the pressure or depth is p, and (b) eos of the bottle
               * (sB, tB, pB); */
              /*  here, eos is always evaluated at the average pressure or
               * depth, (p + */
              /*  pB)/2. */
              loop_ub = Sppc->size[0];
              nij_idx_0 = Sppc->size[0];
              b_loop_ub = Sppc->size[1];
              Sppc_idx_1 = Sppc->size[1];
              Tppc_idx_0 = Tppc->size[0];
              Tppc_idx_1 = Tppc->size[1];
              b_Sppc_size[0] = Sppc->size[0];
              b_Sppc_size[1] = Sppc->size[1];
              for (i1 = 0; i1 < b_loop_ub; i1++) {
                for (i2 = 0; i2 < loop_ub; i2++) {
                  b_Sppc_data[i2 + b_Sppc_size[0] * i1] =
                      Sppc->data[(i2 + nij_idx_0 * i1) +
                                 nij_idx_0 * Sppc_idx_1 * ((int32_T)n - 1)];
                }
              }
              loop_ub = Tppc->size[0];
              b_loop_ub = Tppc->size[1];
              nij_idx_0 = Tppc->size[0];
              for (i1 = 0; i1 < b_loop_ub; i1++) {
                for (i2 = 0; i2 < loop_ub; i2++) {
                  b_Tppc_data[i2 + nij_idx_0 * i1] =
                      Tppc->data[(i2 + Tppc_idx_0 * i1) +
                                 Tppc_idx_0 * Tppc_idx_1 * ((int32_T)n - 1)];
                }
              }
              ppc_val2(Pn_data, Pn_size, b_Sppc_data, b_Sppc_size, b_Tppc_data,
                       ub, &b_s, &b_t);
              /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
              /*  */
              /*  */
              /*  rho = eoscg_densjmd95_bsq(s,t,z)                             [
               * kg / m^3 ] */
              /*  computes the Boussinesq JMD95 in-situ density, given practical
               * salinity */
              /*  s, potential temperature t, and depth z [m, positive and
               * increasing */
              /*  down].  The depth z is converted into hydrostatic pressure for
               * a given */
              /*  gravitational acceleration and Boussinesq reference density,
               * which are */
              /*  hard-coded into this function (edit variables grav and rhob as
               * needed). */
              /*  */
              /*  */
              /*  This function is derived from densjmd95.m, documented below.
               * Input checks */
              /*  and expansion of variables have been removed; instead
               * automatic expansion */
              /*  (requiring MATLAB 2016b or later) is used. The calculation has
               * also been */
              /*  streamlined by pre-allocating arrays for coeffcients, and
               * modifying the */
              /*  coefficients such that pressure need not be internally
               * converted from */
              /*  dbar to bar. This function is compatible with MATLAB's
               * codegen. */
              /*  */
              /*  The units of s and t are as in densjmd95, documented below. */
              /*   */
              /*  The inputs' sizes are more general: they may be arrays of any
               * dimension */
              /*  and any size, so long as their sizes match each other
               * excepting that any */
              /*  input can have any dimension as a singleton. For example, s
               * and t can be */
              /*  3D arrays of size [nz, nx, ny], while p can be a vector of
               * size [nz,1]. */
              /*  */
              /*  */
              /*  Author(s)       : Geoff Stanley */
              /*  Email           : g.stanley@unsw.edu.au */
              /*  Email           : geoffstanley@gmail.com */
              /*  Version         : 1.0 */
              /*  */
              /*  */
              /*  DENSJMD95    Density of sea water */
              /* =========================================================================
               */
              /*  */
              /*  USAGE:  dens = densjmd95(S,Theta,P) */
              /*  */
              /*  DESCRIPTION: */
              /*     Density of Sea Water using Jackett and McDougall 1995 (JAOT
               * 12) */
              /*     polynomial (modified UNESCO polynomial). */
              /*  */
              /*  INPUT:  (all must have same dimensions) */
              /*    S     = salinity    [psu      (PSS-78)] */
              /*    Theta = potential temperature [degree C (IPTS-68)] */
              /*    P     = pressure    [dbar] */
              /*        (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) ) */
              /*  */
              /*  OUTPUT: */
              /*    dens = density  [kg/m^3] */
              /*  */
              /*  AUTHOR:  Martin Losch 2002-08-09  (mlosch@mit.edu) */
              /*  Jackett and McDougall, 1995, JAOT 12(4), pp. 381-388 */
              /*  created by mlosch on 2002-08-09 */
              /*  $Header: /u/gcmpack/MITgcm/utils/matlab/densjmd95.m,v 1.2
               * 2007/02/17 23:49:43 jmc Exp $ */
              /*  $Name:  $ */
              /*  */
              /*  Convert depth [m] to the hydrostatic pressure [dbar] implied
               * by the */
              /*  following two hard-coded parameters (edit as needed): */
              /*  gravitational acceleration [m /s^2] */
              /*  Boussinesq reference density [kg / m^3] */
              /*  Pascal to dbar conversion [dbar / Pa] */
              /*  depth to pressure conversion [dbar / m] */
              z_tmp = (p->data[(int32_T)m - 1] + ub) / 2.0 * 1.015335;
              /*  Henceforth z is actually pressure [dbar] */
              /*  coefficients nonlinear equation of state in pressure
               * coordinates for */
              /*  1. density of fresh water at p = 0 */
              /*  2. density of sea water at p = 0 */
              /*  coefficients in pressure coordinates for */
              /*  3. secant bulk modulus K of fresh water at p = 0 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  4. secant bulk modulus K of sea water at p = 0 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  5. secant bulk modulus K of sea water at p */
              /*  == original / 10 */
              /*  == original / 10 */
              /*  == original / 10 */
              /*  == original / 10 */
              /*  == original / 10 */
              /*  == original / 10 */
              b_s1o2 = muDoubleScalarSqrt(s->data[(int32_T)m - 1]);
              /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
              /*  */
              /*  */
              /*  rho = eoscg_densjmd95_bsq(s,t,z)                             [
               * kg / m^3 ] */
              /*  computes the Boussinesq JMD95 in-situ density, given practical
               * salinity */
              /*  s, potential temperature t, and depth z [m, positive and
               * increasing */
              /*  down].  The depth z is converted into hydrostatic pressure for
               * a given */
              /*  gravitational acceleration and Boussinesq reference density,
               * which are */
              /*  hard-coded into this function (edit variables grav and rhob as
               * needed). */
              /*  */
              /*  */
              /*  This function is derived from densjmd95.m, documented below.
               * Input checks */
              /*  and expansion of variables have been removed; instead
               * automatic expansion */
              /*  (requiring MATLAB 2016b or later) is used. The calculation has
               * also been */
              /*  streamlined by pre-allocating arrays for coeffcients, and
               * modifying the */
              /*  coefficients such that pressure need not be internally
               * converted from */
              /*  dbar to bar. This function is compatible with MATLAB's
               * codegen. */
              /*  */
              /*  The units of s and t are as in densjmd95, documented below. */
              /*   */
              /*  The inputs' sizes are more general: they may be arrays of any
               * dimension */
              /*  and any size, so long as their sizes match each other
               * excepting that any */
              /*  input can have any dimension as a singleton. For example, s
               * and t can be */
              /*  3D arrays of size [nz, nx, ny], while p can be a vector of
               * size [nz,1]. */
              /*  */
              /*  */
              /*  Author(s)       : Geoff Stanley */
              /*  Email           : g.stanley@unsw.edu.au */
              /*  Email           : geoffstanley@gmail.com */
              /*  Version         : 1.0 */
              /*  */
              /*  */
              /*  DENSJMD95    Density of sea water */
              /* =========================================================================
               */
              /*  */
              /*  USAGE:  dens = densjmd95(S,Theta,P) */
              /*  */
              /*  DESCRIPTION: */
              /*     Density of Sea Water using Jackett and McDougall 1995 (JAOT
               * 12) */
              /*     polynomial (modified UNESCO polynomial). */
              /*  */
              /*  INPUT:  (all must have same dimensions) */
              /*    S     = salinity    [psu      (PSS-78)] */
              /*    Theta = potential temperature [degree C (IPTS-68)] */
              /*    P     = pressure    [dbar] */
              /*        (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) ) */
              /*  */
              /*  OUTPUT: */
              /*    dens = density  [kg/m^3] */
              /*  */
              /*  AUTHOR:  Martin Losch 2002-08-09  (mlosch@mit.edu) */
              /*  Jackett and McDougall, 1995, JAOT 12(4), pp. 381-388 */
              /*  created by mlosch on 2002-08-09 */
              /*  $Header: /u/gcmpack/MITgcm/utils/matlab/densjmd95.m,v 1.2
               * 2007/02/17 23:49:43 jmc Exp $ */
              /*  $Name:  $ */
              /*  */
              /*  Convert depth [m] to the hydrostatic pressure [dbar] implied
               * by the */
              /*  following two hard-coded parameters (edit as needed): */
              /*  gravitational acceleration [m /s^2] */
              /*  Boussinesq reference density [kg / m^3] */
              /*  Pascal to dbar conversion [dbar / Pa] */
              /*  depth to pressure conversion [dbar / m] */
              /*  Henceforth z is actually pressure [dbar] */
              /*  coefficients nonlinear equation of state in pressure
               * coordinates for */
              /*  1. density of fresh water at p = 0 */
              /*  2. density of sea water at p = 0 */
              /*  coefficients in pressure coordinates for */
              /*  3. secant bulk modulus K of fresh water at p = 0 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  4. secant bulk modulus K of sea water at p = 0 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  == original * 10 */
              /*  5. secant bulk modulus K of sea water at p */
              /*  == original / 10 */
              /*  == original / 10 */
              /*  == original / 10 */
              /*  == original / 10 */
              /*  == original / 10 */
              /*  == original / 10 */
              s1o2 = muDoubleScalarSqrt(b_s);
              fb =
                  ((t->data[(int32_T)m - 1] *
                        (t->data[(int32_T)m - 1] *
                             (t->data[(int32_T)m - 1] *
                                  (t->data[(int32_T)m - 1] *
                                       (t->data[(int32_T)m - 1] * 6.536332E-9 +
                                        -1.120083E-6) +
                                   0.0001001685) +
                              -0.00909529) +
                         0.06793952) +
                    999.842594) +
                   s->data[(int32_T)m - 1] *
                       (((t->data[(int32_T)m - 1] *
                              (t->data[(int32_T)m - 1] *
                                   (t->data[(int32_T)m - 1] *
                                        (t->data[(int32_T)m - 1] * 5.3875E-9 +
                                         -8.2467E-7) +
                                    7.6438E-5) +
                               -0.0040899) +
                          0.824493) +
                         b_s1o2 * (t->data[(int32_T)m - 1] *
                                       (t->data[(int32_T)m - 1] * -1.6546E-6 +
                                        0.00010227) +
                                   -0.00572466)) +
                        s->data[(int32_T)m - 1] * 0.00048314)) /
                      (1.0 -
                       z_tmp /
                           (((t->data[(int32_T)m - 1] *
                                  (t->data[(int32_T)m - 1] *
                                       (t->data[(int32_T)m - 1] *
                                            (t->data[(int32_T)m - 1] *
                                                 -0.0004190253 +
                                             0.09648704) +
                                        -17.06103) +
                                   1444.304) +
                              196593.3) +
                             s->data[(int32_T)m - 1] *
                                 ((t->data[(int32_T)m - 1] *
                                       (t->data[(int32_T)m - 1] *
                                            (t->data[(int32_T)m - 1] *
                                                 -0.0005084188 +
                                             0.06283263) +
                                        -3.101089) +
                                   528.4855) +
                                  b_s1o2 * (t->data[(int32_T)m - 1] *
                                                (t->data[(int32_T)m - 1] *
                                                     -0.004619924 +
                                                 0.09085835) +
                                            3.88664))) +
                            z_tmp *
                                (((t->data[(int32_T)m - 1] *
                                       (t->data[(int32_T)m - 1] *
                                            (t->data[(int32_T)m - 1] *
                                                 1.956415E-6 +
                                             -0.0002984642) +
                                        0.02212276) +
                                   3.186519) +
                                  s->data[(int32_T)m - 1] *
                                      ((t->data[(int32_T)m - 1] *
                                            (t->data[(int32_T)m - 1] *
                                                 2.059331E-7 +
                                             -0.0001847318) +
                                        0.006704388) +
                                       b_s1o2 * 0.0001480266)) +
                                 z_tmp * ((t->data[(int32_T)m - 1] *
                                               (t->data[(int32_T)m - 1] *
                                                    1.39468E-8 +
                                                -1.202016E-6) +
                                           2.102898E-5) +
                                          s->data[(int32_T)m - 1] *
                                              (t->data[(int32_T)m - 1] *
                                                   (t->data[(int32_T)m - 1] *
                                                        6.207323E-11 +
                                                    6.128773E-9) +
                                               -2.040237E-7))))) -
                  ((b_t * (b_t * (b_t * (b_t * (b_t * 6.536332E-9 +
                                                -1.120083E-6) +
                                         0.0001001685) +
                                  -0.00909529) +
                           0.06793952) +
                    999.842594) +
                   b_s * (((b_t * (b_t * (b_t * (b_t * 5.3875E-9 + -8.2467E-7) +
                                          7.6438E-5) +
                                   -0.0040899) +
                            0.824493) +
                           s1o2 * (b_t * (b_t * -1.6546E-6 + 0.00010227) +
                                   -0.00572466)) +
                          b_s * 0.00048314)) /
                      (1.0 -
                       z_tmp /
                           (((b_t * (b_t * (b_t * (b_t * -0.0004190253 +
                                                   0.09648704) +
                                            -17.06103) +
                                     1444.304) +
                              196593.3) +
                             b_s * ((b_t * (b_t * (b_t * -0.0005084188 +
                                                   0.06283263) +
                                            -3.101089) +
                                     528.4855) +
                                    s1o2 * (b_t * (b_t * -0.004619924 +
                                                   0.09085835) +
                                            3.88664))) +
                            z_tmp * (((b_t * (b_t * (b_t * 1.956415E-6 +
                                                     -0.0002984642) +
                                              0.02212276) +
                                       3.186519) +
                                      b_s * ((b_t * (b_t * 2.059331E-7 +
                                                     -0.0001847318) +
                                              0.006704388) +
                                             s1o2 * 0.0001480266)) +
                                     z_tmp * ((b_t * (b_t * 1.39468E-8 +
                                                      -1.202016E-6) +
                                               2.102898E-5) +
                                              b_s * (b_t * (b_t * 6.207323E-11 +
                                                            6.128773E-9) +
                                                     -2.040237E-7)))));
              c = lb;
              fc = fa;
              e = ub - lb;
              b_d = e;
              do {
                exitg1 = 0;
                if (muDoubleScalarAbs(fc) < muDoubleScalarAbs(fb)) {
                  lb = ub;
                  ub = c;
                  c = lb;
                  fa = fb;
                  fb = fc;
                  fc = fa;
                }
                dxp = 4.4408920985006262E-16 * muDoubleScalarAbs(ub) + TOL_P;
                dxm = 0.5 * (c - ub);
                if ((muDoubleScalarAbs(dxm) <= dxp) || (fb == 0.0)) {
                  exitg1 = 1;
                } else {
                  if ((muDoubleScalarAbs(e) < dxp) ||
                      (muDoubleScalarAbs(fa) <= muDoubleScalarAbs(fb))) {
                    e = dxm;
                    b_d = dxm;
                  } else {
                    b_s = fb / fa;
                    if (lb == c) {
                      s1o2 = 2.0 * dxm * b_s;
                      fa = 1.0 - b_s;
                    } else {
                      fa /= fc;
                      b_r = fb / fc;
                      s1o2 = b_s * (2.0 * dxm * fa * (fa - b_r) -
                                    (ub - lb) * (b_r - 1.0));
                      fa = (fa - 1.0) * (b_r - 1.0) * (b_s - 1.0);
                    }
                    if (0.0 < s1o2) {
                      fa = -fa;
                    } else {
                      s1o2 = -s1o2;
                    }
                    b_s = e;
                    e = b_d;
                    if ((2.0 * s1o2 <
                         3.0 * dxm * fa - muDoubleScalarAbs(dxp * fa)) &&
                        (s1o2 < muDoubleScalarAbs(0.5 * b_s * fa))) {
                      b_d = s1o2 / fa;
                    } else {
                      e = dxm;
                      b_d = dxm;
                    }
                  }
                  lb = ub;
                  fa = fb;
                  if (dxp < muDoubleScalarAbs(b_d)) {
                    ub += b_d;
                  } else if (0.0 < dxm) {
                    ub += dxp;
                  } else {
                    ub -= dxp;
                  }
                  /*  Evaluate difference between (a) eos at location on the
                   * cast (S, T, P) */
                  /*  where the pressure or depth is p, and (b) eos of the
                   * bottle (sB, tB, pB); */
                  /*  here, eos is always evaluated at the average pressure or
                   * depth, (p + */
                  /*  pB)/2. */
                  ppc_val2(Pn_data, Pn_size, Sppc_data, Sppc_size, Tppc_data,
                           ub, &b_s, &b_t);
                  /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density.
                   */
                  /*  */
                  /*  */
                  /*  rho = eoscg_densjmd95_bsq(s,t,z) [ kg / m^3 ] */
                  /*  computes the Boussinesq JMD95 in-situ density, given
                   * practical salinity */
                  /*  s, potential temperature t, and depth z [m, positive and
                   * increasing */
                  /*  down].  The depth z is converted into hydrostatic pressure
                   * for a given */
                  /*  gravitational acceleration and Boussinesq reference
                   * density, which are */
                  /*  hard-coded into this function (edit variables grav and
                   * rhob as needed). */
                  /*  */
                  /*  */
                  /*  This function is derived from densjmd95.m, documented
                   * below. Input checks */
                  /*  and expansion of variables have been removed; instead
                   * automatic expansion */
                  /*  (requiring MATLAB 2016b or later) is used. The calculation
                   * has also been */
                  /*  streamlined by pre-allocating arrays for coeffcients, and
                   * modifying the */
                  /*  coefficients such that pressure need not be internally
                   * converted from */
                  /*  dbar to bar. This function is compatible with MATLAB's
                   * codegen. */
                  /*  */
                  /*  The units of s and t are as in densjmd95, documented
                   * below. */
                  /*   */
                  /*  The inputs' sizes are more general: they may be arrays of
                   * any dimension */
                  /*  and any size, so long as their sizes match each other
                   * excepting that any */
                  /*  input can have any dimension as a singleton. For example,
                   * s and t can be */
                  /*  3D arrays of size [nz, nx, ny], while p can be a vector of
                   * size [nz,1]. */
                  /*  */
                  /*  */
                  /*  Author(s)       : Geoff Stanley */
                  /*  Email           : g.stanley@unsw.edu.au */
                  /*  Email           : geoffstanley@gmail.com */
                  /*  Version         : 1.0 */
                  /*  */
                  /*  */
                  /*  DENSJMD95    Density of sea water */
                  /* =========================================================================
                   */
                  /*  */
                  /*  USAGE:  dens = densjmd95(S,Theta,P) */
                  /*  */
                  /*  DESCRIPTION: */
                  /*     Density of Sea Water using Jackett and McDougall 1995
                   * (JAOT 12) */
                  /*     polynomial (modified UNESCO polynomial). */
                  /*  */
                  /*  INPUT:  (all must have same dimensions) */
                  /*    S     = salinity    [psu      (PSS-78)] */
                  /*    Theta = potential temperature [degree C (IPTS-68)] */
                  /*    P     = pressure    [dbar] */
                  /*        (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) )
                   */
                  /*  */
                  /*  OUTPUT: */
                  /*    dens = density  [kg/m^3] */
                  /*  */
                  /*  AUTHOR:  Martin Losch 2002-08-09  (mlosch@mit.edu) */
                  /*  Jackett and McDougall, 1995, JAOT 12(4), pp. 381-388 */
                  /*  created by mlosch on 2002-08-09 */
                  /*  $Header: /u/gcmpack/MITgcm/utils/matlab/densjmd95.m,v 1.2
                   * 2007/02/17 23:49:43 jmc Exp $ */
                  /*  $Name:  $ */
                  /*  */
                  /*  Convert depth [m] to the hydrostatic pressure [dbar]
                   * implied by the */
                  /*  following two hard-coded parameters (edit as needed): */
                  /*  gravitational acceleration [m /s^2] */
                  /*  Boussinesq reference density [kg / m^3] */
                  /*  Pascal to dbar conversion [dbar / Pa] */
                  /*  depth to pressure conversion [dbar / m] */
                  z_tmp = (pB + ub) / 2.0 * 1.015335;
                  /*  Henceforth z is actually pressure [dbar] */
                  /*  coefficients nonlinear equation of state in pressure
                   * coordinates for */
                  /*  1. density of fresh water at p = 0 */
                  /*  2. density of sea water at p = 0 */
                  /*  coefficients in pressure coordinates for */
                  /*  3. secant bulk modulus K of fresh water at p = 0 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  4. secant bulk modulus K of sea water at p = 0 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  5. secant bulk modulus K of sea water at p */
                  /*  == original / 10 */
                  /*  == original / 10 */
                  /*  == original / 10 */
                  /*  == original / 10 */
                  /*  == original / 10 */
                  /*  == original / 10 */
                  b_s1o2 = muDoubleScalarSqrt(sB);
                  /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density.
                   */
                  /*  */
                  /*  */
                  /*  rho = eoscg_densjmd95_bsq(s,t,z) [ kg / m^3 ] */
                  /*  computes the Boussinesq JMD95 in-situ density, given
                   * practical salinity */
                  /*  s, potential temperature t, and depth z [m, positive and
                   * increasing */
                  /*  down].  The depth z is converted into hydrostatic pressure
                   * for a given */
                  /*  gravitational acceleration and Boussinesq reference
                   * density, which are */
                  /*  hard-coded into this function (edit variables grav and
                   * rhob as needed). */
                  /*  */
                  /*  */
                  /*  This function is derived from densjmd95.m, documented
                   * below. Input checks */
                  /*  and expansion of variables have been removed; instead
                   * automatic expansion */
                  /*  (requiring MATLAB 2016b or later) is used. The calculation
                   * has also been */
                  /*  streamlined by pre-allocating arrays for coeffcients, and
                   * modifying the */
                  /*  coefficients such that pressure need not be internally
                   * converted from */
                  /*  dbar to bar. This function is compatible with MATLAB's
                   * codegen. */
                  /*  */
                  /*  The units of s and t are as in densjmd95, documented
                   * below. */
                  /*   */
                  /*  The inputs' sizes are more general: they may be arrays of
                   * any dimension */
                  /*  and any size, so long as their sizes match each other
                   * excepting that any */
                  /*  input can have any dimension as a singleton. For example,
                   * s and t can be */
                  /*  3D arrays of size [nz, nx, ny], while p can be a vector of
                   * size [nz,1]. */
                  /*  */
                  /*  */
                  /*  Author(s)       : Geoff Stanley */
                  /*  Email           : g.stanley@unsw.edu.au */
                  /*  Email           : geoffstanley@gmail.com */
                  /*  Version         : 1.0 */
                  /*  */
                  /*  */
                  /*  DENSJMD95    Density of sea water */
                  /* =========================================================================
                   */
                  /*  */
                  /*  USAGE:  dens = densjmd95(S,Theta,P) */
                  /*  */
                  /*  DESCRIPTION: */
                  /*     Density of Sea Water using Jackett and McDougall 1995
                   * (JAOT 12) */
                  /*     polynomial (modified UNESCO polynomial). */
                  /*  */
                  /*  INPUT:  (all must have same dimensions) */
                  /*    S     = salinity    [psu      (PSS-78)] */
                  /*    Theta = potential temperature [degree C (IPTS-68)] */
                  /*    P     = pressure    [dbar] */
                  /*        (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) )
                   */
                  /*  */
                  /*  OUTPUT: */
                  /*    dens = density  [kg/m^3] */
                  /*  */
                  /*  AUTHOR:  Martin Losch 2002-08-09  (mlosch@mit.edu) */
                  /*  Jackett and McDougall, 1995, JAOT 12(4), pp. 381-388 */
                  /*  created by mlosch on 2002-08-09 */
                  /*  $Header: /u/gcmpack/MITgcm/utils/matlab/densjmd95.m,v 1.2
                   * 2007/02/17 23:49:43 jmc Exp $ */
                  /*  $Name:  $ */
                  /*  */
                  /*  Convert depth [m] to the hydrostatic pressure [dbar]
                   * implied by the */
                  /*  following two hard-coded parameters (edit as needed): */
                  /*  gravitational acceleration [m /s^2] */
                  /*  Boussinesq reference density [kg / m^3] */
                  /*  Pascal to dbar conversion [dbar / Pa] */
                  /*  depth to pressure conversion [dbar / m] */
                  /*  Henceforth z is actually pressure [dbar] */
                  /*  coefficients nonlinear equation of state in pressure
                   * coordinates for */
                  /*  1. density of fresh water at p = 0 */
                  /*  2. density of sea water at p = 0 */
                  /*  coefficients in pressure coordinates for */
                  /*  3. secant bulk modulus K of fresh water at p = 0 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  4. secant bulk modulus K of sea water at p = 0 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  == original * 10 */
                  /*  5. secant bulk modulus K of sea water at p */
                  /*  == original / 10 */
                  /*  == original / 10 */
                  /*  == original / 10 */
                  /*  == original / 10 */
                  /*  == original / 10 */
                  /*  == original / 10 */
                  s1o2 = muDoubleScalarSqrt(b_s);
                  fb =
                      ((tB * (tB * (tB * (tB * (tB * 6.536332E-9 +
                                                -1.120083E-6) +
                                          0.0001001685) +
                                    -0.00909529) +
                              0.06793952) +
                        999.842594) +
                       sB * (((tB * (tB * (tB * (tB * 5.3875E-9 + -8.2467E-7) +
                                           7.6438E-5) +
                                     -0.0040899) +
                               0.824493) +
                              b_s1o2 * (tB * (tB * -1.6546E-6 + 0.00010227) +
                                        -0.00572466)) +
                             sB * 0.00048314)) /
                          (1.0 -
                           z_tmp /
                               (((tB * (tB * (tB * (tB * -0.0004190253 +
                                                    0.09648704) +
                                              -17.06103) +
                                        1444.304) +
                                  196593.3) +
                                 sB * ((tB * (tB * (tB * -0.0005084188 +
                                                    0.06283263) +
                                              -3.101089) +
                                        528.4855) +
                                       b_s1o2 * (tB * (tB * -0.004619924 +
                                                       0.09085835) +
                                                 3.88664))) +
                                z_tmp *
                                    (((tB * (tB * (tB * 1.956415E-6 +
                                                   -0.0002984642) +
                                             0.02212276) +
                                       3.186519) +
                                      sB * ((tB * (tB * 2.059331E-7 +
                                                   -0.0001847318) +
                                             0.006704388) +
                                            b_s1o2 * 0.0001480266)) +
                                     z_tmp * ((tB * (tB * 1.39468E-8 +
                                                     -1.202016E-6) +
                                               2.102898E-5) +
                                              sB * (tB * (tB * 6.207323E-11 +
                                                          6.128773E-9) +
                                                    -2.040237E-7))))) -
                      ((b_t * (b_t * (b_t * (b_t * (b_t * 6.536332E-9 +
                                                    -1.120083E-6) +
                                             0.0001001685) +
                                      -0.00909529) +
                               0.06793952) +
                        999.842594) +
                       b_s * (((b_t * (b_t * (b_t * (b_t * 5.3875E-9 +
                                                     -8.2467E-7) +
                                              7.6438E-5) +
                                       -0.0040899) +
                                0.824493) +
                               s1o2 * (b_t * (b_t * -1.6546E-6 + 0.00010227) +
                                       -0.00572466)) +
                              b_s * 0.00048314)) /
                          (1.0 -
                           z_tmp /
                               (((b_t * (b_t * (b_t * (b_t * -0.0004190253 +
                                                       0.09648704) +
                                                -17.06103) +
                                         1444.304) +
                                  196593.3) +
                                 b_s * ((b_t * (b_t * (b_t * -0.0005084188 +
                                                       0.06283263) +
                                                -3.101089) +
                                         528.4855) +
                                        s1o2 * (b_t * (b_t * -0.004619924 +
                                                       0.09085835) +
                                                3.88664))) +
                                z_tmp *
                                    (((b_t * (b_t * (b_t * 1.956415E-6 +
                                                     -0.0002984642) +
                                              0.02212276) +
                                       3.186519) +
                                      b_s * ((b_t * (b_t * 2.059331E-7 +
                                                     -0.0001847318) +
                                              0.006704388) +
                                             s1o2 * 0.0001480266)) +
                                     z_tmp * ((b_t * (b_t * 1.39468E-8 +
                                                      -1.202016E-6) +
                                               2.102898E-5) +
                                              b_s * (b_t * (b_t * 6.207323E-11 +
                                                            6.128773E-9) +
                                                     -2.040237E-7)))));
                  if (((0.0 < fb) && (0.0 < fc)) ||
                      ((fb <= 0.0) && (fc <= 0.0))) {
                    c = lb;
                    fc = fa;
                    e = ub - lb;
                    b_d = e;
                  }
                }
              } while (exitg1 == 0);
              /*  Interpolate S and T onto the updated surface */
            } else {
              ub = rtNaN;
            }
          } else {
            ub = rtNaN;
          }
          p->data[(int32_T)n - 1] = ub;
          if ((!muDoubleScalarIsInf(p->data[(int32_T)n - 1])) &&
              (!muDoubleScalarIsNaN(p->data[(int32_T)n - 1])) &&
              (p->data[(int32_T)n - 1] > ML->data[(int32_T)n - 1])) {
            /*  The NTP connection was successful, and its location on the */
            /*  neighbouring cast is below the mixed layer. */
            loop_ub = Sppc->size[0];
            nij_idx_0 = Sppc->size[0];
            b_loop_ub = Sppc->size[1];
            Sppc_idx_1 = Sppc->size[1];
            Tppc_idx_0 = Tppc->size[0];
            Tppc_idx_1 = Tppc->size[1];
            Sppc_size[0] = Sppc->size[0];
            Sppc_size[1] = Sppc->size[1];
            for (i1 = 0; i1 < b_loop_ub; i1++) {
              for (i2 = 0; i2 < loop_ub; i2++) {
                Sppc_data[i2 + Sppc_size[0] * i1] =
                    Sppc->data[(i2 + nij_idx_0 * i1) +
                               nij_idx_0 * Sppc_idx_1 * ((int32_T)n - 1)];
              }
            }
            loop_ub = Tppc->size[0];
            b_loop_ub = Tppc->size[1];
            nij_idx_0 = Tppc->size[0];
            for (i1 = 0; i1 < b_loop_ub; i1++) {
              for (i2 = 0; i2 < loop_ub; i2++) {
                Tppc_data[i2 + nij_idx_0 * i1] =
                    Tppc->data[(i2 + Tppc_idx_0 * i1) +
                               Tppc_idx_0 * Tppc_idx_1 * ((int32_T)n - 1)];
              }
            }
            ppc_val2(Pn_data, Pn_size, Sppc_data, Sppc_size, Tppc_data,
                     p->data[(int32_T)n - 1], &s->data[(int32_T)n - 1],
                     &t->data[(int32_T)n - 1]);
            (*qt)++;
            /*  Add n to queue */
            qu->data[(int32_T)*qt - 1] = n;
            G->data[(int32_T)n - 1] = false;
            /*  mark n as discovered */
            dry->data[(int32_T)n - 1] = false;
            (*freshly_wet)++;
            /*  augment counter of freshly wet casts */
          }
        }
      }
    }
  }
  emxFree_boolean_T(&dry);
  emxFree_boolean_T(&G);
  emlrtHeapReferenceStackLeaveFcnR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (bfs_conncomp1_wet.c) */
