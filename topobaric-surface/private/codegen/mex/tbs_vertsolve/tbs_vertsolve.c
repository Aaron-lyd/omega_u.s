/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * tbs_vertsolve.c
 *
 * Code generation for function 'tbs_vertsolve'
 *
 */

/* Include files */
#include "tbs_vertsolve.h"
#include "ppc_val2.h"
#include "rt_nonfinite.h"
#include "tbs_vertsolve_data.h"
#include "tbs_vertsolve_emxutil.h"
#include "tbs_vertsolve_types.h"
#include "mwmathutil.h"

/* Function Declarations */
static real_T myfcn(real_T p, const real_T Sppc_data[],
                    const int32_T Sppc_size[2], const real_T Tppc_data[],
                    const emxArray_real_T *P, const real_T d_fn_data[],
                    int32_T d_fn_size, real_T s0, real_T t0);

/* Function Definitions */
static real_T myfcn(real_T p, const real_T Sppc_data[],
                    const int32_T Sppc_size[2], const real_T Tppc_data[],
                    const emxArray_real_T *P, const real_T d_fn_data[],
                    int32_T d_fn_size, real_T s0, real_T t0)
{
  real_T distL;
  real_T distR;
  real_T s;
  real_T s1o2;
  real_T slope;
  real_T t;
  real_T xp;
  real_T y;
  int32_T k;
  /*  The difference between delta evaluated (a) using the local branch of the
   */
  /*  multivalued function, and (b) using the equation of state with the local
   */
  /*  water properties. */
  /*  Evaluate water properties on the surface */
  ppc_val2(P, Sppc_data, Sppc_size, Tppc_data, p, &s, &t);
  /*  Evaluate the delta difference */
  /* PVALLIN  Evaluate polynomial, linearly extrapolating outside domain. */
  /*  */
  /*  */
  /*  y = pvallin(f,x) */
  /*  evaluates the polynomial f at location site x if x is in the domain of f,
   */
  /*  or linearly extrapolates f to x if x is outside the domain of f. The */
  /*  domain of f is f(1) to f(2). If f(1) <= x <= f(2), then */
  /*  y = f(3)*(x-f(1))^(N-1) + ... + f(N+1)*(x-f(1)) + f(N+2,l). */
  /*  If x < f(1), y is linearly extrapolated from the left end of f. */
  /*  If f(2) < x , y is linearly extrapolated from the right end of f. */
  /*  */
  /*  */
  /*  --- Input: */
  /*  f [N+2, 1]: polynomial of order N, the valid in a domain f(1) to f(2) and
   */
  /*              with coefficients f(3:N+2) starting with the highest order */
  /*              monomial */
  /*  x [1  , 1]: evaluation site */
  /*  */
  /*  */
  /*  --- Output: */
  /*  y [1  , 1]: f evaluated at x */
  /*  Author(s) : Geoff Stanley */
  /*  Email     : g.stanley@unsw.edu.au */
  /*  Email     : geoffstanley@gmail.com */
  /*  Np2-2 is the polynomial's order. (e.g. Np2==4 for linear) */
  /*  Evaluate local coordinate: */
  distL = p - d_fn_data[0];
  if (distL < 0.0) {
    /*  Interpolate linearly from left. */
    /*    val@left + slope@left * (dist from left) */
    y = d_fn_data[d_fn_size - 1] + d_fn_data[d_fn_size - 2] * distL;
  } else {
    xp = 1.0;
    y = d_fn_data[d_fn_size - 1];
    distR = p - d_fn_data[1];
    if (distR > 0.0) {
      /*  Interpolate linearly from right. */
      /*  Change to evaluate at right, and change to local coordinates */
      distL = d_fn_data[1] - d_fn_data[0];
      /*  Evaluate polynomial AND slope at right: */
      slope = 0.0;
      for (k = 0; k <= d_fn_size - 4; k++) {
        s1o2 = d_fn_data[(d_fn_size - k) - 2];
        slope += ((real_T)k + 1.0) * s1o2 * xp;
        /*  xp == (f(2)-f(1)) ^ (k-1) */
        xp *= distL;
        y += s1o2 * xp;
        /*  xp == (f(2)-f(1)) ^ k */
      }
      /*    val@right + slope@right * (dist from right) */
      y += slope * distR;
    } else {
      /*  Evaluate polynomial: */
      for (k = 0; k <= d_fn_size - 4; k++) {
        xp *= distL;
        /*  xp == (x-f(1)) ^ k */
        y += d_fn_data[(d_fn_size - k) - 2] * xp;
      }
    }
  }
  /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
  /*  */
  /*  */
  /*  rho = eoscg_densjmd95_bsq(s,t,z)                             [ kg / m^3 ]
   */
  /*  computes the Boussinesq JMD95 in-situ density, given practical salinity */
  /*  s, potential temperature t, and depth z [m, positive and increasing */
  /*  down].  The depth z is converted into hydrostatic pressure for a given */
  /*  gravitational acceleration and Boussinesq reference density, which are */
  /*  hard-coded into this function (edit variables grav and rhob as needed). */
  /*  */
  /*  */
  /*  This function is derived from densjmd95.m, documented below. Input checks
   */
  /*  and expansion of variables have been removed; instead automatic expansion
   */
  /*  (requiring MATLAB 2016b or later) is used. The calculation has also been
   */
  /*  streamlined by pre-allocating arrays for coeffcients, and modifying the */
  /*  coefficients such that pressure need not be internally converted from */
  /*  dbar to bar. This function is compatible with MATLAB's codegen. */
  /*  */
  /*  The units of s and t are as in densjmd95, documented below. */
  /*   */
  /*  The inputs' sizes are more general: they may be arrays of any dimension */
  /*  and any size, so long as their sizes match each other excepting that any
   */
  /*  input can have any dimension as a singleton. For example, s and t can be
   */
  /*  3D arrays of size [nz, nx, ny], while p can be a vector of size [nz,1]. */
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
  /*     Density of Sea Water using Jackett and McDougall 1995 (JAOT 12) */
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
  /*  $Header: /u/gcmpack/MITgcm/utils/matlab/densjmd95.m,v 1.2 2007/02/17
   * 23:49:43 jmc Exp $ */
  /*  $Name:  $ */
  /*  */
  /*  Convert depth [m] to the hydrostatic pressure [dbar] implied by the */
  /*  following two hard-coded parameters (edit as needed): */
  /*  gravitational acceleration [m /s^2] */
  /*  Boussinesq reference density [kg / m^3] */
  /*  Pascal to dbar conversion [dbar / Pa] */
  /*  depth to pressure conversion [dbar / m] */
  distL = p * 1.015335;
  /*  Henceforth z is actually pressure [dbar] */
  /*  coefficients nonlinear equation of state in pressure coordinates for */
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
  slope = muDoubleScalarSqrt(s);
  /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
  /*  */
  /*  */
  /*  rho = eoscg_densjmd95_bsq(s,t,z)                             [ kg / m^3 ]
   */
  /*  computes the Boussinesq JMD95 in-situ density, given practical salinity */
  /*  s, potential temperature t, and depth z [m, positive and increasing */
  /*  down].  The depth z is converted into hydrostatic pressure for a given */
  /*  gravitational acceleration and Boussinesq reference density, which are */
  /*  hard-coded into this function (edit variables grav and rhob as needed). */
  /*  */
  /*  */
  /*  This function is derived from densjmd95.m, documented below. Input checks
   */
  /*  and expansion of variables have been removed; instead automatic expansion
   */
  /*  (requiring MATLAB 2016b or later) is used. The calculation has also been
   */
  /*  streamlined by pre-allocating arrays for coeffcients, and modifying the */
  /*  coefficients such that pressure need not be internally converted from */
  /*  dbar to bar. This function is compatible with MATLAB's codegen. */
  /*  */
  /*  The units of s and t are as in densjmd95, documented below. */
  /*   */
  /*  The inputs' sizes are more general: they may be arrays of any dimension */
  /*  and any size, so long as their sizes match each other excepting that any
   */
  /*  input can have any dimension as a singleton. For example, s and t can be
   */
  /*  3D arrays of size [nz, nx, ny], while p can be a vector of size [nz,1]. */
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
  /*     Density of Sea Water using Jackett and McDougall 1995 (JAOT 12) */
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
  /*  $Header: /u/gcmpack/MITgcm/utils/matlab/densjmd95.m,v 1.2 2007/02/17
   * 23:49:43 jmc Exp $ */
  /*  $Name:  $ */
  /*  */
  /*  Convert depth [m] to the hydrostatic pressure [dbar] implied by the */
  /*  following two hard-coded parameters (edit as needed): */
  /*  gravitational acceleration [m /s^2] */
  /*  Boussinesq reference density [kg / m^3] */
  /*  Pascal to dbar conversion [dbar / Pa] */
  /*  depth to pressure conversion [dbar / m] */
  /*  Henceforth z is actually pressure [dbar] */
  /*  coefficients nonlinear equation of state in pressure coordinates for */
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
  s1o2 = muDoubleScalarSqrt(s0);
  return y -
         (((t * (t * (t * (t * (t * 6.536332E-9 + -1.120083E-6) +
                           0.0001001685) +
                      -0.00909529) +
                 0.06793952) +
            999.842594) +
           s * (((t * (t * (t * (t * 5.3875E-9 + -8.2467E-7) + 7.6438E-5) +
                       -0.0040899) +
                  0.824493) +
                 slope * (t * (t * -1.6546E-6 + 0.00010227) + -0.00572466)) +
                s * 0.00048314)) /
              (1.0 -
               distL / (((t * (t * (t * (t * -0.0004190253 + 0.09648704) +
                                    -17.06103) +
                               1444.304) +
                          196593.3) +
                         s * ((t * (t * (t * -0.0005084188 + 0.06283263) +
                                    -3.101089) +
                               528.4855) +
                              slope * (t * (t * -0.004619924 + 0.09085835) +
                                       3.88664))) +
                        distL * (((t * (t * (t * 1.956415E-6 + -0.0002984642) +
                                        0.02212276) +
                                   3.186519) +
                                  s * ((t * (t * 2.059331E-7 + -0.0001847318) +
                                        0.006704388) +
                                       slope * 0.0001480266)) +
                                 distL * ((t * (t * 1.39468E-8 + -1.202016E-6) +
                                           2.102898E-5) +
                                          s * (t * (t * 6.207323E-11 +
                                                    6.128773E-9) +
                                               -2.040237E-7))))) -
          ((t0 * (t0 * (t0 * (t0 * (t0 * 6.536332E-9 + -1.120083E-6) +
                              0.0001001685) +
                        -0.00909529) +
                  0.06793952) +
            999.842594) +
           s0 * (((t0 * (t0 * (t0 * (t0 * 5.3875E-9 + -8.2467E-7) + 7.6438E-5) +
                         -0.0040899) +
                   0.824493) +
                  s1o2 * (t0 * (t0 * -1.6546E-6 + 0.00010227) + -0.00572466)) +
                 s0 * 0.00048314)) /
              (1.0 -
               distL /
                   (((t0 * (t0 * (t0 * (t0 * -0.0004190253 + 0.09648704) +
                                  -17.06103) +
                            1444.304) +
                      196593.3) +
                     s0 * ((t0 * (t0 * (t0 * -0.0005084188 + 0.06283263) +
                                  -3.101089) +
                            528.4855) +
                           s1o2 * (t0 * (t0 * -0.004619924 + 0.09085835) +
                                   3.88664))) +
                    distL *
                        (((t0 * (t0 * (t0 * 1.956415E-6 + -0.0002984642) +
                                 0.02212276) +
                           3.186519) +
                          s0 * ((t0 * (t0 * 2.059331E-7 + -0.0001847318) +
                                 0.006704388) +
                                s1o2 * 0.0001480266)) +
                         distL * ((t0 * (t0 * 1.39468E-8 + -1.202016E-6) +
                                   2.102898E-5) +
                                  s0 * (t0 * (t0 * 6.207323E-11 + 6.128773E-9) +
                                        -2.040237E-7))))));
}

void tbs_vertsolve(const emxArray_real_T *Sppc, const emxArray_real_T *Tppc,
                   const emxArray_real_T *P, const emxArray_real_T *BotK,
                   emxArray_real_T *s, emxArray_real_T *t, emxArray_real_T *p,
                   const emxArray_real_T *branchmap,
                   const emxArray_real_T *d_fn, real_T s0, real_T t0,
                   real_T tolp, real_T DP)
{
  emxArray_real_T *Pn;
  real_T Sppc_data[392];
  real_T Sppcn_data[392];
  real_T Tppc_data[392];
  real_T Tppcn_data[392];
  real_T varargin_4_data[5];
  real_T b_p;
  real_T b_s;
  real_T d;
  real_T dxm;
  real_T dxp;
  real_T e;
  real_T fa;
  real_T fc;
  real_T lb;
  real_T m;
  real_T q;
  real_T r;
  real_T tol;
  real_T ub;
  int32_T Sppc_size[2];
  int32_T Sppcn_size[2];
  int32_T Sppc_idx_0;
  int32_T Sppc_idx_1;
  int32_T Tppc_idx_0;
  int32_T Tppc_idx_1;
  int32_T b_loop_ub;
  int32_T c_loop_ub;
  int32_T exitg1;
  int32_T i;
  int32_T i1;
  int32_T i2;
  int32_T loop_ub;
  int32_T n;
  boolean_T Pmat;
  boolean_T fapos;
  boolean_T fbpos;
  boolean_T guard1 = false;
  emlrtHeapReferenceStackEnterFcnR2012b(emlrtRootTLSGlobal);
  /* TBS_VERTSOLVE  Root finding of pressure or depth that matches equation of
   */
  /*                state with multivalued function. */
  /*  */
  /*  */
  /*  [p, s, t] = tbs_vertsolve(Sppc, Tppc, P, BotK, s, t, p, branch, d_fn, s0,
   * t0, tolp, DP) */
  /*  finds the pressure or depth p (within tolerance tolp) and its associated
   */
  /*  practical / Absolute salinity s and potential / Conservative temperature
   */
  /*  t, of a surface on which delta equals that determined by the multivalued
   */
  /*  function d_fn, in an ocean whose practical / Absolute salinity and */
  /*  potential / Conservative temperature as functions of pressure or depth P
   */
  /*  are given by piecewise polynomials whose coefficients are Sppc and Tppc,
   */
  /*  and whose knots are P.  The number of valid data points in each water */
  /*  column is given by BotK.  The equation of state is given by eos.m in the
   */
  /*  path, taking S, T, P as its 3 inputs. delta is the in-situ density */
  /*  anomaly or specific volume anomaly, defined as eos(s,t,p) - eos(s0,t0,p)
   */
  /*  where s,t are S,T interpolated from P to p, and s0, t0 are reference */
  /*  values.  In each water column I, the branch of the multivalued function */
  /*  for delta as a function of p is given by d_fn(:,b) where b = */
  /*  branchmap(I).  This branch is a polynomial of p between d_fn(1,b) and */
  /*  d_fn(2,b) with coefficients given by d_fn(3:end,b), and a linear */
  /*  extension of this polynomial outside this domain. Solutions are sought in
   */
  /*  the interval [A,B] where A = max(P(1,I), d_fn(1,b) - DP) and B = */
  /*  min(P(BotK(I),I), d_fn(2,b) + DP), i.e. within the valid P range of the */
  /*  local water column AND within DP of the polynomial domain of the local */
  /*  branch. (If P is a vector, the I indexing of P is dropped.)  Limiting */
  /*  this search domain with DP helps prevent the surface from jumping between
   */
  /*  multiple solutions in weakly stratified waters such as around Antarctica.
   */
  /*  The inputs s and t are not used, but provided so these variables may be */
  /*  manipulated in-place. */
  /*  */
  /*  */
  /*  --- Input: */
  /*  Sppc [O, K-1, N]: coefficients for piecewise polynomial for practical */
  /*                    / Absolute Salinity in terms of P */
  /*  Tppc [O, K-1, N]: coefficients for piecewise polynomial for potential */
  /*                    / Conservative Temperature in terms of P */
  /*  P [K, N]: knots for the pressure [dbar] or depth [m] of the casts */
  /*  BotK [1, N]: number of valid data points on each cast */
  /*  s [1, N]: initial practical / Absolute salinity on the initial surface */
  /*  t [1, N]: initial potential / Conservative temperature on the initial
   * surface */
  /*  p [1, N]: initial pressure [dbar] or depth [m] on the initial surface */
  /*  branchmap [1, N]: branch that each cast belongs to. */
  /*                    Must have 1 <= branch(i) <= A, for all 1 <= i <= N. */
  /*  d_fn [D+3, A]: domain and coefficients for each branch's polynomial */
  /*              for delta as a function of p, the b'th branch being defined */
  /*              by d_fn(:,b), which can be evaluated by pvaln or pvallin. */
  /*  s0 [1, 1]: reference S value for delta */
  /*  t0 [1, 1]: reference T value for delta */
  /*  tolp [1, 1]: tolerance on pressure [dbar] or depth [m] for root finding
   * solver */
  /*  DP [1, 1]: seek solutions in the domain of the local d_fn branch expanded
   */
  /*             by DP in both directions. */
  /*  */
  /*  Note: O is the order of the piecewise polynomials down each cast */
  /*        K is the maximum number of knots in these piecewise polynomials, */
  /*            i.e. the maximum number of bottles in any cast */
  /*        N is the number of water columns (possibly including land). */
  /*        A is the number of branches for the multivalued function d_fn. */
  /*        D is the degree of each branch's polynomial */
  /*  */
  /*  Note: P must increase along its first dimension. */
  /*  */
  /*  Note: P can have size [K, 1], in which case it is used for each cast. */
  /*  */
  /*  Note: variables can actually be higher dimensional, e.g. N = [ni, nj], */
  /*        and p can be any dimensional matrix, so long as it has N elements */
  /*        in total. */
  /*  */
  /*  Note: BotK should be given by */
  /*            BotK = squeeze(sum(isfinite(S), 1)); */
  /*  */
  /*  */
  /*  --- Output: */
  /*  p [same as input p]: pressure or depth of the updated surface */
  /*  s [same as input p]: practical / Absolute salinity of the updated surface
   */
  /*  t [same as input p]: potential / Conservative temperature of the updated
   * surface */
  /*  Author(s) : Geoff Stanley */
  /*  Email     : g.stanley@unsw.edu.au */
  /*  Email     : geoffstanley@gmail.com */
  Pmat = (((P->size[0] != 1) && (P->size[1] != 1)) || (P->size[2] != 1));
  /*  Loop over each cast */
  i = p->size[0] * p->size[1];
  emxInit_real_T(&Pn, 1, true);
  for (n = 0; n < i; n++) {
    if ((branchmap->data[n] > 0.0) && (BotK->data[n] > 1.0)) {
      /*  Select this water column */
      if (1.0 > BotK->data[n] - 1.0) {
        loop_ub = 0;
      } else {
        loop_ub = (int32_T)(BotK->data[n] - 1.0);
      }
      b_loop_ub = Sppc->size[0];
      Sppc_idx_0 = Sppc->size[0];
      Sppc_idx_1 = Sppc->size[1];
      Sppcn_size[0] = Sppc->size[0];
      Sppcn_size[1] = loop_ub;
      for (i1 = 0; i1 < loop_ub; i1++) {
        for (i2 = 0; i2 < b_loop_ub; i2++) {
          Sppcn_data[i2 + Sppcn_size[0] * i1] =
              Sppc->data[(i2 + Sppc_idx_0 * i1) + Sppc_idx_0 * Sppc_idx_1 * n];
        }
      }
      if (1.0 > BotK->data[n] - 1.0) {
        b_loop_ub = 0;
      } else {
        b_loop_ub = (int32_T)(BotK->data[n] - 1.0);
      }
      c_loop_ub = Tppc->size[0];
      Tppc_idx_0 = Tppc->size[0];
      Tppc_idx_1 = Tppc->size[1];
      Sppc_idx_1 = Tppc->size[0];
      for (i1 = 0; i1 < b_loop_ub; i1++) {
        for (i2 = 0; i2 < c_loop_ub; i2++) {
          Tppcn_data[i2 + Sppc_idx_1 * i1] =
              Tppc->data[(i2 + Tppc_idx_0 * i1) + Tppc_idx_0 * Tppc_idx_1 * n];
        }
      }
      if (Pmat) {
        Sppc_idx_1 = P->size[0];
        Sppc_idx_0 = (int32_T)BotK->data[n];
        i1 = Pn->size[0];
        Pn->size[0] = Sppc_idx_0;
        emxEnsureCapacity_real_T(Pn, i1);
        for (i1 = 0; i1 < Sppc_idx_0; i1++) {
          Pn->data[i1] = P->data[i1 + Sppc_idx_1 * n];
        }
      } else {
        c_loop_ub = (int32_T)muDoubleScalarFloor(BotK->data[n] - 1.0);
        i1 = Pn->size[0];
        Pn->size[0] = c_loop_ub + 1;
        emxEnsureCapacity_real_T(Pn, i1);
        for (i1 = 0; i1 <= c_loop_ub; i1++) {
          Pn->data[i1] = P->data[i1];
        }
        /*  .' is for codegen, so P and (1:k).' both column vectors */
      }
      fa = muDoubleScalarMax(
          Pn->data[0],
          d_fn->data[d_fn->size[0] * ((int32_T)branchmap->data[n] - 1)] - DP);
      b_p = muDoubleScalarMin(
          Pn->data[(int32_T)BotK->data[n] - 1],
          d_fn->data[d_fn->size[0] * ((int32_T)branchmap->data[n] - 1) + 1] +
              DP);
      /*  Search for a sign-change, expanding outward from an initial guess */
      q = p->data[n];
      c_loop_ub = d_fn->size[0];
      Sppc_idx_1 = d_fn->size[0];
      for (i1 = 0; i1 < c_loop_ub; i1++) {
        varargin_4_data[i1] =
            d_fn->data[i1 + d_fn->size[0] * ((int32_T)branchmap->data[n] - 1)];
      }
      /* FZERO_GUESS_TO_BOUNDS  Search for a sign change bounding a zero of a */
      /*                        univariate function, expanding geometrically */
      /*                        outward from an initial guess. */
      /*  */
      /*  */
      /*  [a, b] = fzero_guess_to_bounds(f, x) */
      /*  finds a < b such that f(a) and f(b) have different sign*, meaning a */
      /*  solution exists within the interval [a,b].  The bounds a,b are
       * expanded */
      /*  outward in geometric progression from an initial guess for the root of
       * f */
      /*  at x. If f evaluates to NaN at any point during the search, then a =
       * nan */
      /*  and b = nan are immediately returned.  If the function is genuinely */
      /*  single-signed, or even if it is not but its values of opposite sign
       * are */
      /*  skipped over, it is possible to enter an infinite loop.  Calling the
       */
      /*  function in this form is therefore not recommended unless you know the
       */
      /*  function will not result in such an infinite loop. */
      /*  */
      /*  [a, b] = fzero_guess_to_bounds(f, x, A, B) */
      /*  as above, but limits [a,b] to remain inside the subset [A, B].  If x
       * is */
      /*  outside of [A, B], it is immediately moved into this range. If no */
      /*  sign-change is found within [A, B], then a = nan and b = nan are */
      /*  returned.  Note, as above, it is possible that a sign-change is
       * skipped */
      /*  over as f is only evaluated at finitely many x values. */
      /*  */
      /*  [a,b] = fzero_guess_to_bounds(f, x, A, B, ...) */
      /*  passes all additional arguments to the function f.  */
      /*  */
      /*  * Note: for computational speed, herein the "sign" of 0 is considered
       * the */
      /*  same as the sign of a negative number. */
      /*  */
      /*  This function is compatible with MATLAB's code generation. */
      /*  */
      /*  */
      /*  --- Input: */
      /*    f       : handle to a function that accepts a real scalar as its
       * first */
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
      /*            : 22/07/2020 - fix infinite loop in bounded case, arising
       * from machine precision rounding */
      /*  Geometrically expand from the guess x, until a sign change is found */
      /*  Handle bad inputs */
      fapos = muDoubleScalarIsNaN(fa);
      if (fapos || muDoubleScalarIsNaN(b_p) || muDoubleScalarIsNaN(q)) {
        lb = rtNaN;
        ub = rtNaN;
      } else if (myfcn(fa, Sppcn_data, Sppcn_size, Tppcn_data, Pn,
                       varargin_4_data, Sppc_idx_1, s0, t0) == 0.0) {
        lb = fa;
        ub = fa;
      } else {
        r = myfcn(b_p, Sppcn_data, Sppcn_size, Tppcn_data, Pn, varargin_4_data,
                  Sppc_idx_1, s0, t0);
        if (r == 0.0) {
          lb = b_p;
          ub = b_p;
        } else {
          q = muDoubleScalarMin(muDoubleScalarMax(q, fa), b_p);
          /*  bounds are given */
          dxp = (b_p - q) / 50.0;
          dxm = (q - fa) / 50.0;
          /*  Set a = x, except when x is so close to A that machine roundoff
           * makes dxm identically 0, */
          /*  which would lead to an infinite loop below.  In this case, set a =
           * A. */
          if (dxm == 0.0) {
            lb = fa;
          } else {
            lb = q;
          }
          fapos = (myfcn(lb, Sppcn_data, Sppcn_size, Tppcn_data, Pn,
                         varargin_4_data, Sppc_idx_1, s0, t0) > 0.0);
          /*  Similarly, set b = x, except for machine precision problems. */
          if (dxp == 0.0) {
            ub = b_p;
            fbpos = (r > 0.0);
          } else {
            ub = q;
            if (dxm == 0.0) {
              fbpos = fapos;
              /*  since a = b = x */
            } else {
              fbpos = (myfcn(q, Sppcn_data, Sppcn_size, Tppcn_data, Pn,
                             varargin_4_data, Sppc_idx_1, s0, t0) > 0.0);
            }
          }
          do {
            exitg1 = 0;
            guard1 = false;
            if (lb > fa) {
              /*  Move a left, and test for a sign change */
              dxm *= 1.4142135623730949;
              lb = muDoubleScalarMax(q - dxm, fa);
              fapos = (myfcn(lb, Sppcn_data, Sppcn_size, Tppcn_data, Pn,
                             varargin_4_data, Sppc_idx_1, s0, t0) > 0.0);
              if (fapos != fbpos) {
                /*  fa and fb have different signs */
                exitg1 = 1;
              } else {
                guard1 = true;
              }
            } else if (ub == b_p) {
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
              if (ub < b_p) {
                /*  Move b right, and test for a sign change */
                dxp *= 1.4142135623730949;
                ub = muDoubleScalarMin(q + dxp, b_p);
                fbpos = (myfcn(ub, Sppcn_data, Sppcn_size, Tppcn_data, Pn,
                               varargin_4_data, Sppc_idx_1, s0, t0) > 0.0);
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
      if (!muDoubleScalarIsNaN(lb)) {
        /*  A sign change was discovered, so a root exists in the interval. */
        /*  Solve the nonlinear root-finding problem using Brent's method */
        /* FZERO_BRENT  Find a root of a univariate function within a given
         * interval */
        /*              using Brent's method */
        /*  */
        /*  x = fzero_brent(f,a,b,t)  */
        /*  finds x in the interval [a,b] satisfying |x-y| <= t/2 where f(y) =
         * 0. */
        /*  f(a) and f(b) must have opposite signs. */
        /*  */
        /*  ... = fzero_brent(f,lb,ub,t,x,...) passes additional inputs directly
         */
        /*  to f. */
        /*  */
        /*  This function is compatible with MATLAB's code generation -- so long
         * f is */
        /*  similarly compatible. Many root-finding problems can be solved with
         */
        /*  fzero_brent by writing another function which calls fzero_brent
         * inside a */
        /*  for loop, and this can be made fast by using code generation on that
         */
        /*  wrapper function. Note that f will only be called with a scalar
         * input as */
        /*  its first argument; codegen knows this, and might strip out
         * unnecessary */
        /*  code from the function definition underlying f. */
        /*  --- Example: */
        /*  Simple bisection between bounds of -0.5 and 1.5 would fail to find
         * the */
        /*  root. Starting with a guess of .85 and expanding outwards finds a
         * root. */
        /*  This example is shown graphically on the MATLAB File Exchange. */
        /*  */
        /*  c = poly([0 1]); % a polynomial with roots at 0 and 1 */
        /*  f = @(x) polyval(c,x); % encapsulate extra parameters */
        /*  [a,b] = fzero_guess_to_bounds(f, .85, -.5, 1.5); */
        /*  root = fzero_brent(f, a, b, .05); */
        /*  */
        /*   Discussion: */
        /*  */
        /*     The interval [A,B] must be a change of sign interval for F. */
        /*     That is, F(A) and F(B) must be of opposite signs.  Then */
        /*     assuming that F is continuous implies the existence of at least
         */
        /*     one value C between A and B for which F(C) = 0. */
        /*  */
        /*     The location of the zero is determined to within an accuracy */
        /*     of 6 * EPS * abs ( C ) + 2 * T, where EPS is the machine epsilon.
         */
        /*  */
        /*     Thanks to Thomas Secretin for pointing out a transcription error
         * in the */
        /*     setting of the value of P, 11 February 2013. */
        /*  */
        /*     Additional parameters given by varargin are passed directly to F.
         */
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
        /*     Minor changes (passing extra arguments) by Geoff Stanley. */
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
        /*     Input, real A, B, the endpoints of the change of sign interval.
         */
        /*  */
        /*     Input, real T, a positive error tolerance. */
        /*  */
        /*     Input, real value = F ( x ), the name of a user-supplied */
        /*     function which evaluates the function whose zero is being sought.
         */
        /*  */
        /*     Output, real VALUE, the estimated value of a zero of */
        /*     the function F. */
        /*  */
        c_loop_ub = Sppc->size[0];
        Sppc_idx_0 = Sppc->size[0];
        Sppc_idx_1 = Sppc->size[1];
        Tppc_idx_0 = Tppc->size[0];
        Tppc_idx_1 = Tppc->size[1];
        Sppc_size[0] = Sppc->size[0];
        Sppc_size[1] = loop_ub;
        for (i1 = 0; i1 < loop_ub; i1++) {
          for (i2 = 0; i2 < c_loop_ub; i2++) {
            Sppc_data[i2 + Sppc_size[0] * i1] =
                Sppc->data[(i2 + Sppc_idx_0 * i1) +
                           Sppc_idx_0 * Sppc_idx_1 * n];
          }
        }
        c_loop_ub = Tppc->size[0];
        Sppc_idx_0 = Tppc->size[0];
        for (i1 = 0; i1 < b_loop_ub; i1++) {
          for (i2 = 0; i2 < c_loop_ub; i2++) {
            Tppc_data[i2 + Sppc_idx_0 * i1] =
                Tppc->data[(i2 + Tppc_idx_0 * i1) +
                           Tppc_idx_0 * Tppc_idx_1 * n];
          }
        }
        c_loop_ub = d_fn->size[0];
        Sppc_idx_1 = d_fn->size[0];
        for (i1 = 0; i1 < c_loop_ub; i1++) {
          varargin_4_data[i1] =
              d_fn->data[i1 +
                         d_fn->size[0] * ((int32_T)branchmap->data[n] - 1)];
        }
        fa = myfcn(lb, Sppc_data, Sppc_size, Tppc_data, Pn, varargin_4_data,
                   Sppc_idx_1, s0, t0);
        c_loop_ub = Sppc->size[0];
        Sppc_idx_0 = Sppc->size[0];
        Sppc_idx_1 = Sppc->size[1];
        Tppc_idx_0 = Tppc->size[0];
        Tppc_idx_1 = Tppc->size[1];
        Sppc_size[0] = Sppc->size[0];
        Sppc_size[1] = loop_ub;
        for (i1 = 0; i1 < loop_ub; i1++) {
          for (i2 = 0; i2 < c_loop_ub; i2++) {
            Sppc_data[i2 + Sppc_size[0] * i1] =
                Sppc->data[(i2 + Sppc_idx_0 * i1) +
                           Sppc_idx_0 * Sppc_idx_1 * n];
          }
        }
        c_loop_ub = Tppc->size[0];
        Sppc_idx_0 = Tppc->size[0];
        for (i1 = 0; i1 < b_loop_ub; i1++) {
          for (i2 = 0; i2 < c_loop_ub; i2++) {
            Tppc_data[i2 + Sppc_idx_0 * i1] =
                Tppc->data[(i2 + Tppc_idx_0 * i1) +
                           Tppc_idx_0 * Tppc_idx_1 * n];
          }
        }
        c_loop_ub = d_fn->size[0];
        Sppc_idx_1 = d_fn->size[0];
        for (i1 = 0; i1 < c_loop_ub; i1++) {
          varargin_4_data[i1] =
              d_fn->data[i1 +
                         d_fn->size[0] * ((int32_T)branchmap->data[n] - 1)];
        }
        dxp = myfcn(ub, Sppc_data, Sppc_size, Tppc_data, Pn, varargin_4_data,
                    Sppc_idx_1, s0, t0);
        dxm = lb;
        fc = fa;
        e = ub - lb;
        d = e;
        do {
          exitg1 = 0;
          if (muDoubleScalarAbs(fc) < muDoubleScalarAbs(dxp)) {
            lb = ub;
            ub = dxm;
            dxm = lb;
            fa = dxp;
            dxp = fc;
            fc = fa;
          }
          tol = 4.4408920985006262E-16 * muDoubleScalarAbs(ub) + tolp;
          m = 0.5 * (dxm - ub);
          if ((muDoubleScalarAbs(m) <= tol) || (dxp == 0.0)) {
            exitg1 = 1;
          } else {
            if ((muDoubleScalarAbs(e) < tol) ||
                (muDoubleScalarAbs(fa) <= muDoubleScalarAbs(dxp))) {
              e = m;
              d = m;
            } else {
              b_s = dxp / fa;
              if (lb == dxm) {
                b_p = 2.0 * m * b_s;
                q = 1.0 - b_s;
              } else {
                q = fa / fc;
                r = dxp / fc;
                b_p = b_s * (2.0 * m * q * (q - r) - (ub - lb) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (b_s - 1.0);
              }
              if (0.0 < b_p) {
                q = -q;
              } else {
                b_p = -b_p;
              }
              b_s = e;
              e = d;
              if ((2.0 * b_p < 3.0 * m * q - muDoubleScalarAbs(tol * q)) &&
                  (b_p < muDoubleScalarAbs(0.5 * b_s * q))) {
                d = b_p / q;
              } else {
                e = m;
                d = m;
              }
            }
            lb = ub;
            fa = dxp;
            if (tol < muDoubleScalarAbs(d)) {
              ub += d;
            } else if (0.0 < m) {
              ub += tol;
            } else {
              ub -= tol;
            }
            c_loop_ub = d_fn->size[0];
            Sppc_idx_1 = d_fn->size[0];
            for (i1 = 0; i1 < c_loop_ub; i1++) {
              varargin_4_data[i1] =
                  d_fn->data[i1 +
                             d_fn->size[0] * ((int32_T)branchmap->data[n] - 1)];
            }
            dxp = myfcn(ub, Sppcn_data, Sppcn_size, Tppcn_data, Pn,
                        varargin_4_data, Sppc_idx_1, s0, t0);
            if (((0.0 < dxp) && (0.0 < fc)) || ((dxp <= 0.0) && (fc <= 0.0))) {
              dxm = lb;
              fc = fa;
              e = ub - lb;
              d = e;
            }
          }
        } while (exitg1 == 0);
        p->data[n] = ub;
        /*  Interpolate S and T onto the updated surface */
        c_loop_ub = Sppc->size[0];
        Sppc_idx_0 = Sppc->size[0];
        Sppc_idx_1 = Sppc->size[1];
        Tppc_idx_0 = Tppc->size[0];
        Tppc_idx_1 = Tppc->size[1];
        Sppc_size[0] = Sppc->size[0];
        Sppc_size[1] = loop_ub;
        for (i1 = 0; i1 < loop_ub; i1++) {
          for (i2 = 0; i2 < c_loop_ub; i2++) {
            Sppc_data[i2 + Sppc_size[0] * i1] =
                Sppc->data[(i2 + Sppc_idx_0 * i1) +
                           Sppc_idx_0 * Sppc_idx_1 * n];
          }
        }
        loop_ub = Tppc->size[0];
        Sppc_idx_0 = Tppc->size[0];
        for (i1 = 0; i1 < b_loop_ub; i1++) {
          for (i2 = 0; i2 < loop_ub; i2++) {
            Tppc_data[i2 + Sppc_idx_0 * i1] =
                Tppc->data[(i2 + Tppc_idx_0 * i1) +
                           Tppc_idx_0 * Tppc_idx_1 * n];
          }
        }
        ppc_val2(Pn, Sppc_data, Sppc_size, Tppc_data, p->data[n], &s->data[n],
                 &t->data[n]);
      } else {
        p->data[n] = rtNaN;
        s->data[n] = rtNaN;
        t->data[n] = rtNaN;
      }
    } else {
      /*  This will ensure s,t,p all have the same nan structure */
      p->data[n] = rtNaN;
      s->data[n] = rtNaN;
      t->data[n] = rtNaN;
    }
  }
  emxFree_real_T(&Pn);
  emlrtHeapReferenceStackLeaveFcnR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (tbs_vertsolve.c) */
