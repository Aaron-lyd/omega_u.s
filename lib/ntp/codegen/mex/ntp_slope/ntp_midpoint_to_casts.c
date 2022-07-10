/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * ntp_midpoint_to_casts.c
 *
 * Code generation for function 'ntp_midpoint_to_casts'
 *
 */

/* Include files */
#include "ntp_midpoint_to_casts.h"
#include "ppc_val2.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"

/* Function Definitions */
real_T
ntp_midpoint_to_casts(const real_T Sppc_A_data[], const int32_T Sppc_A_size[2],
                      const real_T Tppc_A_data[], const real_T P_A_data[],
                      int32_T P_A_size, real_T k_A, const real_T Sppc_B_data[],
                      const int32_T Sppc_B_size[2], const real_T Tppc_B_data[],
                      const real_T P_B_data[], int32_T P_B_size, real_T k_B,
                      real_T p_A, real_T p_B, real_T tolp)
{
  real_T b_s1o2;
  real_T c;
  real_T d;
  real_T dp;
  real_T dxm;
  real_T dxp;
  real_T e;
  real_T fa;
  real_T fc;
  real_T lb;
  real_T m;
  real_T p_avg;
  real_T q;
  real_T r;
  real_T s;
  real_T s1o2;
  real_T s_B;
  real_T t_B;
  real_T z_tmp;
  int32_T exitg1;
  boolean_T fapos;
  boolean_T fbpos;
  boolean_T guard1 = false;
  /* NTP_MIDPOINT_TO_CASTS  Find the Neutral Tangent Plane from a midpoint */
  /*                        between a pair of water columns */
  /*  */
  /*  */
  /*  dp = ntp_midpoint_to_casts(Sppc_A, Tppc_A, P_A, k_A, Sppc_B, Tppc_B, P_B,
   * k_B, p_A, p_B, tolp) */
  /*  finds the P difference, dp, that satisfy, with accuracy tolp, */
  /*      eos(S_A(p_avg + dp / 2), T_A(p_avg + dp / 2), p_avg) */
  /*    = eos(S_B(p_avg - dp / 2), T_B(p_avg - dp / 2), p_avg) */
  /*  where */
  /*    S_A, T_A are the S, T profiles on cast A, */
  /*    S_B, T_B are the S, T profiles on cast B, */
  /*    p_avg = (p_A + p_B) / 2  is the midpoint z. */
  /*  Here, S_A and T_A have been written as functions of p, but more */
  /*  specifically they are piecewise polynomial interpolants with coefficients
   */
  /*  given by Sppc_A and Tppc_A with knots at P_A(1:k_A).  That is, */
  /*    [S_A(z), T_A(z)] = ppc_val2(P_A, Sppc_A, Tppc_A, z). */
  /*  Similarly for S_B and T_B. */
  /*  The equation of state is given by eos.m on the MATLAB path. */
  /*  */
  /*  */
  /*  --- Input: */
  /*  Sppc_A, Sppc_B [O, K-1]: coefficients for piecewise polynomial for */
  /*      practical / Absolute Salinity in terms of P, on casts A and B.  Here,
   */
  /*      O is the polynomial order and K is the number of data points */
  /*      (including NaN's) on the casts. */
  /*  Tppc_A, Tppc_B [O, K-1]: as above, but for potential / Conservative */
  /*      Temperature. */
  /*  P_A, P_B [K, 1]: pressure or depth data points on casts A and B. */
  /*  k_A, k_B [1, 1]: number of valid (non-NaN) salinity / temperature data  */
  /*    points on casts A and B.  Specifically, if S_A were the salinity data */
  /*    on cast A, then S_A(1:k_A) must be non-NaN.  This translates into */
  /*    requiring Sppc_A(:,1:k_A-1) be non-NaN.  Similarly for the others. */
  /*  tolp [1, 1]: tolerance for solving the level of neutral buoyancy (same */
  /*      units as P). */
  /*  p_A, p_B [1, 1]: initial guess for the pressure or depth that the NTP */
  /*      intersects casts A and B.  For example, the pressure or depth of a */
  /*      certain potential density surface.  Though this is used as an initial
   */
  /*      guess for the non-linear root finding solver, the main purpose of p_A
   */
  /*      and p_B is to define the average pressure or depth, p_avg. */
  /*  */
  /*  Note: physical units for Sppc_*, Tppc_*, P_*, p_*  are determined by
   * eos.m. */
  /*  */
  /*  Note: P_A and P_B must increase monotonically. */
  /*  */
  /*  */
  /*  --- Output: */
  /*  dp: pressure or depth difference between where NTP intersects casts A and
   */
  /*      B.  Specifically, the NTP intersects cast A at p_avg + dp / 2 and */
  /*      cast B at p_avg - dp / 2. */
  /*  */
  /*  */
  /*  --- See Also: */
  /*  ntp_bottle_to_cast */
  /*  ntp_slope */
  /*  ntp_slope_error */
  /*  ppc_linterp, ppc_pchip */
  /*  Author(s) : Geoff Stanley */
  /*  Email     : g.stanley@unsw.edu.au */
  /*  Email     : geoffstanley@gmail.com */
  if ((k_A > 1.0) && (k_B > 1.0)) {
    /*  Casts have at least two data points.  A solution is possible. */
    p_avg = (p_A + p_B) / 2.0;
    fa = p_A - p_B;
    /*  Upper and lower bounds of dp must satisfy */
    /*  (a) P_A(1) <= p_avg + dp/2 <= P_A(k_A) */
    /*  (b) P_B(1) <= p_avg - dp/2 <= P_B(k_B) */
    /*  ... therefore ... */
    /*     max( P_A(1) - p_avg, p_avg - P_B(k_B) ) */
    /*  <= dp/2 */
    /*  <= min( P_A(k_A) - p_avg, p_avg - P_B(1) ) */
    q = 2.0 * muDoubleScalarMax(P_A_data[0] - p_avg,
                                p_avg - P_B_data[(int32_T)k_B - 1]);
    r = 2.0 * muDoubleScalarMin(P_A_data[(int32_T)k_A - 1] - p_avg,
                                p_avg - P_B_data[0]);
    /*  Search for a sign change, expanding outward from an initial guess */
    /* FZERO_GUESS_TO_BOUNDS  Search for a sign change bounding a zero of a */
    /*                        univariate function, expanding geometrically */
    /*                        outward from an initial guess. */
    /*  */
    /*  */
    /*  [a, b] = fzero_guess_to_bounds(f, x) */
    /*  finds a < b such that f(a) and f(b) have different sign*, meaning a */
    /*  solution exists within the interval [a,b].  The bounds a,b are expanded
     */
    /*  outward in geometric progression from an initial guess for the root of f
     */
    /*  at x. If f evaluates to NaN at any point during the search, then a = nan
     */
    /*  and b = nan are immediately returned.  If the function is genuinely */
    /*  single-signed, or even if it is not but its values of opposite sign are
     */
    /*  skipped over, it is possible to enter an infinite loop.  Calling the */
    /*  function in this form is therefore not recommended unless you know the
     */
    /*  function will not result in such an infinite loop. */
    /*  */
    /*  [a, b] = fzero_guess_to_bounds(f, x, A, B) */
    /*  as above, but limits [a,b] to remain inside the subset [A, B].  If x is
     */
    /*  outside of [A, B], it is immediately moved into this range. If no */
    /*  sign-change is found within [A, B], then a = nan and b = nan are */
    /*  returned.  Note, as above, it is possible that a sign-change is skipped
     */
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
    /*    f       : handle to a function that accepts a real scalar as its first
     */
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
    /*            : 22/07/2020 - fix infinite loop in bounded case, arising from
     * machine precision rounding */
    /*  Geometrically expand from the guess x, until a sign change is found */
    /*  Handle bad inputs */
    fapos = muDoubleScalarIsNaN(q);
    if (fapos || muDoubleScalarIsNaN(r) || muDoubleScalarIsNaN(fa)) {
      lb = rtNaN;
      dp = rtNaN;
    } else {
      /*  Evaluate difference between (a) eos at location on cast A with */
      /*  properties (S_A, T_A, P_A) where the pressure or depth is p_avg +
       * dp/2, */
      /*  and (b) eos at location on cast B with properties (S_B, T_B, P_B)
       * where */
      /*  the pressure or depth is p_avg - dp/2. Here, eos is always evaluated
       * at */
      /*  the mid-pressure or mid-depth, p_avg. */
      ppc_val2(P_A_data, P_A_size, Sppc_A_data, Sppc_A_size, Tppc_A_data,
               p_avg + 0.5 * q, &m, &s);
      ppc_val2(P_B_data, P_B_size, Sppc_B_data, Sppc_B_size, Tppc_B_data,
               p_avg - 0.5 * q, &s_B, &t_B);
      /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
      /*  */
      /*  */
      /*  rho = eoscg_densjmd95_bsq(s,t,z)                             [ kg /
       * m^3 ] */
      /*  computes the Boussinesq JMD95 in-situ density, given practical
       * salinity */
      /*  s, potential temperature t, and depth z [m, positive and increasing */
      /*  down].  The depth z is converted into hydrostatic pressure for a given
       */
      /*  gravitational acceleration and Boussinesq reference density, which are
       */
      /*  hard-coded into this function (edit variables grav and rhob as
       * needed). */
      /*  */
      /*  */
      /*  This function is derived from densjmd95.m, documented below. Input
       * checks */
      /*  and expansion of variables have been removed; instead automatic
       * expansion */
      /*  (requiring MATLAB 2016b or later) is used. The calculation has also
       * been */
      /*  streamlined by pre-allocating arrays for coeffcients, and modifying
       * the */
      /*  coefficients such that pressure need not be internally converted from
       */
      /*  dbar to bar. This function is compatible with MATLAB's codegen. */
      /*  */
      /*  The units of s and t are as in densjmd95, documented below. */
      /*   */
      /*  The inputs' sizes are more general: they may be arrays of any
       * dimension */
      /*  and any size, so long as their sizes match each other excepting that
       * any */
      /*  input can have any dimension as a singleton. For example, s and t can
       * be */
      /*  3D arrays of size [nz, nx, ny], while p can be a vector of size
       * [nz,1]. */
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
      z_tmp = p_avg * 1.015335;
      /*  Henceforth z is actually pressure [dbar] */
      /*  coefficients nonlinear equation of state in pressure coordinates for
       */
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
      s1o2 = muDoubleScalarSqrt(m);
      /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
      /*  */
      /*  */
      /*  rho = eoscg_densjmd95_bsq(s,t,z)                             [ kg /
       * m^3 ] */
      /*  computes the Boussinesq JMD95 in-situ density, given practical
       * salinity */
      /*  s, potential temperature t, and depth z [m, positive and increasing */
      /*  down].  The depth z is converted into hydrostatic pressure for a given
       */
      /*  gravitational acceleration and Boussinesq reference density, which are
       */
      /*  hard-coded into this function (edit variables grav and rhob as
       * needed). */
      /*  */
      /*  */
      /*  This function is derived from densjmd95.m, documented below. Input
       * checks */
      /*  and expansion of variables have been removed; instead automatic
       * expansion */
      /*  (requiring MATLAB 2016b or later) is used. The calculation has also
       * been */
      /*  streamlined by pre-allocating arrays for coeffcients, and modifying
       * the */
      /*  coefficients such that pressure need not be internally converted from
       */
      /*  dbar to bar. This function is compatible with MATLAB's codegen. */
      /*  */
      /*  The units of s and t are as in densjmd95, documented below. */
      /*   */
      /*  The inputs' sizes are more general: they may be arrays of any
       * dimension */
      /*  and any size, so long as their sizes match each other excepting that
       * any */
      /*  input can have any dimension as a singleton. For example, s and t can
       * be */
      /*  3D arrays of size [nz, nx, ny], while p can be a vector of size
       * [nz,1]. */
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
      /*  coefficients nonlinear equation of state in pressure coordinates for
       */
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
      b_s1o2 = muDoubleScalarSqrt(s_B);
      if (((s * (s * (s * (s * (s * 6.536332E-9 + -1.120083E-6) +
                           0.0001001685) +
                      -0.00909529) +
                 0.06793952) +
            999.842594) +
           m * (((s * (s * (s * (s * 5.3875E-9 + -8.2467E-7) + 7.6438E-5) +
                       -0.0040899) +
                  0.824493) +
                 s1o2 * (s * (s * -1.6546E-6 + 0.00010227) + -0.00572466)) +
                m * 0.00048314)) /
                  (1.0 -
                   z_tmp /
                       (((s * (s * (s * (s * -0.0004190253 + 0.09648704) +
                                    -17.06103) +
                               1444.304) +
                          196593.3) +
                         m * ((s * (s * (s * -0.0005084188 + 0.06283263) +
                                    -3.101089) +
                               528.4855) +
                              s1o2 * (s * (s * -0.004619924 + 0.09085835) +
                                      3.88664))) +
                        z_tmp * (((s * (s * (s * 1.956415E-6 + -0.0002984642) +
                                        0.02212276) +
                                   3.186519) +
                                  m * ((s * (s * 2.059331E-7 + -0.0001847318) +
                                        0.006704388) +
                                       s1o2 * 0.0001480266)) +
                                 z_tmp * ((s * (s * 1.39468E-8 + -1.202016E-6) +
                                           2.102898E-5) +
                                          m * (s * (s * 6.207323E-11 +
                                                    6.128773E-9) +
                                               -2.040237E-7))))) -
              ((t_B * (t_B * (t_B * (t_B * (t_B * 6.536332E-9 + -1.120083E-6) +
                                     0.0001001685) +
                              -0.00909529) +
                       0.06793952) +
                999.842594) +
               s_B * (((t_B * (t_B * (t_B * (t_B * 5.3875E-9 + -8.2467E-7) +
                                      7.6438E-5) +
                               -0.0040899) +
                        0.824493) +
                       b_s1o2 * (t_B * (t_B * -1.6546E-6 + 0.00010227) +
                                 -0.00572466)) +
                      s_B * 0.00048314)) /
                  (1.0 -
                   z_tmp /
                       (((t_B *
                              (t_B * (t_B * (t_B * -0.0004190253 + 0.09648704) +
                                      -17.06103) +
                               1444.304) +
                          196593.3) +
                         s_B *
                             ((t_B * (t_B * (t_B * -0.0005084188 + 0.06283263) +
                                      -3.101089) +
                               528.4855) +
                              b_s1o2 *
                                  (t_B * (t_B * -0.004619924 + 0.09085835) +
                                   3.88664))) +
                        z_tmp *
                            (((t_B *
                                   (t_B * (t_B * 1.956415E-6 + -0.0002984642) +
                                    0.02212276) +
                               3.186519) +
                              s_B *
                                  ((t_B * (t_B * 2.059331E-7 + -0.0001847318) +
                                    0.006704388) +
                                   b_s1o2 * 0.0001480266)) +
                             z_tmp * ((t_B * (t_B * 1.39468E-8 + -1.202016E-6) +
                                       2.102898E-5) +
                                      s_B * (t_B * (t_B * 6.207323E-11 +
                                                    6.128773E-9) +
                                             -2.040237E-7))))) ==
          0.0) {
        lb = q;
        dp = q;
      } else {
        /*  Evaluate difference between (a) eos at location on cast A with */
        /*  properties (S_A, T_A, P_A) where the pressure or depth is p_avg +
         * dp/2, */
        /*  and (b) eos at location on cast B with properties (S_B, T_B, P_B)
         * where */
        /*  the pressure or depth is p_avg - dp/2. Here, eos is always evaluated
         * at */
        /*  the mid-pressure or mid-depth, p_avg. */
        ppc_val2(P_A_data, P_A_size, Sppc_A_data, Sppc_A_size, Tppc_A_data,
                 p_avg + 0.5 * r, &m, &s);
        ppc_val2(P_B_data, P_B_size, Sppc_B_data, Sppc_B_size, Tppc_B_data,
                 p_avg - 0.5 * r, &s_B, &t_B);
        /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
        /*  */
        /*  */
        /*  rho = eoscg_densjmd95_bsq(s,t,z)                             [ kg /
         * m^3 ] */
        /*  computes the Boussinesq JMD95 in-situ density, given practical
         * salinity */
        /*  s, potential temperature t, and depth z [m, positive and increasing
         */
        /*  down].  The depth z is converted into hydrostatic pressure for a
         * given */
        /*  gravitational acceleration and Boussinesq reference density, which
         * are */
        /*  hard-coded into this function (edit variables grav and rhob as
         * needed). */
        /*  */
        /*  */
        /*  This function is derived from densjmd95.m, documented below. Input
         * checks */
        /*  and expansion of variables have been removed; instead automatic
         * expansion */
        /*  (requiring MATLAB 2016b or later) is used. The calculation has also
         * been */
        /*  streamlined by pre-allocating arrays for coeffcients, and modifying
         * the */
        /*  coefficients such that pressure need not be internally converted
         * from */
        /*  dbar to bar. This function is compatible with MATLAB's codegen. */
        /*  */
        /*  The units of s and t are as in densjmd95, documented below. */
        /*   */
        /*  The inputs' sizes are more general: they may be arrays of any
         * dimension */
        /*  and any size, so long as their sizes match each other excepting that
         * any */
        /*  input can have any dimension as a singleton. For example, s and t
         * can be */
        /*  3D arrays of size [nz, nx, ny], while p can be a vector of size
         * [nz,1]. */
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
        /*     Density of Sea Water using Jackett and McDougall 1995 (JAOT 12)
         */
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
        /*  Convert depth [m] to the hydrostatic pressure [dbar] implied by the
         */
        /*  following two hard-coded parameters (edit as needed): */
        /*  gravitational acceleration [m /s^2] */
        /*  Boussinesq reference density [kg / m^3] */
        /*  Pascal to dbar conversion [dbar / Pa] */
        /*  depth to pressure conversion [dbar / m] */
        /*  Henceforth z is actually pressure [dbar] */
        /*  coefficients nonlinear equation of state in pressure coordinates for
         */
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
        s1o2 = muDoubleScalarSqrt(m);
        /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
        /*  */
        /*  */
        /*  rho = eoscg_densjmd95_bsq(s,t,z)                             [ kg /
         * m^3 ] */
        /*  computes the Boussinesq JMD95 in-situ density, given practical
         * salinity */
        /*  s, potential temperature t, and depth z [m, positive and increasing
         */
        /*  down].  The depth z is converted into hydrostatic pressure for a
         * given */
        /*  gravitational acceleration and Boussinesq reference density, which
         * are */
        /*  hard-coded into this function (edit variables grav and rhob as
         * needed). */
        /*  */
        /*  */
        /*  This function is derived from densjmd95.m, documented below. Input
         * checks */
        /*  and expansion of variables have been removed; instead automatic
         * expansion */
        /*  (requiring MATLAB 2016b or later) is used. The calculation has also
         * been */
        /*  streamlined by pre-allocating arrays for coeffcients, and modifying
         * the */
        /*  coefficients such that pressure need not be internally converted
         * from */
        /*  dbar to bar. This function is compatible with MATLAB's codegen. */
        /*  */
        /*  The units of s and t are as in densjmd95, documented below. */
        /*   */
        /*  The inputs' sizes are more general: they may be arrays of any
         * dimension */
        /*  and any size, so long as their sizes match each other excepting that
         * any */
        /*  input can have any dimension as a singleton. For example, s and t
         * can be */
        /*  3D arrays of size [nz, nx, ny], while p can be a vector of size
         * [nz,1]. */
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
        /*     Density of Sea Water using Jackett and McDougall 1995 (JAOT 12)
         */
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
        /*  Convert depth [m] to the hydrostatic pressure [dbar] implied by the
         */
        /*  following two hard-coded parameters (edit as needed): */
        /*  gravitational acceleration [m /s^2] */
        /*  Boussinesq reference density [kg / m^3] */
        /*  Pascal to dbar conversion [dbar / Pa] */
        /*  depth to pressure conversion [dbar / m] */
        /*  Henceforth z is actually pressure [dbar] */
        /*  coefficients nonlinear equation of state in pressure coordinates for
         */
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
        b_s1o2 = muDoubleScalarSqrt(s_B);
        if (((s * (s * (s * (s * (s * 6.536332E-9 + -1.120083E-6) +
                             0.0001001685) +
                        -0.00909529) +
                   0.06793952) +
              999.842594) +
             m * (((s * (s * (s * (s * 5.3875E-9 + -8.2467E-7) + 7.6438E-5) +
                         -0.0040899) +
                    0.824493) +
                   s1o2 * (s * (s * -1.6546E-6 + 0.00010227) + -0.00572466)) +
                  m * 0.00048314)) /
                    (1.0 -
                     z_tmp /
                         (((s * (s * (s * (s * -0.0004190253 + 0.09648704) +
                                      -17.06103) +
                                 1444.304) +
                            196593.3) +
                           m * ((s * (s * (s * -0.0005084188 + 0.06283263) +
                                      -3.101089) +
                                 528.4855) +
                                s1o2 * (s * (s * -0.004619924 + 0.09085835) +
                                        3.88664))) +
                          z_tmp *
                              (((s * (s * (s * 1.956415E-6 + -0.0002984642) +
                                      0.02212276) +
                                 3.186519) +
                                m * ((s * (s * 2.059331E-7 + -0.0001847318) +
                                      0.006704388) +
                                     s1o2 * 0.0001480266)) +
                               z_tmp *
                                   ((s * (s * 1.39468E-8 + -1.202016E-6) +
                                     2.102898E-5) +
                                    m * (s * (s * 6.207323E-11 + 6.128773E-9) +
                                         -2.040237E-7))))) -
                ((t_B *
                      (t_B * (t_B * (t_B * (t_B * 6.536332E-9 + -1.120083E-6) +
                                     0.0001001685) +
                              -0.00909529) +
                       0.06793952) +
                  999.842594) +
                 s_B * (((t_B * (t_B * (t_B * (t_B * 5.3875E-9 + -8.2467E-7) +
                                        7.6438E-5) +
                                 -0.0040899) +
                          0.824493) +
                         b_s1o2 * (t_B * (t_B * -1.6546E-6 + 0.00010227) +
                                   -0.00572466)) +
                        s_B * 0.00048314)) /
                    (1.0 -
                     z_tmp /
                         (((t_B * (t_B * (t_B * (t_B * -0.0004190253 +
                                                 0.09648704) +
                                          -17.06103) +
                                   1444.304) +
                            196593.3) +
                           s_B * ((t_B * (t_B * (t_B * -0.0005084188 +
                                                 0.06283263) +
                                          -3.101089) +
                                   528.4855) +
                                  b_s1o2 *
                                      (t_B * (t_B * -0.004619924 + 0.09085835) +
                                       3.88664))) +
                          z_tmp * (((t_B * (t_B * (t_B * 1.956415E-6 +
                                                   -0.0002984642) +
                                            0.02212276) +
                                     3.186519) +
                                    s_B * ((t_B * (t_B * 2.059331E-7 +
                                                   -0.0001847318) +
                                            0.006704388) +
                                           b_s1o2 * 0.0001480266)) +
                                   z_tmp * ((t_B * (t_B * 1.39468E-8 +
                                                    -1.202016E-6) +
                                             2.102898E-5) +
                                            s_B * (t_B * (t_B * 6.207323E-11 +
                                                          6.128773E-9) +
                                                   -2.040237E-7))))) ==
            0.0) {
          lb = r;
          dp = r;
        } else {
          fa = muDoubleScalarMin(muDoubleScalarMax(fa, q), r);
          /*  bounds are given */
          dxp = (r - fa) / 50.0;
          dxm = (fa - q) / 50.0;
          /*  Set a = x, except when x is so close to A that machine roundoff
           * makes dxm identically 0, */
          /*  which would lead to an infinite loop below.  In this case, set a =
           * A. */
          if (dxm == 0.0) {
            lb = q;
          } else {
            lb = fa;
          }
          /*  Evaluate difference between (a) eos at location on cast A with */
          /*  properties (S_A, T_A, P_A) where the pressure or depth is p_avg +
           * dp/2, */
          /*  and (b) eos at location on cast B with properties (S_B, T_B, P_B)
           * where */
          /*  the pressure or depth is p_avg - dp/2. Here, eos is always
           * evaluated at */
          /*  the mid-pressure or mid-depth, p_avg. */
          ppc_val2(P_A_data, P_A_size, Sppc_A_data, Sppc_A_size, Tppc_A_data,
                   p_avg + 0.5 * lb, &m, &s);
          ppc_val2(P_B_data, P_B_size, Sppc_B_data, Sppc_B_size, Tppc_B_data,
                   p_avg - 0.5 * lb, &s_B, &t_B);
          /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
          /*  */
          /*  */
          /*  rho = eoscg_densjmd95_bsq(s,t,z)                             [ kg
           * / m^3 ] */
          /*  computes the Boussinesq JMD95 in-situ density, given practical
           * salinity */
          /*  s, potential temperature t, and depth z [m, positive and
           * increasing */
          /*  down].  The depth z is converted into hydrostatic pressure for a
           * given */
          /*  gravitational acceleration and Boussinesq reference density, which
           * are */
          /*  hard-coded into this function (edit variables grav and rhob as
           * needed). */
          /*  */
          /*  */
          /*  This function is derived from densjmd95.m, documented below. Input
           * checks */
          /*  and expansion of variables have been removed; instead automatic
           * expansion */
          /*  (requiring MATLAB 2016b or later) is used. The calculation has
           * also been */
          /*  streamlined by pre-allocating arrays for coeffcients, and
           * modifying the */
          /*  coefficients such that pressure need not be internally converted
           * from */
          /*  dbar to bar. This function is compatible with MATLAB's codegen. */
          /*  */
          /*  The units of s and t are as in densjmd95, documented below. */
          /*   */
          /*  The inputs' sizes are more general: they may be arrays of any
           * dimension */
          /*  and any size, so long as their sizes match each other excepting
           * that any */
          /*  input can have any dimension as a singleton. For example, s and t
           * can be */
          /*  3D arrays of size [nz, nx, ny], while p can be a vector of size
           * [nz,1]. */
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
          /*     Density of Sea Water using Jackett and McDougall 1995 (JAOT 12)
           */
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
          /*  Convert depth [m] to the hydrostatic pressure [dbar] implied by
           * the */
          /*  following two hard-coded parameters (edit as needed): */
          /*  gravitational acceleration [m /s^2] */
          /*  Boussinesq reference density [kg / m^3] */
          /*  Pascal to dbar conversion [dbar / Pa] */
          /*  depth to pressure conversion [dbar / m] */
          /*  Henceforth z is actually pressure [dbar] */
          /*  coefficients nonlinear equation of state in pressure coordinates
           * for */
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
          s1o2 = muDoubleScalarSqrt(m);
          /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
          /*  */
          /*  */
          /*  rho = eoscg_densjmd95_bsq(s,t,z)                             [ kg
           * / m^3 ] */
          /*  computes the Boussinesq JMD95 in-situ density, given practical
           * salinity */
          /*  s, potential temperature t, and depth z [m, positive and
           * increasing */
          /*  down].  The depth z is converted into hydrostatic pressure for a
           * given */
          /*  gravitational acceleration and Boussinesq reference density, which
           * are */
          /*  hard-coded into this function (edit variables grav and rhob as
           * needed). */
          /*  */
          /*  */
          /*  This function is derived from densjmd95.m, documented below. Input
           * checks */
          /*  and expansion of variables have been removed; instead automatic
           * expansion */
          /*  (requiring MATLAB 2016b or later) is used. The calculation has
           * also been */
          /*  streamlined by pre-allocating arrays for coeffcients, and
           * modifying the */
          /*  coefficients such that pressure need not be internally converted
           * from */
          /*  dbar to bar. This function is compatible with MATLAB's codegen. */
          /*  */
          /*  The units of s and t are as in densjmd95, documented below. */
          /*   */
          /*  The inputs' sizes are more general: they may be arrays of any
           * dimension */
          /*  and any size, so long as their sizes match each other excepting
           * that any */
          /*  input can have any dimension as a singleton. For example, s and t
           * can be */
          /*  3D arrays of size [nz, nx, ny], while p can be a vector of size
           * [nz,1]. */
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
          /*     Density of Sea Water using Jackett and McDougall 1995 (JAOT 12)
           */
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
          /*  Convert depth [m] to the hydrostatic pressure [dbar] implied by
           * the */
          /*  following two hard-coded parameters (edit as needed): */
          /*  gravitational acceleration [m /s^2] */
          /*  Boussinesq reference density [kg / m^3] */
          /*  Pascal to dbar conversion [dbar / Pa] */
          /*  depth to pressure conversion [dbar / m] */
          /*  Henceforth z is actually pressure [dbar] */
          /*  coefficients nonlinear equation of state in pressure coordinates
           * for */
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
          b_s1o2 = muDoubleScalarSqrt(s_B);
          fapos =
              (((s * (s * (s * (s * (s * 6.536332E-9 + -1.120083E-6) +
                                0.0001001685) +
                           -0.00909529) +
                      0.06793952) +
                 999.842594) +
                m * (((s * (s * (s * (s * 5.3875E-9 + -8.2467E-7) + 7.6438E-5) +
                            -0.0040899) +
                       0.824493) +
                      s1o2 *
                          (s * (s * -1.6546E-6 + 0.00010227) + -0.00572466)) +
                     m * 0.00048314)) /
                       (1.0 -
                        z_tmp /
                            (((s * (s * (s * (s * -0.0004190253 + 0.09648704) +
                                         -17.06103) +
                                    1444.304) +
                               196593.3) +
                              m * ((s * (s * (s * -0.0005084188 + 0.06283263) +
                                         -3.101089) +
                                    528.4855) +
                                   s1o2 * (s * (s * -0.004619924 + 0.09085835) +
                                           3.88664))) +
                             z_tmp *
                                 (((s * (s * (s * 1.956415E-6 + -0.0002984642) +
                                         0.02212276) +
                                    3.186519) +
                                   m * ((s * (s * 2.059331E-7 + -0.0001847318) +
                                         0.006704388) +
                                        s1o2 * 0.0001480266)) +
                                  z_tmp *
                                      ((s * (s * 1.39468E-8 + -1.202016E-6) +
                                        2.102898E-5) +
                                       m * (s * (s * 6.207323E-11 +
                                                 6.128773E-9) +
                                            -2.040237E-7))))) -
                   ((t_B * (t_B * (t_B * (t_B * (t_B * 6.536332E-9 +
                                                 -1.120083E-6) +
                                          0.0001001685) +
                                   -0.00909529) +
                            0.06793952) +
                     999.842594) +
                    s_B *
                        (((t_B * (t_B * (t_B * (t_B * 5.3875E-9 + -8.2467E-7) +
                                         7.6438E-5) +
                                  -0.0040899) +
                           0.824493) +
                          b_s1o2 * (t_B * (t_B * -1.6546E-6 + 0.00010227) +
                                    -0.00572466)) +
                         s_B * 0.00048314)) /
                       (1.0 -
                        z_tmp /
                            (((t_B * (t_B * (t_B * (t_B * -0.0004190253 +
                                                    0.09648704) +
                                             -17.06103) +
                                      1444.304) +
                               196593.3) +
                              s_B * ((t_B * (t_B * (t_B * -0.0005084188 +
                                                    0.06283263) +
                                             -3.101089) +
                                      528.4855) +
                                     b_s1o2 * (t_B * (t_B * -0.004619924 +
                                                      0.09085835) +
                                               3.88664))) +
                             z_tmp *
                                 (((t_B * (t_B * (t_B * 1.956415E-6 +
                                                  -0.0002984642) +
                                           0.02212276) +
                                    3.186519) +
                                   s_B * ((t_B * (t_B * 2.059331E-7 +
                                                  -0.0001847318) +
                                           0.006704388) +
                                          b_s1o2 * 0.0001480266)) +
                                  z_tmp * ((t_B * (t_B * 1.39468E-8 +
                                                   -1.202016E-6) +
                                            2.102898E-5) +
                                           s_B * (t_B * (t_B * 6.207323E-11 +
                                                         6.128773E-9) +
                                                  -2.040237E-7))))) >
               0.0);
          /*  Similarly, set b = x, except for machine precision problems. */
          if (dxp == 0.0) {
            dp = r;
            /*  Evaluate difference between (a) eos at location on cast A with
             */
            /*  properties (S_A, T_A, P_A) where the pressure or depth is p_avg
             * + dp/2, */
            /*  and (b) eos at location on cast B with properties (S_B, T_B,
             * P_B) where */
            /*  the pressure or depth is p_avg - dp/2. Here, eos is always
             * evaluated at */
            /*  the mid-pressure or mid-depth, p_avg. */
            ppc_val2(P_A_data, P_A_size, Sppc_A_data, Sppc_A_size, Tppc_A_data,
                     p_avg + 0.5 * r, &m, &s);
            ppc_val2(P_B_data, P_B_size, Sppc_B_data, Sppc_B_size, Tppc_B_data,
                     p_avg - 0.5 * r, &s_B, &t_B);
            /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
            /*  */
            /*  */
            /*  rho = eoscg_densjmd95_bsq(s,t,z)                             [
             * kg / m^3 ] */
            /*  computes the Boussinesq JMD95 in-situ density, given practical
             * salinity */
            /*  s, potential temperature t, and depth z [m, positive and
             * increasing */
            /*  down].  The depth z is converted into hydrostatic pressure for a
             * given */
            /*  gravitational acceleration and Boussinesq reference density,
             * which are */
            /*  hard-coded into this function (edit variables grav and rhob as
             * needed). */
            /*  */
            /*  */
            /*  This function is derived from densjmd95.m, documented below.
             * Input checks */
            /*  and expansion of variables have been removed; instead automatic
             * expansion */
            /*  (requiring MATLAB 2016b or later) is used. The calculation has
             * also been */
            /*  streamlined by pre-allocating arrays for coeffcients, and
             * modifying the */
            /*  coefficients such that pressure need not be internally converted
             * from */
            /*  dbar to bar. This function is compatible with MATLAB's codegen.
             */
            /*  */
            /*  The units of s and t are as in densjmd95, documented below. */
            /*   */
            /*  The inputs' sizes are more general: they may be arrays of any
             * dimension */
            /*  and any size, so long as their sizes match each other excepting
             * that any */
            /*  input can have any dimension as a singleton. For example, s and
             * t can be */
            /*  3D arrays of size [nz, nx, ny], while p can be a vector of size
             * [nz,1]. */
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
            /*  Convert depth [m] to the hydrostatic pressure [dbar] implied by
             * the */
            /*  following two hard-coded parameters (edit as needed): */
            /*  gravitational acceleration [m /s^2] */
            /*  Boussinesq reference density [kg / m^3] */
            /*  Pascal to dbar conversion [dbar / Pa] */
            /*  depth to pressure conversion [dbar / m] */
            /*  Henceforth z is actually pressure [dbar] */
            /*  coefficients nonlinear equation of state in pressure coordinates
             * for */
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
            s1o2 = muDoubleScalarSqrt(m);
            /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
            /*  */
            /*  */
            /*  rho = eoscg_densjmd95_bsq(s,t,z)                             [
             * kg / m^3 ] */
            /*  computes the Boussinesq JMD95 in-situ density, given practical
             * salinity */
            /*  s, potential temperature t, and depth z [m, positive and
             * increasing */
            /*  down].  The depth z is converted into hydrostatic pressure for a
             * given */
            /*  gravitational acceleration and Boussinesq reference density,
             * which are */
            /*  hard-coded into this function (edit variables grav and rhob as
             * needed). */
            /*  */
            /*  */
            /*  This function is derived from densjmd95.m, documented below.
             * Input checks */
            /*  and expansion of variables have been removed; instead automatic
             * expansion */
            /*  (requiring MATLAB 2016b or later) is used. The calculation has
             * also been */
            /*  streamlined by pre-allocating arrays for coeffcients, and
             * modifying the */
            /*  coefficients such that pressure need not be internally converted
             * from */
            /*  dbar to bar. This function is compatible with MATLAB's codegen.
             */
            /*  */
            /*  The units of s and t are as in densjmd95, documented below. */
            /*   */
            /*  The inputs' sizes are more general: they may be arrays of any
             * dimension */
            /*  and any size, so long as their sizes match each other excepting
             * that any */
            /*  input can have any dimension as a singleton. For example, s and
             * t can be */
            /*  3D arrays of size [nz, nx, ny], while p can be a vector of size
             * [nz,1]. */
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
            /*  Convert depth [m] to the hydrostatic pressure [dbar] implied by
             * the */
            /*  following two hard-coded parameters (edit as needed): */
            /*  gravitational acceleration [m /s^2] */
            /*  Boussinesq reference density [kg / m^3] */
            /*  Pascal to dbar conversion [dbar / Pa] */
            /*  depth to pressure conversion [dbar / m] */
            /*  Henceforth z is actually pressure [dbar] */
            /*  coefficients nonlinear equation of state in pressure coordinates
             * for */
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
            b_s1o2 = muDoubleScalarSqrt(s_B);
            fbpos =
                (((s * (s * (s * (s * (s * 6.536332E-9 + -1.120083E-6) +
                                  0.0001001685) +
                             -0.00909529) +
                        0.06793952) +
                   999.842594) +
                  m * (((s * (s * (s * (s * 5.3875E-9 + -8.2467E-7) +
                                   7.6438E-5) +
                              -0.0040899) +
                         0.824493) +
                        s1o2 *
                            (s * (s * -1.6546E-6 + 0.00010227) + -0.00572466)) +
                       m * 0.00048314)) /
                         (1.0 -
                          z_tmp /
                              (((s * (s * (s * (s * -0.0004190253 +
                                                0.09648704) +
                                           -17.06103) +
                                      1444.304) +
                                 196593.3) +
                                m * ((s * (s * (s * -0.0005084188 +
                                                0.06283263) +
                                           -3.101089) +
                                      528.4855) +
                                     s1o2 *
                                         (s * (s * -0.004619924 + 0.09085835) +
                                          3.88664))) +
                               z_tmp * (((s * (s * (s * 1.956415E-6 +
                                                    -0.0002984642) +
                                               0.02212276) +
                                          3.186519) +
                                         m * ((s * (s * 2.059331E-7 +
                                                    -0.0001847318) +
                                               0.006704388) +
                                              s1o2 * 0.0001480266)) +
                                        z_tmp * ((s * (s * 1.39468E-8 +
                                                       -1.202016E-6) +
                                                  2.102898E-5) +
                                                 m * (s * (s * 6.207323E-11 +
                                                           6.128773E-9) +
                                                      -2.040237E-7))))) -
                     ((t_B * (t_B * (t_B * (t_B * (t_B * 6.536332E-9 +
                                                   -1.120083E-6) +
                                            0.0001001685) +
                                     -0.00909529) +
                              0.06793952) +
                       999.842594) +
                      s_B * (((t_B * (t_B * (t_B * (t_B * 5.3875E-9 +
                                                    -8.2467E-7) +
                                             7.6438E-5) +
                                      -0.0040899) +
                               0.824493) +
                              b_s1o2 * (t_B * (t_B * -1.6546E-6 + 0.00010227) +
                                        -0.00572466)) +
                             s_B * 0.00048314)) /
                         (1.0 -
                          z_tmp /
                              (((t_B * (t_B * (t_B * (t_B * -0.0004190253 +
                                                      0.09648704) +
                                               -17.06103) +
                                        1444.304) +
                                 196593.3) +
                                s_B * ((t_B * (t_B * (t_B * -0.0005084188 +
                                                      0.06283263) +
                                               -3.101089) +
                                        528.4855) +
                                       b_s1o2 * (t_B * (t_B * -0.004619924 +
                                                        0.09085835) +
                                                 3.88664))) +
                               z_tmp *
                                   (((t_B * (t_B * (t_B * 1.956415E-6 +
                                                    -0.0002984642) +
                                             0.02212276) +
                                      3.186519) +
                                     s_B * ((t_B * (t_B * 2.059331E-7 +
                                                    -0.0001847318) +
                                             0.006704388) +
                                            b_s1o2 * 0.0001480266)) +
                                    z_tmp * ((t_B * (t_B * 1.39468E-8 +
                                                     -1.202016E-6) +
                                              2.102898E-5) +
                                             s_B * (t_B * (t_B * 6.207323E-11 +
                                                           6.128773E-9) +
                                                    -2.040237E-7))))) >
                 0.0);
          } else {
            dp = fa;
            if (dxm == 0.0) {
              fbpos = fapos;
              /*  since a = b = x */
            } else {
              /*  Evaluate difference between (a) eos at location on cast A with
               */
              /*  properties (S_A, T_A, P_A) where the pressure or depth is
               * p_avg + dp/2, */
              /*  and (b) eos at location on cast B with properties (S_B, T_B,
               * P_B) where */
              /*  the pressure or depth is p_avg - dp/2. Here, eos is always
               * evaluated at */
              /*  the mid-pressure or mid-depth, p_avg. */
              ppc_val2(P_A_data, P_A_size, Sppc_A_data, Sppc_A_size,
                       Tppc_A_data, p_avg + 0.5 * fa, &m, &s);
              ppc_val2(P_B_data, P_B_size, Sppc_B_data, Sppc_B_size,
                       Tppc_B_data, p_avg - 0.5 * fa, &s_B, &t_B);
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
              s1o2 = muDoubleScalarSqrt(m);
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
              b_s1o2 = muDoubleScalarSqrt(s_B);
              fbpos =
                  (((s * (s * (s * (s * (s * 6.536332E-9 + -1.120083E-6) +
                                    0.0001001685) +
                               -0.00909529) +
                          0.06793952) +
                     999.842594) +
                    m * (((s * (s * (s * (s * 5.3875E-9 + -8.2467E-7) +
                                     7.6438E-5) +
                                -0.0040899) +
                           0.824493) +
                          s1o2 * (s * (s * -1.6546E-6 + 0.00010227) +
                                  -0.00572466)) +
                         m * 0.00048314)) /
                           (1.0 -
                            z_tmp /
                                (((s * (s * (s * (s * -0.0004190253 +
                                                  0.09648704) +
                                             -17.06103) +
                                        1444.304) +
                                   196593.3) +
                                  m * ((s * (s * (s * -0.0005084188 +
                                                  0.06283263) +
                                             -3.101089) +
                                        528.4855) +
                                       s1o2 * (s * (s * -0.004619924 +
                                                    0.09085835) +
                                               3.88664))) +
                                 z_tmp * (((s * (s * (s * 1.956415E-6 +
                                                      -0.0002984642) +
                                                 0.02212276) +
                                            3.186519) +
                                           m * ((s * (s * 2.059331E-7 +
                                                      -0.0001847318) +
                                                 0.006704388) +
                                                s1o2 * 0.0001480266)) +
                                          z_tmp * ((s * (s * 1.39468E-8 +
                                                         -1.202016E-6) +
                                                    2.102898E-5) +
                                                   m * (s * (s * 6.207323E-11 +
                                                             6.128773E-9) +
                                                        -2.040237E-7))))) -
                       ((t_B * (t_B * (t_B * (t_B * (t_B * 6.536332E-9 +
                                                     -1.120083E-6) +
                                              0.0001001685) +
                                       -0.00909529) +
                                0.06793952) +
                         999.842594) +
                        s_B *
                            (((t_B * (t_B * (t_B * (t_B * 5.3875E-9 +
                                                    -8.2467E-7) +
                                             7.6438E-5) +
                                      -0.0040899) +
                               0.824493) +
                              b_s1o2 * (t_B * (t_B * -1.6546E-6 + 0.00010227) +
                                        -0.00572466)) +
                             s_B * 0.00048314)) /
                           (1.0 -
                            z_tmp /
                                (((t_B * (t_B * (t_B * (t_B * -0.0004190253 +
                                                        0.09648704) +
                                                 -17.06103) +
                                          1444.304) +
                                   196593.3) +
                                  s_B * ((t_B * (t_B * (t_B * -0.0005084188 +
                                                        0.06283263) +
                                                 -3.101089) +
                                          528.4855) +
                                         b_s1o2 * (t_B * (t_B * -0.004619924 +
                                                          0.09085835) +
                                                   3.88664))) +
                                 z_tmp *
                                     (((t_B * (t_B * (t_B * 1.956415E-6 +
                                                      -0.0002984642) +
                                               0.02212276) +
                                        3.186519) +
                                       s_B * ((t_B * (t_B * 2.059331E-7 +
                                                      -0.0001847318) +
                                               0.006704388) +
                                              b_s1o2 * 0.0001480266)) +
                                      z_tmp *
                                          ((t_B * (t_B * 1.39468E-8 +
                                                   -1.202016E-6) +
                                            2.102898E-5) +
                                           s_B * (t_B * (t_B * 6.207323E-11 +
                                                         6.128773E-9) +
                                                  -2.040237E-7))))) >
                   0.0);
            }
          }
          do {
            exitg1 = 0;
            guard1 = false;
            if (lb > q) {
              /*  Move a left, and test for a sign change */
              dxm *= 1.4142135623730949;
              lb = muDoubleScalarMax(fa - dxm, q);
              /*  Evaluate difference between (a) eos at location on cast A with
               */
              /*  properties (S_A, T_A, P_A) where the pressure or depth is
               * p_avg + dp/2, */
              /*  and (b) eos at location on cast B with properties (S_B, T_B,
               * P_B) where */
              /*  the pressure or depth is p_avg - dp/2. Here, eos is always
               * evaluated at */
              /*  the mid-pressure or mid-depth, p_avg. */
              ppc_val2(P_A_data, P_A_size, Sppc_A_data, Sppc_A_size,
                       Tppc_A_data, p_avg + 0.5 * lb, &m, &s);
              ppc_val2(P_B_data, P_B_size, Sppc_B_data, Sppc_B_size,
                       Tppc_B_data, p_avg - 0.5 * lb, &s_B, &t_B);
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
              s1o2 = muDoubleScalarSqrt(m);
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
              b_s1o2 = muDoubleScalarSqrt(s_B);
              fapos =
                  (((s * (s * (s * (s * (s * 6.536332E-9 + -1.120083E-6) +
                                    0.0001001685) +
                               -0.00909529) +
                          0.06793952) +
                     999.842594) +
                    m * (((s * (s * (s * (s * 5.3875E-9 + -8.2467E-7) +
                                     7.6438E-5) +
                                -0.0040899) +
                           0.824493) +
                          s1o2 * (s * (s * -1.6546E-6 + 0.00010227) +
                                  -0.00572466)) +
                         m * 0.00048314)) /
                           (1.0 -
                            z_tmp /
                                (((s * (s * (s * (s * -0.0004190253 +
                                                  0.09648704) +
                                             -17.06103) +
                                        1444.304) +
                                   196593.3) +
                                  m * ((s * (s * (s * -0.0005084188 +
                                                  0.06283263) +
                                             -3.101089) +
                                        528.4855) +
                                       s1o2 * (s * (s * -0.004619924 +
                                                    0.09085835) +
                                               3.88664))) +
                                 z_tmp * (((s * (s * (s * 1.956415E-6 +
                                                      -0.0002984642) +
                                                 0.02212276) +
                                            3.186519) +
                                           m * ((s * (s * 2.059331E-7 +
                                                      -0.0001847318) +
                                                 0.006704388) +
                                                s1o2 * 0.0001480266)) +
                                          z_tmp * ((s * (s * 1.39468E-8 +
                                                         -1.202016E-6) +
                                                    2.102898E-5) +
                                                   m * (s * (s * 6.207323E-11 +
                                                             6.128773E-9) +
                                                        -2.040237E-7))))) -
                       ((t_B * (t_B * (t_B * (t_B * (t_B * 6.536332E-9 +
                                                     -1.120083E-6) +
                                              0.0001001685) +
                                       -0.00909529) +
                                0.06793952) +
                         999.842594) +
                        s_B *
                            (((t_B * (t_B * (t_B * (t_B * 5.3875E-9 +
                                                    -8.2467E-7) +
                                             7.6438E-5) +
                                      -0.0040899) +
                               0.824493) +
                              b_s1o2 * (t_B * (t_B * -1.6546E-6 + 0.00010227) +
                                        -0.00572466)) +
                             s_B * 0.00048314)) /
                           (1.0 -
                            z_tmp /
                                (((t_B * (t_B * (t_B * (t_B * -0.0004190253 +
                                                        0.09648704) +
                                                 -17.06103) +
                                          1444.304) +
                                   196593.3) +
                                  s_B * ((t_B * (t_B * (t_B * -0.0005084188 +
                                                        0.06283263) +
                                                 -3.101089) +
                                          528.4855) +
                                         b_s1o2 * (t_B * (t_B * -0.004619924 +
                                                          0.09085835) +
                                                   3.88664))) +
                                 z_tmp *
                                     (((t_B * (t_B * (t_B * 1.956415E-6 +
                                                      -0.0002984642) +
                                               0.02212276) +
                                        3.186519) +
                                       s_B * ((t_B * (t_B * 2.059331E-7 +
                                                      -0.0001847318) +
                                               0.006704388) +
                                              b_s1o2 * 0.0001480266)) +
                                      z_tmp *
                                          ((t_B * (t_B * 1.39468E-8 +
                                                   -1.202016E-6) +
                                            2.102898E-5) +
                                           s_B * (t_B * (t_B * 6.207323E-11 +
                                                         6.128773E-9) +
                                                  -2.040237E-7))))) >
                   0.0);
              if (fapos != fbpos) {
                /*  fa and fb have different signs */
                exitg1 = 1;
              } else {
                guard1 = true;
              }
            } else if (dp == r) {
              /*  also a == A, so cannot expand anymore */
              if (fapos == fbpos) {
                /*  no sign change found */
                lb = rtNaN;
                dp = rtNaN;
              } else {
                /*  one last test for sign change */
              }
              exitg1 = 1;
            } else {
              guard1 = true;
            }
            if (guard1) {
              if (dp < r) {
                /*  Move b right, and test for a sign change */
                dxp *= 1.4142135623730949;
                dp = muDoubleScalarMin(fa + dxp, r);
                /*  Evaluate difference between (a) eos at location on cast A
                 * with */
                /*  properties (S_A, T_A, P_A) where the pressure or depth is
                 * p_avg + dp/2, */
                /*  and (b) eos at location on cast B with properties (S_B, T_B,
                 * P_B) where */
                /*  the pressure or depth is p_avg - dp/2. Here, eos is always
                 * evaluated at */
                /*  the mid-pressure or mid-depth, p_avg. */
                ppc_val2(P_A_data, P_A_size, Sppc_A_data, Sppc_A_size,
                         Tppc_A_data, p_avg + 0.5 * dp, &m, &s);
                ppc_val2(P_B_data, P_B_size, Sppc_B_data, Sppc_B_size,
                         Tppc_B_data, p_avg - 0.5 * dp, &s_B, &t_B);
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
                s1o2 = muDoubleScalarSqrt(m);
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
                b_s1o2 = muDoubleScalarSqrt(s_B);
                fbpos =
                    (((s * (s * (s * (s * (s * 6.536332E-9 + -1.120083E-6) +
                                      0.0001001685) +
                                 -0.00909529) +
                            0.06793952) +
                       999.842594) +
                      m * (((s * (s * (s * (s * 5.3875E-9 + -8.2467E-7) +
                                       7.6438E-5) +
                                  -0.0040899) +
                             0.824493) +
                            s1o2 * (s * (s * -1.6546E-6 + 0.00010227) +
                                    -0.00572466)) +
                           m * 0.00048314)) /
                             (1.0 -
                              z_tmp /
                                  (((s * (s * (s * (s * -0.0004190253 +
                                                    0.09648704) +
                                               -17.06103) +
                                          1444.304) +
                                     196593.3) +
                                    m * ((s * (s * (s * -0.0005084188 +
                                                    0.06283263) +
                                               -3.101089) +
                                          528.4855) +
                                         s1o2 * (s * (s * -0.004619924 +
                                                      0.09085835) +
                                                 3.88664))) +
                                   z_tmp *
                                       (((s * (s * (s * 1.956415E-6 +
                                                    -0.0002984642) +
                                               0.02212276) +
                                          3.186519) +
                                         m * ((s * (s * 2.059331E-7 +
                                                    -0.0001847318) +
                                               0.006704388) +
                                              s1o2 * 0.0001480266)) +
                                        z_tmp * ((s * (s * 1.39468E-8 +
                                                       -1.202016E-6) +
                                                  2.102898E-5) +
                                                 m * (s * (s * 6.207323E-11 +
                                                           6.128773E-9) +
                                                      -2.040237E-7))))) -
                         ((t_B * (t_B * (t_B * (t_B * (t_B * 6.536332E-9 +
                                                       -1.120083E-6) +
                                                0.0001001685) +
                                         -0.00909529) +
                                  0.06793952) +
                           999.842594) +
                          s_B * (((t_B * (t_B * (t_B * (t_B * 5.3875E-9 +
                                                        -8.2467E-7) +
                                                 7.6438E-5) +
                                          -0.0040899) +
                                   0.824493) +
                                  b_s1o2 *
                                      (t_B * (t_B * -1.6546E-6 + 0.00010227) +
                                       -0.00572466)) +
                                 s_B * 0.00048314)) /
                             (1.0 -
                              z_tmp /
                                  (((t_B * (t_B * (t_B * (t_B * -0.0004190253 +
                                                          0.09648704) +
                                                   -17.06103) +
                                            1444.304) +
                                     196593.3) +
                                    s_B * ((t_B * (t_B * (t_B * -0.0005084188 +
                                                          0.06283263) +
                                                   -3.101089) +
                                            528.4855) +
                                           b_s1o2 * (t_B * (t_B * -0.004619924 +
                                                            0.09085835) +
                                                     3.88664))) +
                                   z_tmp *
                                       (((t_B * (t_B * (t_B * 1.956415E-6 +
                                                        -0.0002984642) +
                                                 0.02212276) +
                                          3.186519) +
                                         s_B * ((t_B * (t_B * 2.059331E-7 +
                                                        -0.0001847318) +
                                                 0.006704388) +
                                                b_s1o2 * 0.0001480266)) +
                                        z_tmp *
                                            ((t_B * (t_B * 1.39468E-8 +
                                                     -1.202016E-6) +
                                              2.102898E-5) +
                                             s_B * (t_B * (t_B * 6.207323E-11 +
                                                           6.128773E-9) +
                                                    -2.040237E-7))))) >
                     0.0);
                if (fapos != fbpos) {
                  /*  fa and fb have different signs */
                  exitg1 = 1;
                }
              } else if (lb == q) {
                /*  also b == B, so cannot expand anymore */
                if (fapos == fbpos) {
                  /*  no sign change found */
                  lb = rtNaN;
                  dp = rtNaN;
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
      /*  A sign change was discovered, so a root exists in the interval. */
      /*  Solve the nonlinear root-finding problem using Brent's method */
      /* FZERO_BRENT  Find a root of a univariate function within a given
       * interval */
      /*              using Brent's method */
      /*  */
      /*  x = fzero_brent(f,a,b,t)  */
      /*  finds x in the interval [a,b] satisfying |x-y| <= t/2 where f(y) = 0.
       */
      /*  f(a) and f(b) must have opposite signs. */
      /*  */
      /*  ... = fzero_brent(f,lb,ub,t,x,...) passes additional inputs directly
       */
      /*  to f. */
      /*  */
      /*  This function is compatible with MATLAB's code generation -- so long f
       * is */
      /*  similarly compatible. Many root-finding problems can be solved with */
      /*  fzero_brent by writing another function which calls fzero_brent inside
       * a */
      /*  for loop, and this can be made fast by using code generation on that
       */
      /*  wrapper function. Note that f will only be called with a scalar input
       * as */
      /*  its first argument; codegen knows this, and might strip out
       * unnecessary */
      /*  code from the function definition underlying f. */
      /*  --- Example: */
      /*  Simple bisection between bounds of -0.5 and 1.5 would fail to find the
       */
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
      /*     assuming that F is continuous implies the existence of at least */
      /*     one value C between A and B for which F(C) = 0. */
      /*  */
      /*     The location of the zero is determined to within an accuracy */
      /*     of 6 * EPS * abs ( C ) + 2 * T, where EPS is the machine epsilon.
       */
      /*  */
      /*     Thanks to Thomas Secretin for pointing out a transcription error in
       * the */
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
      /*     Input, real A, B, the endpoints of the change of sign interval. */
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
      /*  Evaluate difference between (a) eos at location on cast A with */
      /*  properties (S_A, T_A, P_A) where the pressure or depth is p_avg +
       * dp/2, */
      /*  and (b) eos at location on cast B with properties (S_B, T_B, P_B)
       * where */
      /*  the pressure or depth is p_avg - dp/2. Here, eos is always evaluated
       * at */
      /*  the mid-pressure or mid-depth, p_avg. */
      ppc_val2(P_A_data, P_A_size, Sppc_A_data, Sppc_A_size, Tppc_A_data,
               p_avg + 0.5 * lb, &m, &s);
      ppc_val2(P_B_data, P_B_size, Sppc_B_data, Sppc_B_size, Tppc_B_data,
               p_avg - 0.5 * lb, &s_B, &t_B);
      /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
      /*  */
      /*  */
      /*  rho = eoscg_densjmd95_bsq(s,t,z)                             [ kg /
       * m^3 ] */
      /*  computes the Boussinesq JMD95 in-situ density, given practical
       * salinity */
      /*  s, potential temperature t, and depth z [m, positive and increasing */
      /*  down].  The depth z is converted into hydrostatic pressure for a given
       */
      /*  gravitational acceleration and Boussinesq reference density, which are
       */
      /*  hard-coded into this function (edit variables grav and rhob as
       * needed). */
      /*  */
      /*  */
      /*  This function is derived from densjmd95.m, documented below. Input
       * checks */
      /*  and expansion of variables have been removed; instead automatic
       * expansion */
      /*  (requiring MATLAB 2016b or later) is used. The calculation has also
       * been */
      /*  streamlined by pre-allocating arrays for coeffcients, and modifying
       * the */
      /*  coefficients such that pressure need not be internally converted from
       */
      /*  dbar to bar. This function is compatible with MATLAB's codegen. */
      /*  */
      /*  The units of s and t are as in densjmd95, documented below. */
      /*   */
      /*  The inputs' sizes are more general: they may be arrays of any
       * dimension */
      /*  and any size, so long as their sizes match each other excepting that
       * any */
      /*  input can have any dimension as a singleton. For example, s and t can
       * be */
      /*  3D arrays of size [nz, nx, ny], while p can be a vector of size
       * [nz,1]. */
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
      z_tmp = p_avg * 1.015335;
      /*  Henceforth z is actually pressure [dbar] */
      /*  coefficients nonlinear equation of state in pressure coordinates for
       */
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
      s1o2 = muDoubleScalarSqrt(m);
      /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
      /*  */
      /*  */
      /*  rho = eoscg_densjmd95_bsq(s,t,z)                             [ kg /
       * m^3 ] */
      /*  computes the Boussinesq JMD95 in-situ density, given practical
       * salinity */
      /*  s, potential temperature t, and depth z [m, positive and increasing */
      /*  down].  The depth z is converted into hydrostatic pressure for a given
       */
      /*  gravitational acceleration and Boussinesq reference density, which are
       */
      /*  hard-coded into this function (edit variables grav and rhob as
       * needed). */
      /*  */
      /*  */
      /*  This function is derived from densjmd95.m, documented below. Input
       * checks */
      /*  and expansion of variables have been removed; instead automatic
       * expansion */
      /*  (requiring MATLAB 2016b or later) is used. The calculation has also
       * been */
      /*  streamlined by pre-allocating arrays for coeffcients, and modifying
       * the */
      /*  coefficients such that pressure need not be internally converted from
       */
      /*  dbar to bar. This function is compatible with MATLAB's codegen. */
      /*  */
      /*  The units of s and t are as in densjmd95, documented below. */
      /*   */
      /*  The inputs' sizes are more general: they may be arrays of any
       * dimension */
      /*  and any size, so long as their sizes match each other excepting that
       * any */
      /*  input can have any dimension as a singleton. For example, s and t can
       * be */
      /*  3D arrays of size [nz, nx, ny], while p can be a vector of size
       * [nz,1]. */
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
      /*  coefficients nonlinear equation of state in pressure coordinates for
       */
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
      b_s1o2 = muDoubleScalarSqrt(s_B);
      fa = ((s * (s * (s * (s * (s * 6.536332E-9 + -1.120083E-6) +
                            0.0001001685) +
                       -0.00909529) +
                  0.06793952) +
             999.842594) +
            m * (((s * (s * (s * (s * 5.3875E-9 + -8.2467E-7) + 7.6438E-5) +
                        -0.0040899) +
                   0.824493) +
                  s1o2 * (s * (s * -1.6546E-6 + 0.00010227) + -0.00572466)) +
                 m * 0.00048314)) /
               (1.0 -
                z_tmp /
                    (((s * (s * (s * (s * -0.0004190253 + 0.09648704) +
                                 -17.06103) +
                            1444.304) +
                       196593.3) +
                      m * ((s * (s * (s * -0.0005084188 + 0.06283263) +
                                 -3.101089) +
                            528.4855) +
                           s1o2 * (s * (s * -0.004619924 + 0.09085835) +
                                   3.88664))) +
                     z_tmp *
                         (((s * (s * (s * 1.956415E-6 + -0.0002984642) +
                                 0.02212276) +
                            3.186519) +
                           m * ((s * (s * 2.059331E-7 + -0.0001847318) +
                                 0.006704388) +
                                s1o2 * 0.0001480266)) +
                          z_tmp * ((s * (s * 1.39468E-8 + -1.202016E-6) +
                                    2.102898E-5) +
                                   m * (s * (s * 6.207323E-11 + 6.128773E-9) +
                                        -2.040237E-7))))) -
           ((t_B * (t_B * (t_B * (t_B * (t_B * 6.536332E-9 + -1.120083E-6) +
                                  0.0001001685) +
                           -0.00909529) +
                    0.06793952) +
             999.842594) +
            s_B * (((t_B * (t_B * (t_B * (t_B * 5.3875E-9 + -8.2467E-7) +
                                   7.6438E-5) +
                            -0.0040899) +
                     0.824493) +
                    b_s1o2 *
                        (t_B * (t_B * -1.6546E-6 + 0.00010227) + -0.00572466)) +
                   s_B * 0.00048314)) /
               (1.0 -
                z_tmp /
                    (((t_B * (t_B * (t_B * (t_B * -0.0004190253 + 0.09648704) +
                                     -17.06103) +
                              1444.304) +
                       196593.3) +
                      s_B * ((t_B * (t_B * (t_B * -0.0005084188 + 0.06283263) +
                                     -3.101089) +
                              528.4855) +
                             b_s1o2 * (t_B * (t_B * -0.004619924 + 0.09085835) +
                                       3.88664))) +
                     z_tmp *
                         (((t_B * (t_B * (t_B * 1.956415E-6 + -0.0002984642) +
                                   0.02212276) +
                            3.186519) +
                           s_B * ((t_B * (t_B * 2.059331E-7 + -0.0001847318) +
                                   0.006704388) +
                                  b_s1o2 * 0.0001480266)) +
                          z_tmp *
                              ((t_B * (t_B * 1.39468E-8 + -1.202016E-6) +
                                2.102898E-5) +
                               s_B * (t_B * (t_B * 6.207323E-11 + 6.128773E-9) +
                                      -2.040237E-7)))));
      /*  Evaluate difference between (a) eos at location on cast A with */
      /*  properties (S_A, T_A, P_A) where the pressure or depth is p_avg +
       * dp/2, */
      /*  and (b) eos at location on cast B with properties (S_B, T_B, P_B)
       * where */
      /*  the pressure or depth is p_avg - dp/2. Here, eos is always evaluated
       * at */
      /*  the mid-pressure or mid-depth, p_avg. */
      ppc_val2(P_A_data, P_A_size, Sppc_A_data, Sppc_A_size, Tppc_A_data,
               p_avg + 0.5 * dp, &m, &s);
      ppc_val2(P_B_data, P_B_size, Sppc_B_data, Sppc_B_size, Tppc_B_data,
               p_avg - 0.5 * dp, &s_B, &t_B);
      /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
      /*  */
      /*  */
      /*  rho = eoscg_densjmd95_bsq(s,t,z)                             [ kg /
       * m^3 ] */
      /*  computes the Boussinesq JMD95 in-situ density, given practical
       * salinity */
      /*  s, potential temperature t, and depth z [m, positive and increasing */
      /*  down].  The depth z is converted into hydrostatic pressure for a given
       */
      /*  gravitational acceleration and Boussinesq reference density, which are
       */
      /*  hard-coded into this function (edit variables grav and rhob as
       * needed). */
      /*  */
      /*  */
      /*  This function is derived from densjmd95.m, documented below. Input
       * checks */
      /*  and expansion of variables have been removed; instead automatic
       * expansion */
      /*  (requiring MATLAB 2016b or later) is used. The calculation has also
       * been */
      /*  streamlined by pre-allocating arrays for coeffcients, and modifying
       * the */
      /*  coefficients such that pressure need not be internally converted from
       */
      /*  dbar to bar. This function is compatible with MATLAB's codegen. */
      /*  */
      /*  The units of s and t are as in densjmd95, documented below. */
      /*   */
      /*  The inputs' sizes are more general: they may be arrays of any
       * dimension */
      /*  and any size, so long as their sizes match each other excepting that
       * any */
      /*  input can have any dimension as a singleton. For example, s and t can
       * be */
      /*  3D arrays of size [nz, nx, ny], while p can be a vector of size
       * [nz,1]. */
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
      /*  coefficients nonlinear equation of state in pressure coordinates for
       */
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
      s1o2 = muDoubleScalarSqrt(m);
      /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
      /*  */
      /*  */
      /*  rho = eoscg_densjmd95_bsq(s,t,z)                             [ kg /
       * m^3 ] */
      /*  computes the Boussinesq JMD95 in-situ density, given practical
       * salinity */
      /*  s, potential temperature t, and depth z [m, positive and increasing */
      /*  down].  The depth z is converted into hydrostatic pressure for a given
       */
      /*  gravitational acceleration and Boussinesq reference density, which are
       */
      /*  hard-coded into this function (edit variables grav and rhob as
       * needed). */
      /*  */
      /*  */
      /*  This function is derived from densjmd95.m, documented below. Input
       * checks */
      /*  and expansion of variables have been removed; instead automatic
       * expansion */
      /*  (requiring MATLAB 2016b or later) is used. The calculation has also
       * been */
      /*  streamlined by pre-allocating arrays for coeffcients, and modifying
       * the */
      /*  coefficients such that pressure need not be internally converted from
       */
      /*  dbar to bar. This function is compatible with MATLAB's codegen. */
      /*  */
      /*  The units of s and t are as in densjmd95, documented below. */
      /*   */
      /*  The inputs' sizes are more general: they may be arrays of any
       * dimension */
      /*  and any size, so long as their sizes match each other excepting that
       * any */
      /*  input can have any dimension as a singleton. For example, s and t can
       * be */
      /*  3D arrays of size [nz, nx, ny], while p can be a vector of size
       * [nz,1]. */
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
      /*  coefficients nonlinear equation of state in pressure coordinates for
       */
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
      b_s1o2 = muDoubleScalarSqrt(s_B);
      dxp =
          ((s * (s * (s * (s * (s * 6.536332E-9 + -1.120083E-6) +
                           0.0001001685) +
                      -0.00909529) +
                 0.06793952) +
            999.842594) +
           m * (((s * (s * (s * (s * 5.3875E-9 + -8.2467E-7) + 7.6438E-5) +
                       -0.0040899) +
                  0.824493) +
                 s1o2 * (s * (s * -1.6546E-6 + 0.00010227) + -0.00572466)) +
                m * 0.00048314)) /
              (1.0 -
               z_tmp /
                   (((s * (s * (s * (s * -0.0004190253 + 0.09648704) +
                                -17.06103) +
                           1444.304) +
                      196593.3) +
                     m * ((s * (s * (s * -0.0005084188 + 0.06283263) +
                                -3.101089) +
                           528.4855) +
                          s1o2 * (s * (s * -0.004619924 + 0.09085835) +
                                  3.88664))) +
                    z_tmp *
                        (((s * (s * (s * 1.956415E-6 + -0.0002984642) +
                                0.02212276) +
                           3.186519) +
                          m * ((s * (s * 2.059331E-7 + -0.0001847318) +
                                0.006704388) +
                               s1o2 * 0.0001480266)) +
                         z_tmp * ((s * (s * 1.39468E-8 + -1.202016E-6) +
                                   2.102898E-5) +
                                  m * (s * (s * 6.207323E-11 + 6.128773E-9) +
                                       -2.040237E-7))))) -
          ((t_B * (t_B * (t_B * (t_B * (t_B * 6.536332E-9 + -1.120083E-6) +
                                 0.0001001685) +
                          -0.00909529) +
                   0.06793952) +
            999.842594) +
           s_B * (((t_B * (t_B * (t_B * (t_B * 5.3875E-9 + -8.2467E-7) +
                                  7.6438E-5) +
                           -0.0040899) +
                    0.824493) +
                   b_s1o2 *
                       (t_B * (t_B * -1.6546E-6 + 0.00010227) + -0.00572466)) +
                  s_B * 0.00048314)) /
              (1.0 -
               z_tmp /
                   (((t_B * (t_B * (t_B * (t_B * -0.0004190253 + 0.09648704) +
                                    -17.06103) +
                             1444.304) +
                      196593.3) +
                     s_B * ((t_B * (t_B * (t_B * -0.0005084188 + 0.06283263) +
                                    -3.101089) +
                             528.4855) +
                            b_s1o2 * (t_B * (t_B * -0.004619924 + 0.09085835) +
                                      3.88664))) +
                    z_tmp *
                        (((t_B * (t_B * (t_B * 1.956415E-6 + -0.0002984642) +
                                  0.02212276) +
                           3.186519) +
                          s_B * ((t_B * (t_B * 2.059331E-7 + -0.0001847318) +
                                  0.006704388) +
                                 b_s1o2 * 0.0001480266)) +
                         z_tmp *
                             ((t_B * (t_B * 1.39468E-8 + -1.202016E-6) +
                               2.102898E-5) +
                              s_B * (t_B * (t_B * 6.207323E-11 + 6.128773E-9) +
                                     -2.040237E-7)))));
      c = lb;
      fc = fa;
      e = dp - lb;
      d = e;
      do {
        exitg1 = 0;
        if (muDoubleScalarAbs(fc) < muDoubleScalarAbs(dxp)) {
          lb = dp;
          dp = c;
          c = lb;
          fa = dxp;
          dxp = fc;
          fc = fa;
        }
        dxm = 4.4408920985006262E-16 * muDoubleScalarAbs(dp) + tolp;
        m = 0.5 * (c - dp);
        if ((muDoubleScalarAbs(m) <= dxm) || (dxp == 0.0)) {
          exitg1 = 1;
        } else {
          if ((muDoubleScalarAbs(e) < dxm) ||
              (muDoubleScalarAbs(fa) <= muDoubleScalarAbs(dxp))) {
            e = m;
            d = m;
          } else {
            s = dxp / fa;
            if (lb == c) {
              fa = 2.0 * m * s;
              q = 1.0 - s;
            } else {
              q = fa / fc;
              r = dxp / fc;
              fa = s * (2.0 * m * q * (q - r) - (dp - lb) * (r - 1.0));
              q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (0.0 < fa) {
              q = -q;
            } else {
              fa = -fa;
            }
            s = e;
            e = d;
            if ((2.0 * fa < 3.0 * m * q - muDoubleScalarAbs(dxm * q)) &&
                (fa < muDoubleScalarAbs(0.5 * s * q))) {
              d = fa / q;
            } else {
              e = m;
              d = m;
            }
          }
          lb = dp;
          fa = dxp;
          if (dxm < muDoubleScalarAbs(d)) {
            dp += d;
          } else if (0.0 < m) {
            dp += dxm;
          } else {
            dp -= dxm;
          }
          /*  Evaluate difference between (a) eos at location on cast A with */
          /*  properties (S_A, T_A, P_A) where the pressure or depth is p_avg +
           * dp/2, */
          /*  and (b) eos at location on cast B with properties (S_B, T_B, P_B)
           * where */
          /*  the pressure or depth is p_avg - dp/2. Here, eos is always
           * evaluated at */
          /*  the mid-pressure or mid-depth, p_avg. */
          ppc_val2(P_A_data, P_A_size, Sppc_A_data, Sppc_A_size, Tppc_A_data,
                   p_avg + 0.5 * dp, &m, &s);
          ppc_val2(P_B_data, P_B_size, Sppc_B_data, Sppc_B_size, Tppc_B_data,
                   p_avg - 0.5 * dp, &s_B, &t_B);
          /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
          /*  */
          /*  */
          /*  rho = eoscg_densjmd95_bsq(s,t,z)                             [ kg
           * / m^3 ] */
          /*  computes the Boussinesq JMD95 in-situ density, given practical
           * salinity */
          /*  s, potential temperature t, and depth z [m, positive and
           * increasing */
          /*  down].  The depth z is converted into hydrostatic pressure for a
           * given */
          /*  gravitational acceleration and Boussinesq reference density, which
           * are */
          /*  hard-coded into this function (edit variables grav and rhob as
           * needed). */
          /*  */
          /*  */
          /*  This function is derived from densjmd95.m, documented below. Input
           * checks */
          /*  and expansion of variables have been removed; instead automatic
           * expansion */
          /*  (requiring MATLAB 2016b or later) is used. The calculation has
           * also been */
          /*  streamlined by pre-allocating arrays for coeffcients, and
           * modifying the */
          /*  coefficients such that pressure need not be internally converted
           * from */
          /*  dbar to bar. This function is compatible with MATLAB's codegen. */
          /*  */
          /*  The units of s and t are as in densjmd95, documented below. */
          /*   */
          /*  The inputs' sizes are more general: they may be arrays of any
           * dimension */
          /*  and any size, so long as their sizes match each other excepting
           * that any */
          /*  input can have any dimension as a singleton. For example, s and t
           * can be */
          /*  3D arrays of size [nz, nx, ny], while p can be a vector of size
           * [nz,1]. */
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
          /*     Density of Sea Water using Jackett and McDougall 1995 (JAOT 12)
           */
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
          /*  Convert depth [m] to the hydrostatic pressure [dbar] implied by
           * the */
          /*  following two hard-coded parameters (edit as needed): */
          /*  gravitational acceleration [m /s^2] */
          /*  Boussinesq reference density [kg / m^3] */
          /*  Pascal to dbar conversion [dbar / Pa] */
          /*  depth to pressure conversion [dbar / m] */
          /*  Henceforth z is actually pressure [dbar] */
          /*  coefficients nonlinear equation of state in pressure coordinates
           * for */
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
          s1o2 = muDoubleScalarSqrt(m);
          /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
          /*  */
          /*  */
          /*  rho = eoscg_densjmd95_bsq(s,t,z)                             [ kg
           * / m^3 ] */
          /*  computes the Boussinesq JMD95 in-situ density, given practical
           * salinity */
          /*  s, potential temperature t, and depth z [m, positive and
           * increasing */
          /*  down].  The depth z is converted into hydrostatic pressure for a
           * given */
          /*  gravitational acceleration and Boussinesq reference density, which
           * are */
          /*  hard-coded into this function (edit variables grav and rhob as
           * needed). */
          /*  */
          /*  */
          /*  This function is derived from densjmd95.m, documented below. Input
           * checks */
          /*  and expansion of variables have been removed; instead automatic
           * expansion */
          /*  (requiring MATLAB 2016b or later) is used. The calculation has
           * also been */
          /*  streamlined by pre-allocating arrays for coeffcients, and
           * modifying the */
          /*  coefficients such that pressure need not be internally converted
           * from */
          /*  dbar to bar. This function is compatible with MATLAB's codegen. */
          /*  */
          /*  The units of s and t are as in densjmd95, documented below. */
          /*   */
          /*  The inputs' sizes are more general: they may be arrays of any
           * dimension */
          /*  and any size, so long as their sizes match each other excepting
           * that any */
          /*  input can have any dimension as a singleton. For example, s and t
           * can be */
          /*  3D arrays of size [nz, nx, ny], while p can be a vector of size
           * [nz,1]. */
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
          /*     Density of Sea Water using Jackett and McDougall 1995 (JAOT 12)
           */
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
          /*  Convert depth [m] to the hydrostatic pressure [dbar] implied by
           * the */
          /*  following two hard-coded parameters (edit as needed): */
          /*  gravitational acceleration [m /s^2] */
          /*  Boussinesq reference density [kg / m^3] */
          /*  Pascal to dbar conversion [dbar / Pa] */
          /*  depth to pressure conversion [dbar / m] */
          /*  Henceforth z is actually pressure [dbar] */
          /*  coefficients nonlinear equation of state in pressure coordinates
           * for */
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
          b_s1o2 = muDoubleScalarSqrt(s_B);
          dxp =
              ((s * (s * (s * (s * (s * 6.536332E-9 + -1.120083E-6) +
                               0.0001001685) +
                          -0.00909529) +
                     0.06793952) +
                999.842594) +
               m * (((s * (s * (s * (s * 5.3875E-9 + -8.2467E-7) + 7.6438E-5) +
                           -0.0040899) +
                      0.824493) +
                     s1o2 * (s * (s * -1.6546E-6 + 0.00010227) + -0.00572466)) +
                    m * 0.00048314)) /
                  (1.0 -
                   z_tmp /
                       (((s * (s * (s * (s * -0.0004190253 + 0.09648704) +
                                    -17.06103) +
                               1444.304) +
                          196593.3) +
                         m * ((s * (s * (s * -0.0005084188 + 0.06283263) +
                                    -3.101089) +
                               528.4855) +
                              s1o2 * (s * (s * -0.004619924 + 0.09085835) +
                                      3.88664))) +
                        z_tmp * (((s * (s * (s * 1.956415E-6 + -0.0002984642) +
                                        0.02212276) +
                                   3.186519) +
                                  m * ((s * (s * 2.059331E-7 + -0.0001847318) +
                                        0.006704388) +
                                       s1o2 * 0.0001480266)) +
                                 z_tmp * ((s * (s * 1.39468E-8 + -1.202016E-6) +
                                           2.102898E-5) +
                                          m * (s * (s * 6.207323E-11 +
                                                    6.128773E-9) +
                                               -2.040237E-7))))) -
              ((t_B * (t_B * (t_B * (t_B * (t_B * 6.536332E-9 + -1.120083E-6) +
                                     0.0001001685) +
                              -0.00909529) +
                       0.06793952) +
                999.842594) +
               s_B * (((t_B * (t_B * (t_B * (t_B * 5.3875E-9 + -8.2467E-7) +
                                      7.6438E-5) +
                               -0.0040899) +
                        0.824493) +
                       b_s1o2 * (t_B * (t_B * -1.6546E-6 + 0.00010227) +
                                 -0.00572466)) +
                      s_B * 0.00048314)) /
                  (1.0 -
                   z_tmp /
                       (((t_B *
                              (t_B * (t_B * (t_B * -0.0004190253 + 0.09648704) +
                                      -17.06103) +
                               1444.304) +
                          196593.3) +
                         s_B *
                             ((t_B * (t_B * (t_B * -0.0005084188 + 0.06283263) +
                                      -3.101089) +
                               528.4855) +
                              b_s1o2 *
                                  (t_B * (t_B * -0.004619924 + 0.09085835) +
                                   3.88664))) +
                        z_tmp *
                            (((t_B *
                                   (t_B * (t_B * 1.956415E-6 + -0.0002984642) +
                                    0.02212276) +
                               3.186519) +
                              s_B *
                                  ((t_B * (t_B * 2.059331E-7 + -0.0001847318) +
                                    0.006704388) +
                                   b_s1o2 * 0.0001480266)) +
                             z_tmp * ((t_B * (t_B * 1.39468E-8 + -1.202016E-6) +
                                       2.102898E-5) +
                                      s_B * (t_B * (t_B * 6.207323E-11 +
                                                    6.128773E-9) +
                                             -2.040237E-7)))));
          if (((0.0 < dxp) && (0.0 < fc)) || ((dxp <= 0.0) && (fc <= 0.0))) {
            c = lb;
            fc = fa;
            e = dp - lb;
            d = e;
          }
        }
      } while (exitg1 == 0);
    } else {
      dp = rtNaN;
    }
  } else {
    dp = rtNaN;
  }
  return dp;
}

/* End of code generation (ntp_midpoint_to_casts.c) */
