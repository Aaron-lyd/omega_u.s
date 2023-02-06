/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * omega_vertsolve.c
 *
 * Code generation for function 'omega_vertsolve'
 *
 */

/* Include files */
#include "omega_vertsolve.h"
#include "omega_vertsolve_data.h"
#include "omega_vertsolve_emxutil.h"
#include "omega_vertsolve_types.h"
#include "ppc_val2.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"

/* Function Definitions */
void omega_vertsolve(const emxArray_real_T *Sppc, const emxArray_real_T *Tppc,
                     const emxArray_real_T *P, const emxArray_real_T *BotK,
                     emxArray_real_T *s, emxArray_real_T *t, emxArray_real_T *p,
                     real_T tolp, const emxArray_real_T *phi)
{
  emxArray_real_T *Pn;
  const real_T *BotK_data;
  const real_T *P_data;
  const real_T *Sppc_data;
  const real_T *Tppc_data;
  const real_T *phi_data;
  real_T b_p;
  real_T b_s;
  real_T b_t;
  real_T d;
  real_T *Pn_data;
  real_T *p_data;
  real_T *s_data;
  real_T *t_data;
  int32_T Sppcn_size[2];
  int32_T i;
  int32_T i1;
  int32_T i2;
  int32_T n;
  boolean_T Pmat;
  phi_data = phi->data;
  p_data = p->data;
  t_data = t->data;
  s_data = s->data;
  BotK_data = BotK->data;
  P_data = P->data;
  Tppc_data = Tppc->data;
  Sppc_data = Sppc->data;
  emlrtHeapReferenceStackEnterFcnR2012b(emlrtRootTLSGlobal);
  /* OMEGA_VERTSOLVE  Root finding of a new surface with a specified a density
   */
  /*                  difference from the current surface */
  /*  */
  /*  [s, t, p] = omega_vertsolve(Sppc, Tppc, P, BotK, s0, t0, p0, tolp, phi) */
  /*  determines pressures p that satisfy */
  /*    |eos(S_n(p(n)), T_n(p(n)), p0(n)) + phi(n) - eos(s0(n), t0(n), p0(n))| <
   * tolp */
  /*  where S_n(p') and T_n(p') are interpolants whose coefficients are given */
  /*  by Sppc(:,:,n) and Tppc(:,:,n) in water column n. These interpolants  */
  /*  determine salinities s and temperatures t on the surface, namely */
  /*    s(n) = S_n(p(n)) and t(n) = T_n(p(n)). */
  /*  The function eos.m determines either the in-situ density or the specific
   */
  /*  volume. */
  /*  */
  /*  --- Input */
  /*  Sppc [O, K-1, N]: coefficients for piecewise polynomial for practical */
  /*                    / Absolute Salinity in terms of P */
  /*  Tppc [O, K-1, N]: coefficients for piecewise polynomial for potential */
  /*                    / Conservative Temperature in terms of P */
  /*  P [K, N]: knots for the pressure or depth of the casts */
  /*  BotK [1, N]: number of valid data points on each cast */
  /*  s [1, N]: practical / Absolute salinity on the initial surface */
  /*  t [1, N]: potential / Conservative temperature on the initial surface */
  /*  p [1, N]: pressure or depth on the initial surface */
  /*  tolp [1,1]: precision of solution in pressure or depth */
  /*  phi [1,N]: the desired in-situ density or specific volume change of the */
  /*             surface */
  /*  */
  /*  Note: O is the order of the piecewise polynomials */
  /*        K is the maximum number of knots in these piecewise polynomials, */
  /*            i.e. the maximum number of bottles in any cast */
  /*        N is the number of water columns (possibly including land). */
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
  /*  --- Output */
  /*  p [same as input p]: pressure or depth of the updated surface */
  /*  s [same as input p]: practical / Absolute salinity of the updated surface
   */
  /*  t [same as input p]: potential / Conservative temperature of the updated
   * surface */
  /*  */
  /*  --- Units */
  /*  The units of s, t, p, , T, P, tolp, and phi are determined by the */
  /*  function eos.m. */
  /*  Author(s) : Geoff Stanley */
  /*  Email     : g.stanley@unsw.edu.au */
  /*  Email     : geoffstanley@gmail.com */
  /*  Inputs s0, t0, and p0 are named s, t, p so operations are done in-place.
   */
  Pmat = (((P->size[0] != 1) && (P->size[1] != 1)) || (P->size[2] != 1));
  /*  Loop over each water column */
  i = p->size[0] * p->size[1];
  emxInit_real_T(&Pn, 1);
  for (n = 0; n < i; n++) {
    d = phi_data[n];
    if ((!muDoubleScalarIsNaN(phi_data[n])) && (BotK_data[n] > 1.0)) {
      real_T Sppcn_data[392];
      real_T Tppcn_data[392];
      real_T dxm;
      real_T dxp;
      real_T fb;
      real_T lb;
      real_T m;
      real_T rn;
      real_T s1o2;
      real_T tol;
      real_T ub;
      int32_T Sppc_idx_0;
      int32_T Sppc_idx_1;
      int32_T Tppcn_size_idx_0;
      int32_T b_loop_ub;
      int32_T exitg1;
      int32_T loop_ub;
      boolean_T fapos;
      /*  Select this water column */
      if (BotK_data[n] - 1.0 < 1.0) {
        loop_ub = 0;
      } else {
        loop_ub = (int32_T)(BotK_data[n] - 1.0);
      }
      Sppc_idx_0 = Sppc->size[0];
      Sppc_idx_1 = Sppc->size[1];
      Sppcn_size[0] = Sppc->size[0];
      Sppcn_size[1] = loop_ub;
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_loop_ub = Sppc->size[0];
        for (i2 = 0; i2 < b_loop_ub; i2++) {
          Sppcn_data[i2 + Sppcn_size[0] * i1] =
              Sppc_data[(i2 + Sppc_idx_0 * i1) + Sppc_idx_0 * Sppc_idx_1 * n];
        }
      }
      if (BotK_data[n] - 1.0 < 1.0) {
        loop_ub = 0;
      } else {
        loop_ub = (int32_T)(BotK_data[n] - 1.0);
      }
      Sppc_idx_0 = Tppc->size[0];
      Sppc_idx_1 = Tppc->size[1];
      Tppcn_size_idx_0 = Tppc->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_loop_ub = Tppc->size[0];
        for (i2 = 0; i2 < b_loop_ub; i2++) {
          Tppcn_data[i2 + Tppcn_size_idx_0 * i1] =
              Tppc_data[(i2 + Sppc_idx_0 * i1) + Sppc_idx_0 * Sppc_idx_1 * n];
        }
      }
      if (Pmat) {
        loop_ub = P->size[0];
        Sppc_idx_0 = (int32_T)BotK_data[n];
        i1 = Pn->size[0];
        Pn->size[0] = Sppc_idx_0;
        emxEnsureCapacity_real_T(Pn, i1);
        Pn_data = Pn->data;
        for (i1 = 0; i1 < Sppc_idx_0; i1++) {
          Pn_data[i1] = P_data[i1 + loop_ub * n];
        }
      } else {
        loop_ub = (int32_T)(BotK_data[n] - 1.0);
        i1 = Pn->size[0];
        Pn->size[0] = loop_ub + 1;
        emxEnsureCapacity_real_T(Pn, i1);
        Pn_data = Pn->data;
        for (i1 = 0; i1 <= loop_ub; i1++) {
          Pn_data[i1] = P_data[i1];
        }
        /*  .' is for codegen, so P and (1:k).' both column vectors */
      }
      /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
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
      /*  Author(s)       : Geoff Stanley */
      /*  Email           : g.stanley@unsw.edu.au */
      /*  Email           : geoffstanley@gmail.com */
      /*  Version         : 1.0 */
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
      b_p = p_data[n] * 1.015335;
      /*  Henceforth z is actually pressure [dbar] */
      /*  coefficients nonlinear equation of state in pressure coordinates for
       */
      /*  1. density of fresh water at p = 0 */
      /*  2. density of sea water at p = 0 */
      /*  coefficients in pressure coordinates for */
      /*  3. secant bulk modulus K of fresh water at p = 0 */
      /*  == original * 10 */
      /*  4. secant bulk modulus K of sea water at p = 0 */
      /*  == original * 10 */
      /*  5. secant bulk modulus K of sea water at p */
      /*  == original / 10 */
      s1o2 = muDoubleScalarSqrt(s_data[n]);
      rn =
          ((t_data[n] *
                (t_data[n] *
                     (t_data[n] *
                          (t_data[n] * (t_data[n] * 6.536332E-9 - 1.120083E-6) +
                           0.0001001685) -
                      0.00909529) +
                 0.06793952) +
            999.842594) +
           s_data[n] *
               (((t_data[n] * (t_data[n] * (t_data[n] * (t_data[n] * 5.3875E-9 -
                                                         8.2467E-7) +
                                            7.6438E-5) -
                               0.0040899) +
                  0.824493) +
                 s1o2 * (t_data[n] * (t_data[n] * -1.6546E-6 + 0.00010227) -
                         0.00572466)) +
                s_data[n] * 0.00048314)) /
          (1.0 -
           b_p /
               (((t_data[n] *
                      (t_data[n] * (t_data[n] * (t_data[n] * -0.0004190253 +
                                                 0.09648704) -
                                    17.06103) +
                       1444.304) +
                  196593.3) +
                 s_data[n] *
                     ((t_data[n] * (t_data[n] * (t_data[n] * -0.0005084188 +
                                                 0.06283263) -
                                    3.101089) +
                       528.4855) +
                      s1o2 *
                          (t_data[n] * (t_data[n] * -0.004619924 + 0.09085835) +
                           3.88664))) +
                b_p *
                    (((t_data[n] * (t_data[n] * (t_data[n] * 1.956415E-6 -
                                                 0.0002984642) +
                                    0.02212276) +
                       3.186519) +
                      s_data[n] * ((t_data[n] * (t_data[n] * 2.059331E-7 -
                                                 0.0001847318) +
                                    0.006704388) +
                                   s1o2 * 0.0001480266)) +
                     b_p *
                         ((t_data[n] * (t_data[n] * 1.39468E-8 - 1.202016E-6) +
                           2.102898E-5) +
                          s_data[n] * (t_data[n] * (t_data[n] * 6.207323E-11 +
                                                    6.128773E-9) -
                                       2.040237E-7)))));
      /*  Search for a sign-change, expanding outward from an initial guess */
      /* [lb, ub] = fzero_guess_to_bounds(@diff_avgx, p(n), Pn(1), Pn(k), ... */
      /*   Sppcn, Tppcn, Pn, s(n), t(n), p(n), d);  % DEV:  testing reference p
       * = average of current and point-to-be */
      fb = p_data[n];
      tol = Pn_data[0];
      m = Pn_data[(int32_T)BotK_data[n] - 1];
      b_p = p_data[n];
      /* FZERO_GUESS_TO_BOUNDS  Search for a sign change bounding a zero of a */
      /*                        univariate function, expanding geometrically */
      /*                        outward from an initial guess. */
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
      /*  --- Input: */
      /*    f       : handle to a function that accepts a real scalar as its
       * first */
      /*              input and returns a real scalar */
      /*    x       : initial scalar guess for a root of f */
      /*    A       : scalar lower bound */
      /*    B       : scalar upper bound */
      /*    varargin: All additional inputs are passed directly to f */
      /*  */
      /*  --- Output: */
      /*    a : lower bound for interval containing a root, scalar */
      /*    b : upper bound for interval containing a root, scalar */
      /*  */
      /*  --- Acknowledgements: */
      /*  Expansion from initial guess inspired by MATLAB's fzero.m. */
      /*  */
      /*  Author    : Geoff Stanley */
      /*  Email     : geoffstanley@gmail.com */
      /*  Version   : 1.0 */
      /*  History   : 01/07/2020 - initial release */
      /*            : 22/07/2020 - fix infinite loop in bounded case, arising
       * from machine precision rounding */
      /*  Geometrically expand from the guess x, until a sign change is found */
      /*  Handle bad inputs */
      fapos = muDoubleScalarIsNaN(Pn_data[0]);
      if (fapos || muDoubleScalarIsNaN(m) || muDoubleScalarIsNaN(fb)) {
        lb = rtNaN;
        ub = rtNaN;
      } else {
        /*  Evaluate difference between (a) eos at location on the cast where
         * the */
        /*  pressure or depth is p, and (b) eos at location on the cast where
         * the */
        /*  pressure or depth is p0 (where the surface currently is) plus the
         * density */
        /*  perturbation d.  Part (b) is precomputed as r0.  Here, eos always */
        /*  evaluated at the pressure or depth of the original position, p0;
         * this is */
        /*  to calculate locally referenced potential density with reference
         * pressure */
        /*  p0. */
        /*  Interpolate S and T to the current pressure or depth */
        ppc_val2(Pn, Sppcn_data, Sppcn_size, Tppcn_data, tol, &b_s, &b_t);
        /*  Calculate the potential density or potential specific volume
         * difference */
        /* out =  eos(s, t, p0) - d - r0 ; */
        /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
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
        /*  Author(s)       : Geoff Stanley */
        /*  Email           : g.stanley@unsw.edu.au */
        /*  Email           : geoffstanley@gmail.com */
        /*  Version         : 1.0 */
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
        b_p *= 1.015335;
        /*  Henceforth z is actually pressure [dbar] */
        /*  coefficients nonlinear equation of state in pressure coordinates for
         */
        /*  1. density of fresh water at p = 0 */
        /*  2. density of sea water at p = 0 */
        /*  coefficients in pressure coordinates for */
        /*  3. secant bulk modulus K of fresh water at p = 0 */
        /*  == original * 10 */
        /*  4. secant bulk modulus K of sea water at p = 0 */
        /*  == original * 10 */
        /*  5. secant bulk modulus K of sea water at p */
        /*  == original / 10 */
        s1o2 = muDoubleScalarSqrt(b_s);
        if ((((b_t * (b_t * (b_t * (b_t * (b_t * 6.536332E-9 - 1.120083E-6) +
                                    0.0001001685) -
                             0.00909529) +
                      0.06793952) +
               999.842594) +
              b_s * (((b_t * (b_t * (b_t * (b_t * 5.3875E-9 - 8.2467E-7) +
                                     7.6438E-5) -
                              0.0040899) +
                       0.824493) +
                      s1o2 * (b_t * (b_t * -1.6546E-6 + 0.00010227) -
                              0.00572466)) +
                     b_s * 0.00048314)) /
                 (1.0 -
                  b_p /
                      (((b_t *
                             (b_t * (b_t * (b_t * -0.0004190253 + 0.09648704) -
                                     17.06103) +
                              1444.304) +
                         196593.3) +
                        b_s *
                            ((b_t * (b_t * (b_t * -0.0005084188 + 0.06283263) -
                                     3.101089) +
                              528.4855) +
                             s1o2 * (b_t * (b_t * -0.004619924 + 0.09085835) +
                                     3.88664))) +
                       b_p *
                           (((b_t * (b_t * (b_t * 1.956415E-6 - 0.0002984642) +
                                     0.02212276) +
                              3.186519) +
                             b_s * ((b_t * (b_t * 2.059331E-7 - 0.0001847318) +
                                     0.006704388) +
                                    s1o2 * 0.0001480266)) +
                            b_p * ((b_t * (b_t * 1.39468E-8 - 1.202016E-6) +
                                    2.102898E-5) +
                                   b_s * (b_t * (b_t * 6.207323E-11 +
                                                 6.128773E-9) -
                                          2.040237E-7))))) -
             rn) -
                d ==
            0.0) {
          lb = tol;
          ub = tol;
        } else {
          /*  Evaluate difference between (a) eos at location on the cast where
           * the */
          /*  pressure or depth is p, and (b) eos at location on the cast where
           * the */
          /*  pressure or depth is p0 (where the surface currently is) plus the
           * density */
          /*  perturbation d.  Part (b) is precomputed as r0.  Here, eos always
           */
          /*  evaluated at the pressure or depth of the original position, p0;
           * this is */
          /*  to calculate locally referenced potential density with reference
           * pressure */
          /*  p0. */
          /*  Interpolate S and T to the current pressure or depth */
          ppc_val2(Pn, Sppcn_data, Sppcn_size, Tppcn_data, m, &b_s, &b_t);
          /*  Calculate the potential density or potential specific volume
           * difference */
          /* out =  eos(s, t, p0) - d - r0 ; */
          /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
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
          /*  Author(s)       : Geoff Stanley */
          /*  Email           : g.stanley@unsw.edu.au */
          /*  Email           : geoffstanley@gmail.com */
          /*  Version         : 1.0 */
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
          /*  4. secant bulk modulus K of sea water at p = 0 */
          /*  == original * 10 */
          /*  5. secant bulk modulus K of sea water at p */
          /*  == original / 10 */
          s1o2 = muDoubleScalarSqrt(b_s);
          if ((((b_t * (b_t * (b_t * (b_t * (b_t * 6.536332E-9 - 1.120083E-6) +
                                      0.0001001685) -
                               0.00909529) +
                        0.06793952) +
                 999.842594) +
                b_s * (((b_t * (b_t * (b_t * (b_t * 5.3875E-9 - 8.2467E-7) +
                                       7.6438E-5) -
                                0.0040899) +
                         0.824493) +
                        s1o2 * (b_t * (b_t * -1.6546E-6 + 0.00010227) -
                                0.00572466)) +
                       b_s * 0.00048314)) /
                   (1.0 -
                    b_p /
                        (((b_t * (b_t * (b_t * (b_t * -0.0004190253 +
                                                0.09648704) -
                                         17.06103) +
                                  1444.304) +
                           196593.3) +
                          b_s *
                              ((b_t *
                                    (b_t * (b_t * -0.0005084188 + 0.06283263) -
                                     3.101089) +
                                528.4855) +
                               s1o2 * (b_t * (b_t * -0.004619924 + 0.09085835) +
                                       3.88664))) +
                         b_p * (((b_t * (b_t * (b_t * 1.956415E-6 -
                                                0.0002984642) +
                                         0.02212276) +
                                  3.186519) +
                                 b_s * ((b_t * (b_t * 2.059331E-7 -
                                                0.0001847318) +
                                         0.006704388) +
                                        s1o2 * 0.0001480266)) +
                                b_p * ((b_t * (b_t * 1.39468E-8 - 1.202016E-6) +
                                        2.102898E-5) +
                                       b_s * (b_t * (b_t * 6.207323E-11 +
                                                     6.128773E-9) -
                                              2.040237E-7))))) -
               rn) -
                  d ==
              0.0) {
            lb = m;
            ub = m;
          } else {
            boolean_T fbpos;
            fb = muDoubleScalarMin(muDoubleScalarMax(fb, tol), m);
            /*  bounds are given */
            dxp = (m - fb) / 50.0;
            dxm = (fb - tol) / 50.0;
            /*  Set a = x, except when x is so close to A that machine roundoff
             * makes dxm identically 0, */
            /*  which would lead to an infinite loop below.  In this case, set a
             * = A. */
            if (dxm == 0.0) {
              lb = tol;
            } else {
              lb = fb;
            }
            /*  Evaluate difference between (a) eos at location on the cast
             * where the */
            /*  pressure or depth is p, and (b) eos at location on the cast
             * where the */
            /*  pressure or depth is p0 (where the surface currently is) plus
             * the density */
            /*  perturbation d.  Part (b) is precomputed as r0.  Here, eos
             * always */
            /*  evaluated at the pressure or depth of the original position, p0;
             * this is */
            /*  to calculate locally referenced potential density with reference
             * pressure */
            /*  p0. */
            /*  Interpolate S and T to the current pressure or depth */
            ppc_val2(Pn, Sppcn_data, Sppcn_size, Tppcn_data, lb, &b_s, &b_t);
            /*  Calculate the potential density or potential specific volume
             * difference */
            /* out =  eos(s, t, p0) - d - r0 ; */
            /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
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
            /*  Author(s)       : Geoff Stanley */
            /*  Email           : g.stanley@unsw.edu.au */
            /*  Email           : geoffstanley@gmail.com */
            /*  Version         : 1.0 */
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
            /*  4. secant bulk modulus K of sea water at p = 0 */
            /*  == original * 10 */
            /*  5. secant bulk modulus K of sea water at p */
            /*  == original / 10 */
            s1o2 = muDoubleScalarSqrt(b_s);
            fapos =
                ((((b_t *
                        (b_t * (b_t * (b_t * (b_t * 6.536332E-9 - 1.120083E-6) +
                                       0.0001001685) -
                                0.00909529) +
                         0.06793952) +
                    999.842594) +
                   b_s * (((b_t * (b_t * (b_t * (b_t * 5.3875E-9 - 8.2467E-7) +
                                          7.6438E-5) -
                                   0.0040899) +
                            0.824493) +
                           s1o2 * (b_t * (b_t * -1.6546E-6 + 0.00010227) -
                                   0.00572466)) +
                          b_s * 0.00048314)) /
                      (1.0 -
                       b_p / (((b_t * (b_t * (b_t * (b_t * -0.0004190253 +
                                                     0.09648704) -
                                              17.06103) +
                                       1444.304) +
                                196593.3) +
                               b_s * ((b_t * (b_t * (b_t * -0.0005084188 +
                                                     0.06283263) -
                                              3.101089) +
                                       528.4855) +
                                      s1o2 * (b_t * (b_t * -0.004619924 +
                                                     0.09085835) +
                                              3.88664))) +
                              b_p * (((b_t * (b_t * (b_t * 1.956415E-6 -
                                                     0.0002984642) +
                                              0.02212276) +
                                       3.186519) +
                                      b_s * ((b_t * (b_t * 2.059331E-7 -
                                                     0.0001847318) +
                                              0.006704388) +
                                             s1o2 * 0.0001480266)) +
                                     b_p * ((b_t * (b_t * 1.39468E-8 -
                                                    1.202016E-6) +
                                             2.102898E-5) +
                                            b_s * (b_t * (b_t * 6.207323E-11 +
                                                          6.128773E-9) -
                                                   2.040237E-7))))) -
                  rn) -
                     d >
                 0.0);
            /*  Similarly, set b = x, except for machine precision problems. */
            if (dxp == 0.0) {
              ub = m;
              /*  Evaluate difference between (a) eos at location on the cast
               * where the */
              /*  pressure or depth is p, and (b) eos at location on the cast
               * where the */
              /*  pressure or depth is p0 (where the surface currently is) plus
               * the density */
              /*  perturbation d.  Part (b) is precomputed as r0.  Here, eos
               * always */
              /*  evaluated at the pressure or depth of the original position,
               * p0; this is */
              /*  to calculate locally referenced potential density with
               * reference pressure */
              /*  p0. */
              /*  Interpolate S and T to the current pressure or depth */
              ppc_val2(Pn, Sppcn_data, Sppcn_size, Tppcn_data, m, &b_s, &b_t);
              /*  Calculate the potential density or potential specific volume
               * difference */
              /* out =  eos(s, t, p0) - d - r0 ; */
              /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
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
              /*  Author(s)       : Geoff Stanley */
              /*  Email           : g.stanley@unsw.edu.au */
              /*  Email           : geoffstanley@gmail.com */
              /*  Version         : 1.0 */
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
              /*  4. secant bulk modulus K of sea water at p = 0 */
              /*  == original * 10 */
              /*  5. secant bulk modulus K of sea water at p */
              /*  == original / 10 */
              s1o2 = muDoubleScalarSqrt(b_s);
              fbpos =
                  ((((b_t * (b_t * (b_t * (b_t * (b_t * 6.536332E-9 -
                                                  1.120083E-6) +
                                           0.0001001685) -
                                    0.00909529) +
                             0.06793952) +
                      999.842594) +
                     b_s *
                         (((b_t * (b_t * (b_t * (b_t * 5.3875E-9 - 8.2467E-7) +
                                          7.6438E-5) -
                                   0.0040899) +
                            0.824493) +
                           s1o2 * (b_t * (b_t * -1.6546E-6 + 0.00010227) -
                                   0.00572466)) +
                          b_s * 0.00048314)) /
                        (1.0 -
                         b_p / (((b_t * (b_t * (b_t * (b_t * -0.0004190253 +
                                                       0.09648704) -
                                                17.06103) +
                                         1444.304) +
                                  196593.3) +
                                 b_s * ((b_t * (b_t * (b_t * -0.0005084188 +
                                                       0.06283263) -
                                                3.101089) +
                                         528.4855) +
                                        s1o2 * (b_t * (b_t * -0.004619924 +
                                                       0.09085835) +
                                                3.88664))) +
                                b_p * (((b_t * (b_t * (b_t * 1.956415E-6 -
                                                       0.0002984642) +
                                                0.02212276) +
                                         3.186519) +
                                        b_s * ((b_t * (b_t * 2.059331E-7 -
                                                       0.0001847318) +
                                                0.006704388) +
                                               s1o2 * 0.0001480266)) +
                                       b_p * ((b_t * (b_t * 1.39468E-8 -
                                                      1.202016E-6) +
                                               2.102898E-5) +
                                              b_s * (b_t * (b_t * 6.207323E-11 +
                                                            6.128773E-9) -
                                                     2.040237E-7))))) -
                    rn) -
                       d >
                   0.0);
            } else {
              ub = fb;
              if (dxm == 0.0) {
                fbpos = fapos;
                /*  since a = b = x */
              } else {
                /*  Evaluate difference between (a) eos at location on the cast
                 * where the */
                /*  pressure or depth is p, and (b) eos at location on the cast
                 * where the */
                /*  pressure or depth is p0 (where the surface currently is)
                 * plus the density */
                /*  perturbation d.  Part (b) is precomputed as r0.  Here, eos
                 * always */
                /*  evaluated at the pressure or depth of the original position,
                 * p0; this is */
                /*  to calculate locally referenced potential density with
                 * reference pressure */
                /*  p0. */
                /*  Interpolate S and T to the current pressure or depth */
                ppc_val2(Pn, Sppcn_data, Sppcn_size, Tppcn_data, fb, &b_s,
                         &b_t);
                /*  Calculate the potential density or potential specific volume
                 * difference */
                /* out =  eos(s, t, p0) - d - r0 ; */
                /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
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
                /*  Author(s)       : Geoff Stanley */
                /*  Email           : g.stanley@unsw.edu.au */
                /*  Email           : geoffstanley@gmail.com */
                /*  Version         : 1.0 */
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
                /*  4. secant bulk modulus K of sea water at p = 0 */
                /*  == original * 10 */
                /*  5. secant bulk modulus K of sea water at p */
                /*  == original / 10 */
                s1o2 = muDoubleScalarSqrt(b_s);
                fbpos =
                    ((((b_t * (b_t * (b_t * (b_t * (b_t * 6.536332E-9 -
                                                    1.120083E-6) +
                                             0.0001001685) -
                                      0.00909529) +
                               0.06793952) +
                        999.842594) +
                       b_s * (((b_t * (b_t * (b_t * (b_t * 5.3875E-9 -
                                                     8.2467E-7) +
                                              7.6438E-5) -
                                       0.0040899) +
                                0.824493) +
                               s1o2 * (b_t * (b_t * -1.6546E-6 + 0.00010227) -
                                       0.00572466)) +
                              b_s * 0.00048314)) /
                          (1.0 -
                           b_p /
                               (((b_t * (b_t * (b_t * (b_t * -0.0004190253 +
                                                       0.09648704) -
                                                17.06103) +
                                         1444.304) +
                                  196593.3) +
                                 b_s * ((b_t * (b_t * (b_t * -0.0005084188 +
                                                       0.06283263) -
                                                3.101089) +
                                         528.4855) +
                                        s1o2 * (b_t * (b_t * -0.004619924 +
                                                       0.09085835) +
                                                3.88664))) +
                                b_p * (((b_t * (b_t * (b_t * 1.956415E-6 -
                                                       0.0002984642) +
                                                0.02212276) +
                                         3.186519) +
                                        b_s * ((b_t * (b_t * 2.059331E-7 -
                                                       0.0001847318) +
                                                0.006704388) +
                                               s1o2 * 0.0001480266)) +
                                       b_p * ((b_t * (b_t * 1.39468E-8 -
                                                      1.202016E-6) +
                                               2.102898E-5) +
                                              b_s * (b_t * (b_t * 6.207323E-11 +
                                                            6.128773E-9) -
                                                     2.040237E-7))))) -
                      rn) -
                         d >
                     0.0);
              }
            }
            boolean_T guard1 = false;
            do {
              exitg1 = 0;
              guard1 = false;
              if (lb > tol) {
                /*  Move a left, and test for a sign change */
                dxm *= 1.4142135623730949;
                lb = muDoubleScalarMax(fb - dxm, tol);
                /*  Evaluate difference between (a) eos at location on the cast
                 * where the */
                /*  pressure or depth is p, and (b) eos at location on the cast
                 * where the */
                /*  pressure or depth is p0 (where the surface currently is)
                 * plus the density */
                /*  perturbation d.  Part (b) is precomputed as r0.  Here, eos
                 * always */
                /*  evaluated at the pressure or depth of the original position,
                 * p0; this is */
                /*  to calculate locally referenced potential density with
                 * reference pressure */
                /*  p0. */
                /*  Interpolate S and T to the current pressure or depth */
                ppc_val2(Pn, Sppcn_data, Sppcn_size, Tppcn_data, lb, &b_s,
                         &b_t);
                /*  Calculate the potential density or potential specific volume
                 * difference */
                /* out =  eos(s, t, p0) - d - r0 ; */
                /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
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
                /*  Author(s)       : Geoff Stanley */
                /*  Email           : g.stanley@unsw.edu.au */
                /*  Email           : geoffstanley@gmail.com */
                /*  Version         : 1.0 */
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
                /*  4. secant bulk modulus K of sea water at p = 0 */
                /*  == original * 10 */
                /*  5. secant bulk modulus K of sea water at p */
                /*  == original / 10 */
                s1o2 = muDoubleScalarSqrt(b_s);
                fapos =
                    ((((b_t * (b_t * (b_t * (b_t * (b_t * 6.536332E-9 -
                                                    1.120083E-6) +
                                             0.0001001685) -
                                      0.00909529) +
                               0.06793952) +
                        999.842594) +
                       b_s * (((b_t * (b_t * (b_t * (b_t * 5.3875E-9 -
                                                     8.2467E-7) +
                                              7.6438E-5) -
                                       0.0040899) +
                                0.824493) +
                               s1o2 * (b_t * (b_t * -1.6546E-6 + 0.00010227) -
                                       0.00572466)) +
                              b_s * 0.00048314)) /
                          (1.0 -
                           b_p /
                               (((b_t * (b_t * (b_t * (b_t * -0.0004190253 +
                                                       0.09648704) -
                                                17.06103) +
                                         1444.304) +
                                  196593.3) +
                                 b_s * ((b_t * (b_t * (b_t * -0.0005084188 +
                                                       0.06283263) -
                                                3.101089) +
                                         528.4855) +
                                        s1o2 * (b_t * (b_t * -0.004619924 +
                                                       0.09085835) +
                                                3.88664))) +
                                b_p * (((b_t * (b_t * (b_t * 1.956415E-6 -
                                                       0.0002984642) +
                                                0.02212276) +
                                         3.186519) +
                                        b_s * ((b_t * (b_t * 2.059331E-7 -
                                                       0.0001847318) +
                                                0.006704388) +
                                               s1o2 * 0.0001480266)) +
                                       b_p * ((b_t * (b_t * 1.39468E-8 -
                                                      1.202016E-6) +
                                               2.102898E-5) +
                                              b_s * (b_t * (b_t * 6.207323E-11 +
                                                            6.128773E-9) -
                                                     2.040237E-7))))) -
                      rn) -
                         d >
                     0.0);
                if ((boolean_T)(fapos ^ fbpos)) {
                  /*  fa and fb have different signs */
                  exitg1 = 1;
                } else {
                  guard1 = true;
                }
              } else if (ub == m) {
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
                if (ub < m) {
                  /*  Move b right, and test for a sign change */
                  dxp *= 1.4142135623730949;
                  ub = muDoubleScalarMin(fb + dxp, m);
                  /*  Evaluate difference between (a) eos at location on the
                   * cast where the */
                  /*  pressure or depth is p, and (b) eos at location on the
                   * cast where the */
                  /*  pressure or depth is p0 (where the surface currently is)
                   * plus the density */
                  /*  perturbation d.  Part (b) is precomputed as r0.  Here, eos
                   * always */
                  /*  evaluated at the pressure or depth of the original
                   * position, p0; this is */
                  /*  to calculate locally referenced potential density with
                   * reference pressure */
                  /*  p0. */
                  /*  Interpolate S and T to the current pressure or depth */
                  ppc_val2(Pn, Sppcn_data, Sppcn_size, Tppcn_data, ub, &b_s,
                           &b_t);
                  /*  Calculate the potential density or potential specific
                   * volume difference */
                  /* out =  eos(s, t, p0) - d - r0 ; */
                  /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density.
                   */
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
                  /*  Author(s)       : Geoff Stanley */
                  /*  Email           : g.stanley@unsw.edu.au */
                  /*  Email           : geoffstanley@gmail.com */
                  /*  Version         : 1.0 */
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
                  /*  4. secant bulk modulus K of sea water at p = 0 */
                  /*  == original * 10 */
                  /*  5. secant bulk modulus K of sea water at p */
                  /*  == original / 10 */
                  s1o2 = muDoubleScalarSqrt(b_s);
                  fbpos =
                      ((((b_t * (b_t * (b_t * (b_t * (b_t * 6.536332E-9 -
                                                      1.120083E-6) +
                                               0.0001001685) -
                                        0.00909529) +
                                 0.06793952) +
                          999.842594) +
                         b_s * (((b_t * (b_t * (b_t * (b_t * 5.3875E-9 -
                                                       8.2467E-7) +
                                                7.6438E-5) -
                                         0.0040899) +
                                  0.824493) +
                                 s1o2 * (b_t * (b_t * -1.6546E-6 + 0.00010227) -
                                         0.00572466)) +
                                b_s * 0.00048314)) /
                            (1.0 -
                             b_p /
                                 (((b_t * (b_t * (b_t * (b_t * -0.0004190253 +
                                                         0.09648704) -
                                                  17.06103) +
                                           1444.304) +
                                    196593.3) +
                                   b_s * ((b_t * (b_t * (b_t * -0.0005084188 +
                                                         0.06283263) -
                                                  3.101089) +
                                           528.4855) +
                                          s1o2 * (b_t * (b_t * -0.004619924 +
                                                         0.09085835) +
                                                  3.88664))) +
                                  b_p *
                                      (((b_t * (b_t * (b_t * 1.956415E-6 -
                                                       0.0002984642) +
                                                0.02212276) +
                                         3.186519) +
                                        b_s * ((b_t * (b_t * 2.059331E-7 -
                                                       0.0001847318) +
                                                0.006704388) +
                                               s1o2 * 0.0001480266)) +
                                       b_p * ((b_t * (b_t * 1.39468E-8 -
                                                      1.202016E-6) +
                                               2.102898E-5) +
                                              b_s * (b_t * (b_t * 6.207323E-11 +
                                                            6.128773E-9) -
                                                     2.040237E-7))))) -
                        rn) -
                           d >
                       0.0);
                  if ((boolean_T)(fapos ^ fbpos)) {
                    /*  fa and fb have different signs */
                    exitg1 = 1;
                  }
                } else if (lb == tol) {
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
        real_T b_d;
        real_T c;
        real_T e;
        real_T fc;
        /*  A sign change was discovered, so a root exists in the interval. */
        /*  Solve the nonlinear root-finding problem using Brent's method */
        /* p(n) = fzero_brent(@diff_avgx, lb, ub, tolp, ... */
        /*   Sppcn, Tppcn, Pn, s(n), t(n), p(n), d);  % DEV:  testing reference
         * p = average of current and point-to-be */
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
        /*  Evaluate difference between (a) eos at location on the cast where
         * the */
        /*  pressure or depth is p, and (b) eos at location on the cast where
         * the */
        /*  pressure or depth is p0 (where the surface currently is) plus the
         * density */
        /*  perturbation d.  Part (b) is precomputed as r0.  Here, eos always */
        /*  evaluated at the pressure or depth of the original position, p0;
         * this is */
        /*  to calculate locally referenced potential density with reference
         * pressure */
        /*  p0. */
        /*  Interpolate S and T to the current pressure or depth */
        ppc_val2(Pn, Sppcn_data, Sppcn_size, Tppcn_data, lb, &b_s, &b_t);
        /*  Calculate the potential density or potential specific volume
         * difference */
        /* out =  eos(s, t, p0) - d - r0 ; */
        /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
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
        /*  Author(s)       : Geoff Stanley */
        /*  Email           : g.stanley@unsw.edu.au */
        /*  Email           : geoffstanley@gmail.com */
        /*  Version         : 1.0 */
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
        b_p = p_data[n] * 1.015335;
        /*  Henceforth z is actually pressure [dbar] */
        /*  coefficients nonlinear equation of state in pressure coordinates for
         */
        /*  1. density of fresh water at p = 0 */
        /*  2. density of sea water at p = 0 */
        /*  coefficients in pressure coordinates for */
        /*  3. secant bulk modulus K of fresh water at p = 0 */
        /*  == original * 10 */
        /*  4. secant bulk modulus K of sea water at p = 0 */
        /*  == original * 10 */
        /*  5. secant bulk modulus K of sea water at p */
        /*  == original / 10 */
        s1o2 = muDoubleScalarSqrt(b_s);
        dxp =
            (((b_t * (b_t * (b_t * (b_t * (b_t * 6.536332E-9 - 1.120083E-6) +
                                    0.0001001685) -
                             0.00909529) +
                      0.06793952) +
               999.842594) +
              b_s * (((b_t * (b_t * (b_t * (b_t * 5.3875E-9 - 8.2467E-7) +
                                     7.6438E-5) -
                              0.0040899) +
                       0.824493) +
                      s1o2 * (b_t * (b_t * -1.6546E-6 + 0.00010227) -
                              0.00572466)) +
                     b_s * 0.00048314)) /
                 (1.0 -
                  b_p /
                      (((b_t *
                             (b_t * (b_t * (b_t * -0.0004190253 + 0.09648704) -
                                     17.06103) +
                              1444.304) +
                         196593.3) +
                        b_s *
                            ((b_t * (b_t * (b_t * -0.0005084188 + 0.06283263) -
                                     3.101089) +
                              528.4855) +
                             s1o2 * (b_t * (b_t * -0.004619924 + 0.09085835) +
                                     3.88664))) +
                       b_p *
                           (((b_t * (b_t * (b_t * 1.956415E-6 - 0.0002984642) +
                                     0.02212276) +
                              3.186519) +
                             b_s * ((b_t * (b_t * 2.059331E-7 - 0.0001847318) +
                                     0.006704388) +
                                    s1o2 * 0.0001480266)) +
                            b_p * ((b_t * (b_t * 1.39468E-8 - 1.202016E-6) +
                                    2.102898E-5) +
                                   b_s * (b_t * (b_t * 6.207323E-11 +
                                                 6.128773E-9) -
                                          2.040237E-7))))) -
             rn) -
            phi_data[n];
        /*  Evaluate difference between (a) eos at location on the cast where
         * the */
        /*  pressure or depth is p, and (b) eos at location on the cast where
         * the */
        /*  pressure or depth is p0 (where the surface currently is) plus the
         * density */
        /*  perturbation d.  Part (b) is precomputed as r0.  Here, eos always */
        /*  evaluated at the pressure or depth of the original position, p0;
         * this is */
        /*  to calculate locally referenced potential density with reference
         * pressure */
        /*  p0. */
        /*  Interpolate S and T to the current pressure or depth */
        ppc_val2(Pn, Sppcn_data, Sppcn_size, Tppcn_data, ub, &b_s, &b_t);
        /*  Calculate the potential density or potential specific volume
         * difference */
        /* out =  eos(s, t, p0) - d - r0 ; */
        /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
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
        /*  Author(s)       : Geoff Stanley */
        /*  Email           : g.stanley@unsw.edu.au */
        /*  Email           : geoffstanley@gmail.com */
        /*  Version         : 1.0 */
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
        b_p = p_data[n] * 1.015335;
        /*  Henceforth z is actually pressure [dbar] */
        /*  coefficients nonlinear equation of state in pressure coordinates for
         */
        /*  1. density of fresh water at p = 0 */
        /*  2. density of sea water at p = 0 */
        /*  coefficients in pressure coordinates for */
        /*  3. secant bulk modulus K of fresh water at p = 0 */
        /*  == original * 10 */
        /*  4. secant bulk modulus K of sea water at p = 0 */
        /*  == original * 10 */
        /*  5. secant bulk modulus K of sea water at p */
        /*  == original / 10 */
        s1o2 = muDoubleScalarSqrt(b_s);
        fb = (((b_t * (b_t * (b_t * (b_t * (b_t * 6.536332E-9 - 1.120083E-6) +
                                     0.0001001685) -
                              0.00909529) +
                       0.06793952) +
                999.842594) +
               b_s * (((b_t * (b_t * (b_t * (b_t * 5.3875E-9 - 8.2467E-7) +
                                      7.6438E-5) -
                               0.0040899) +
                        0.824493) +
                       s1o2 * (b_t * (b_t * -1.6546E-6 + 0.00010227) -
                               0.00572466)) +
                      b_s * 0.00048314)) /
                  (1.0 -
                   b_p /
                       (((b_t *
                              (b_t * (b_t * (b_t * -0.0004190253 + 0.09648704) -
                                      17.06103) +
                               1444.304) +
                          196593.3) +
                         b_s *
                             ((b_t * (b_t * (b_t * -0.0005084188 + 0.06283263) -
                                      3.101089) +
                               528.4855) +
                              s1o2 * (b_t * (b_t * -0.004619924 + 0.09085835) +
                                      3.88664))) +
                        b_p *
                            (((b_t * (b_t * (b_t * 1.956415E-6 - 0.0002984642) +
                                      0.02212276) +
                               3.186519) +
                              b_s * ((b_t * (b_t * 2.059331E-7 - 0.0001847318) +
                                      0.006704388) +
                                     s1o2 * 0.0001480266)) +
                             b_p * ((b_t * (b_t * 1.39468E-8 - 1.202016E-6) +
                                     2.102898E-5) +
                                    b_s * (b_t * (b_t * 6.207323E-11 +
                                                  6.128773E-9) -
                                           2.040237E-7))))) -
              rn) -
             phi_data[n];
        c = lb;
        fc = dxp;
        e = ub - lb;
        b_d = e;
        do {
          exitg1 = 0;
          if (muDoubleScalarAbs(fc) < muDoubleScalarAbs(fb)) {
            lb = ub;
            ub = c;
            c = lb;
            dxp = fb;
            fb = fc;
            fc = dxp;
          }
          tol = 4.4408920985006262E-16 * muDoubleScalarAbs(ub) + tolp;
          m = 0.5 * (c - ub);
          if ((muDoubleScalarAbs(m) <= tol) || (fb == 0.0)) {
            exitg1 = 1;
          } else {
            if ((muDoubleScalarAbs(e) < tol) ||
                (muDoubleScalarAbs(dxp) <= muDoubleScalarAbs(fb))) {
              e = m;
              b_d = m;
            } else {
              b_s = fb / dxp;
              if (lb == c) {
                b_p = 2.0 * m * b_s;
                dxp = 1.0 - b_s;
              } else {
                dxp /= fc;
                dxm = fb / fc;
                b_p = b_s *
                      (2.0 * m * dxp * (dxp - dxm) - (ub - lb) * (dxm - 1.0));
                dxp = (dxp - 1.0) * (dxm - 1.0) * (b_s - 1.0);
              }
              if (b_p > 0.0) {
                dxp = -dxp;
              } else {
                b_p = -b_p;
              }
              b_s = e;
              e = b_d;
              if ((2.0 * b_p < 3.0 * m * dxp - muDoubleScalarAbs(tol * dxp)) &&
                  (b_p < muDoubleScalarAbs(0.5 * b_s * dxp))) {
                b_d = b_p / dxp;
              } else {
                e = m;
                b_d = m;
              }
            }
            lb = ub;
            dxp = fb;
            if (tol < muDoubleScalarAbs(b_d)) {
              ub += b_d;
            } else if (m > 0.0) {
              ub += tol;
            } else {
              ub -= tol;
            }
            /*  Evaluate difference between (a) eos at location on the cast
             * where the */
            /*  pressure or depth is p, and (b) eos at location on the cast
             * where the */
            /*  pressure or depth is p0 (where the surface currently is) plus
             * the density */
            /*  perturbation d.  Part (b) is precomputed as r0.  Here, eos
             * always */
            /*  evaluated at the pressure or depth of the original position, p0;
             * this is */
            /*  to calculate locally referenced potential density with reference
             * pressure */
            /*  p0. */
            /*  Interpolate S and T to the current pressure or depth */
            ppc_val2(Pn, Sppcn_data, Sppcn_size, Tppcn_data, ub, &b_s, &b_t);
            /*  Calculate the potential density or potential specific volume
             * difference */
            /* out =  eos(s, t, p0) - d - r0 ; */
            /* EOSCG_DENSJMD95_BSQ Fast Boussinesq JMD95 in-situ density. */
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
            /*  Author(s)       : Geoff Stanley */
            /*  Email           : g.stanley@unsw.edu.au */
            /*  Email           : geoffstanley@gmail.com */
            /*  Version         : 1.0 */
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
            b_p = p_data[n] * 1.015335;
            /*  Henceforth z is actually pressure [dbar] */
            /*  coefficients nonlinear equation of state in pressure coordinates
             * for */
            /*  1. density of fresh water at p = 0 */
            /*  2. density of sea water at p = 0 */
            /*  coefficients in pressure coordinates for */
            /*  3. secant bulk modulus K of fresh water at p = 0 */
            /*  == original * 10 */
            /*  4. secant bulk modulus K of sea water at p = 0 */
            /*  == original * 10 */
            /*  5. secant bulk modulus K of sea water at p */
            /*  == original / 10 */
            s1o2 = muDoubleScalarSqrt(b_s);
            fb = (((b_t *
                        (b_t * (b_t * (b_t * (b_t * 6.536332E-9 - 1.120083E-6) +
                                       0.0001001685) -
                                0.00909529) +
                         0.06793952) +
                    999.842594) +
                   b_s * (((b_t * (b_t * (b_t * (b_t * 5.3875E-9 - 8.2467E-7) +
                                          7.6438E-5) -
                                   0.0040899) +
                            0.824493) +
                           s1o2 * (b_t * (b_t * -1.6546E-6 + 0.00010227) -
                                   0.00572466)) +
                          b_s * 0.00048314)) /
                      (1.0 -
                       b_p / (((b_t * (b_t * (b_t * (b_t * -0.0004190253 +
                                                     0.09648704) -
                                              17.06103) +
                                       1444.304) +
                                196593.3) +
                               b_s * ((b_t * (b_t * (b_t * -0.0005084188 +
                                                     0.06283263) -
                                              3.101089) +
                                       528.4855) +
                                      s1o2 * (b_t * (b_t * -0.004619924 +
                                                     0.09085835) +
                                              3.88664))) +
                              b_p * (((b_t * (b_t * (b_t * 1.956415E-6 -
                                                     0.0002984642) +
                                              0.02212276) +
                                       3.186519) +
                                      b_s * ((b_t * (b_t * 2.059331E-7 -
                                                     0.0001847318) +
                                              0.006704388) +
                                             s1o2 * 0.0001480266)) +
                                     b_p * ((b_t * (b_t * 1.39468E-8 -
                                                    1.202016E-6) +
                                             2.102898E-5) +
                                            b_s * (b_t * (b_t * 6.207323E-11 +
                                                          6.128773E-9) -
                                                   2.040237E-7))))) -
                  rn) -
                 d;
            if (((fb > 0.0) && (fc > 0.0)) || ((fb <= 0.0) && (fc <= 0.0))) {
              c = lb;
              fc = dxp;
              e = ub - lb;
              b_d = e;
            }
          }
        } while (exitg1 == 0);
        p_data[n] = ub;
        /*  Interpolate S and T onto the updated surface */
        ppc_val2(Pn, Sppcn_data, Sppcn_size, Tppcn_data, p_data[n], &d, &b_p);
        t_data[n] = b_p;
        s_data[n] = d;
      } else {
        p_data[n] = rtNaN;
        s_data[n] = rtNaN;
        t_data[n] = rtNaN;
      }
    } else {
      /*  phi is nan, or only one grid cell so cannot interpolate. */
      /*  This will ensure s,t,p all have the same nan structure */
      p_data[n] = rtNaN;
      s_data[n] = rtNaN;
      t_data[n] = rtNaN;
    }
  }
  emxFree_real_T(&Pn);
  emlrtHeapReferenceStackLeaveFcnR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (omega_vertsolve.c) */
