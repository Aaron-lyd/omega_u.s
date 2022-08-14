/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * ppc_val.c
 *
 * Code generation for function 'ppc_val'
 *
 */

/* Include files */
#include "ppc_val.h"
#include "assertValidSizeArg.h"
#include "ppc_val_data.h"
#include "ppc_val_emxutil.h"
#include "ppc_val_types.h"
#include "prod.h"
#include "rt_nonfinite.h"
#include "squeeze.h"
#include "mwmathutil.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = {
    77,        /* lineNo */
    "ppc_val", /* fcnName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/omega_u.s/lib/ppc/ppc_val.m" /* pathName */
};

static emlrtRSInfo b_emlrtRSI = {
    112,       /* lineNo */
    "ppc_val", /* fcnName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/omega_u.s/lib/ppc/ppc_val.m" /* pathName */
};

static emlrtRSInfo c_emlrtRSI = {
    126,       /* lineNo */
    "ppc_val", /* fcnName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/omega_u.s/lib/ppc/ppc_val.m" /* pathName */
};

static emlrtRSInfo d_emlrtRSI = {
    140,       /* lineNo */
    "ppc_val", /* fcnName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/omega_u.s/lib/ppc/ppc_val.m" /* pathName */
};

static emlrtRSInfo e_emlrtRSI = {
    29,                  /* lineNo */
    "reshapeSizeChecks", /* fcnName */
    "/Applications/MATLAB_R2021a.app/toolbox/eml/eml/+coder/+internal/"
    "reshapeSizeChecks.m" /* pathName */
};

static emlrtRSInfo f_emlrtRSI = {
    109,               /* lineNo */
    "computeDimsData", /* fcnName */
    "/Applications/MATLAB_R2021a.app/toolbox/eml/eml/+coder/+internal/"
    "reshapeSizeChecks.m" /* pathName */
};

static emlrtRSInfo l_emlrtRSI = {
    28,      /* lineNo */
    "colon", /* fcnName */
    "/Applications/MATLAB_R2021a.app/toolbox/eml/lib/matlab/ops/colon.m" /* pathName
                                                                          */
};

static emlrtRSInfo m_emlrtRSI = {
    117,     /* lineNo */
    "colon", /* fcnName */
    "/Applications/MATLAB_R2021a.app/toolbox/eml/lib/matlab/ops/colon.m" /* pathName
                                                                          */
};

static emlrtRSInfo n_emlrtRSI = {
    311,               /* lineNo */
    "eml_float_colon", /* fcnName */
    "/Applications/MATLAB_R2021a.app/toolbox/eml/lib/matlab/ops/colon.m" /* pathName
                                                                          */
};

static emlrtRSInfo o_emlrtRSI = {
    18,                  /* lineNo */
    "reshapeSizeChecks", /* fcnName */
    "/Applications/MATLAB_R2021a.app/toolbox/eml/eml/+coder/+internal/"
    "reshapeSizeChecks.m" /* pathName */
};

static emlrtMCInfo emlrtMCI = {
    73,        /* lineNo */
    1,         /* colNo */
    "ppc_val", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/omega_u.s/lib/ppc/ppc_val.m" /* pName */
};

static emlrtBCInfo emlrtBCI = {
    -1,        /* iFirst */
    -1,        /* iLast */
    127,       /* lineNo */
    21,        /* colNo */
    "y",       /* aName */
    "ppc_val", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/omega_u.s/lib/ppc/ppc_val.m", /* pName */
    0                                             /* checkKind */
};

static emlrtBCInfo b_emlrtBCI = {
    -1,        /* iFirst */
    -1,        /* iLast */
    127,       /* lineNo */
    41,        /* colNo */
    "C",       /* aName */
    "ppc_val", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/omega_u.s/lib/ppc/ppc_val.m", /* pName */
    0                                             /* checkKind */
};

static emlrtBCInfo c_emlrtBCI = {
    -1,        /* iFirst */
    -1,        /* iLast */
    127,       /* lineNo */
    33,        /* colNo */
    "y",       /* aName */
    "ppc_val", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/omega_u.s/lib/ppc/ppc_val.m", /* pName */
    0                                             /* checkKind */
};

static emlrtBCInfo d_emlrtBCI = {
    -1,        /* iFirst */
    -1,        /* iLast */
    113,       /* lineNo */
    17,        /* colNo */
    "y",       /* aName */
    "ppc_val", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/omega_u.s/lib/ppc/ppc_val.m", /* pName */
    0                                             /* checkKind */
};

static emlrtDCInfo emlrtDCI = {
    113,       /* lineNo */
    25,        /* colNo */
    "ppc_val", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/omega_u.s/lib/ppc/ppc_val.m", /* pName */
    1                                             /* checkKind */
};

static emlrtBCInfo e_emlrtBCI = {
    -1,        /* iFirst */
    -1,        /* iLast */
    113,       /* lineNo */
    25,        /* colNo */
    "C",       /* aName */
    "ppc_val", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/omega_u.s/lib/ppc/ppc_val.m", /* pName */
    0                                             /* checkKind */
};

static emlrtBCInfo f_emlrtBCI = {
    -1,        /* iFirst */
    -1,        /* iLast */
    123,       /* lineNo */
    17,        /* colNo */
    "y",       /* aName */
    "ppc_val", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/omega_u.s/lib/ppc/ppc_val.m", /* pName */
    0                                             /* checkKind */
};

static emlrtBCInfo g_emlrtBCI = {
    -1,        /* iFirst */
    -1,        /* iLast */
    89,        /* lineNo */
    15,        /* colNo */
    "x",       /* aName */
    "ppc_val", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/omega_u.s/lib/ppc/ppc_val.m", /* pName */
    0                                             /* checkKind */
};

static emlrtRTEInfo emlrtRTEI = {
    116,               /* lineNo */
    9,                 /* colNo */
    "computeDimsData", /* fName */
    "/Applications/MATLAB_R2021a.app/toolbox/eml/eml/+coder/+internal/"
    "reshapeSizeChecks.m" /* pName */
};

static emlrtRTEInfo b_emlrtRTEI = {
    125,       /* lineNo */
    25,        /* colNo */
    "ppc_val", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/omega_u.s/lib/ppc/ppc_val.m" /* pName */
};

static emlrtRTEInfo e_emlrtRTEI = {
    417,               /* lineNo */
    15,                /* colNo */
    "assert_pmaxsize", /* fName */
    "/Applications/MATLAB_R2021a.app/toolbox/eml/lib/matlab/ops/colon.m" /* pName
                                                                          */
};

static emlrtRTEInfo j_emlrtRTEI = {
    78,        /* lineNo */
    1,         /* colNo */
    "ppc_val", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/omega_u.s/lib/ppc/ppc_val.m" /* pName */
};

static emlrtRTEInfo k_emlrtRTEI = {
    28,      /* lineNo */
    9,       /* colNo */
    "colon", /* fName */
    "/Applications/MATLAB_R2021a.app/toolbox/eml/lib/matlab/ops/colon.m" /* pName
                                                                          */
};

static emlrtRTEInfo l_emlrtRTEI = {
    312,     /* lineNo */
    20,      /* colNo */
    "colon", /* fName */
    "/Applications/MATLAB_R2021a.app/toolbox/eml/lib/matlab/ops/colon.m" /* pName
                                                                          */
};

static emlrtRTEInfo m_emlrtRTEI = {
    1,         /* lineNo */
    10,        /* colNo */
    "ppc_val", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/omega_u.s/lib/ppc/ppc_val.m" /* pName */
};

static emlrtRTEInfo n_emlrtRTEI = {
    126,       /* lineNo */
    30,        /* colNo */
    "ppc_val", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/omega_u.s/lib/ppc/ppc_val.m" /* pName */
};

static emlrtRSInfo q_emlrtRSI = {
    73,        /* lineNo */
    "ppc_val", /* fcnName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/omega_u.s/lib/ppc/ppc_val.m" /* pathName */
};

static emlrtRSInfo r_emlrtRSI = {
    18,            /* lineNo */
    "indexDivide", /* fcnName */
    "/Applications/MATLAB_R2021a.app/toolbox/eml/eml/+coder/+internal/"
    "indexDivide.m" /* pathName */
};

/* Function Declarations */
static int32_T div_s32(const emlrtStack *sp, int32_T numerator,
                       int32_T denominator);

static void error(const emlrtStack *sp, const mxArray *b, const mxArray *c,
                  emlrtMCInfo *location);

/* Function Definitions */
static int32_T div_s32(const emlrtStack *sp, int32_T numerator,
                       int32_T denominator)
{
  int32_T quotient;
  uint32_T b_denominator;
  uint32_T b_numerator;
  if (denominator == 0) {
    emlrtDivisionByZeroErrorR2012b(NULL, (emlrtCTX)sp);
  } else {
    if (numerator < 0) {
      b_numerator = ~(uint32_T)numerator + 1U;
    } else {
      b_numerator = (uint32_T)numerator;
    }
    if (denominator < 0) {
      b_denominator = ~(uint32_T)denominator + 1U;
    } else {
      b_denominator = (uint32_T)denominator;
    }
    b_numerator /= b_denominator;
    if ((numerator < 0) != (denominator < 0)) {
      quotient = -(int32_T)b_numerator;
    } else {
      quotient = (int32_T)b_numerator;
    }
  }
  return quotient;
}

static void error(const emlrtStack *sp, const mxArray *b, const mxArray *c,
                  emlrtMCInfo *location)
{
  const mxArray *pArrays[2];
  pArrays[0] = b;
  pArrays[1] = c;
  emlrtCallMATLABR2012b((emlrtCTX)sp, 0, NULL, 2, &pArrays[0],
                        (const char_T *)"error", true, location);
}

void ppc_val(const emlrtStack *sp, const real_T X[20], const emxArray_real_T *C,
             const real_T x_data[], const int32_T x_size[2], real_T d,
             emxArray_real_T *y)
{
  static const int32_T iv[2] = {1, 34};
  static const char_T u[34] = {'x', ' ', 'm', 'u', 's', 't', ' ', 'h', 'a',
                               'v', 'e', ' ', '%', 'd', ' ', 'c', 'o', 'l',
                               'u', 'm', 'n', 's', ' ', 'o', 'r', ' ', '1',
                               ' ', 'c', 'o', 'l', 'u', 'm', 'n'};
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack st;
  emxArray_real_T f_y;
  emxArray_real_T *d_y;
  emxArray_real_T *e_y;
  const mxArray *b_y;
  const mxArray *c_y;
  const mxArray *m;
  real_T szy[4];
  real_T apnd;
  real_T cdiff;
  real_T ndbl;
  real_T p;
  real_T t;
  real_T xln;
  int32_T num[4];
  int32_T L;
  int32_T MC;
  int32_T O;
  int32_T b_i;
  int32_T calclen;
  int32_T i;
  int32_T i1;
  int32_T k;
  int32_T l;
  int32_T maxdimlen;
  int32_T n;
  int32_T nx;
  int32_T o;
  int32_T xM;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtCTX)sp);
  /* PPC_VAL  Piecewise Polynomial Evaluation, quick */
  /*  */
  /*  */
  /*  y = ppc_val(X, C, x) */
  /*  evaluates the piecewise polynomials whose coefficients are C and whose */
  /*  knots are X, at data sites x.  */
  /*  */
  /*  y = ppc_val(..., d) */
  /*  as above but evaluates the d'th derivative. */
  /*  */
  /*  */
  /*  --- Input: */
  /*  X [    K   x M], knots of the piecewise polynomials */
  /*  C [O x K-1 x N], coefficients of the first piecewise polynomial */
  /*  x [    L   x M], evaluation sites */
  /*  d [1, 1]       , order of derivative to evaluate; defaults to 0 if not */
  /*                   given, which means to evaluate with no differentiation.
   */
  /*  */
  /*  */
  /*  --- Output: */
  /*  y [L x N], the first piecewise polynomial evaluated at x */
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
  /*  If L == 1, this dimension is squeeze()'d out of y. */
  /*  */
  /*  */
  /*  --- Acknowledgements: */
  /*  This code is adapted from MATLAB's ppval.m */
  /*  Author    : Geoff Stanley */
  /*  Email     : geoffstanley@gmail.com */
  /*  Version   : 1.0 */
  /*  History   : 24/10/2019 - initial release */
  /*  Set default to evaluate, not differentiate: */
  O = (int8_T)C->size[0];
  /*  Order of the piecewise polynomial */
  /*  number of knots of the piecewise polynomials */
  L = (int8_T)x_size[0];
  /*  number of levels to interpolate */
  xM = (int8_T)x_size[1];
  /*  Handle case when L should be 1. */
  if ((int8_T)x_size[0] * (int8_T)x_size[1] == 4225) {
    L = 1;
    xM = 4225;
  }
  szy[0] = L;
  szy[1] = 65.0;
  szy[2] = 65.0;
  szy[3] = 1.0;
  /*  Add a trailing 1, to ensure 2 dimensions at least */
  if ((xM != 1) && (xM != 4225)) {
    b_y = NULL;
    m = emlrtCreateCharArray(2, &iv[0]);
    emlrtInitCharArrayR2013a((emlrtCTX)sp, 34, m, &u[0]);
    emlrtAssign(&b_y, m);
    c_y = NULL;
    m = emlrtCreateDoubleScalar(4225.0);
    emlrtAssign(&c_y, m);
    st.site = &q_emlrtRSI;
    error(&st, b_y, c_y, &emlrtMCI);
  }
  st.site = &emlrtRSI;
  nx = x_size[0] * x_size[1];
  b_st.site = &e_emlrtRSI;
  c_st.site = &f_emlrtRSI;
  if (L > 0) {
    c_st.site = &r_emlrtRSI;
    calclen = div_s32(&c_st, nx, L);
    if (calclen > nx) {
      emlrtErrorWithMessageIdR2018a(&b_st, &emlrtRTEI,
                                    "Coder:builtins:AssertionFailed",
                                    "Coder:builtins:AssertionFailed", 0);
    }
  } else {
    calclen = 0;
  }
  n = x_size[0];
  if (x_size[1] > x_size[0]) {
    n = x_size[1];
  }
  maxdimlen = muIntScalarMax_sint32(nx, n);
  if (L > maxdimlen) {
    emlrtErrorWithMessageIdR2018a(&st, &c_emlrtRTEI,
                                  "Coder:toolbox:reshape_emptyReshapeLimit",
                                  "Coder:toolbox:reshape_emptyReshapeLimit", 0);
  }
  if (calclen > maxdimlen) {
    emlrtErrorWithMessageIdR2018a(&st, &c_emlrtRTEI,
                                  "Coder:toolbox:reshape_emptyReshapeLimit",
                                  "Coder:toolbox:reshape_emptyReshapeLimit", 0);
  }
  i = L * calclen;
  if (i != nx) {
    emlrtErrorWithMessageIdR2018a(
        &st, &d_emlrtRTEI, "Coder:MATLAB:getReshapeDims_notSameNumel",
        "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
  }
  emxInit_real_T(&st, &d_y, 2, &m_emlrtRTEI, true);
  i1 = d_y->size[0] * d_y->size[1];
  d_y->size[0] = L;
  d_y->size[1] = 4225;
  emxEnsureCapacity_real_T(sp, d_y, i1, &j_emlrtRTEI);
  maxdimlen = L * 4225;
  for (i1 = 0; i1 < maxdimlen; i1++) {
    d_y->data[i1] = rtNaN;
  }
  /*  Evaluate each piecewise polynomial */
  xM = (xM > 1);
  /*  used for linear indexing */
  /*  used for linear indexing */
  MC = (int8_T)C->size[0] * 19;
  /*  used for linear indexing */
  /*  used for linear indexing */
  emxInit_real_T(sp, &e_y, 2, &n_emlrtRTEI, true);
  for (n = 0; n < 4225; n++) {
    for (l = 0; l < L; l++) {
      i1 = n * L;
      k = (l + i1 * xM) + 1;
      if (k > i) {
        emlrtDynamicBoundsCheckR2012b(k, 1, i, &g_emlrtBCI, (emlrtCTX)sp);
      }
      xln = x_data[k - 1];
      /*  x(l,n) or x(l,1) as appropriate */
      if ((!muDoubleScalarIsNaN(xln)) && (!muDoubleScalarIsNaN(X[0])) &&
          (!(xln < X[0])) && (!(xln > X[19]))) {
        /*  Leftmost binary search to find i such that: */
        /*  i = 1                      if xln <= X(1), or */
        /*  i = M                      if X(M) < xln */
        /*  X(i-1,n) < xln <= X(i,n)   otherwise */
        /*  We use leftmost so that NaN's at the end of X are treated as though
         * they are Inf. */
        b_i = 1;
        /*  Result will be >= 1 always */
        k = 20;
        /*  Result will be <= M always */
        while (b_i < k) {
          maxdimlen = (int32_T)muDoubleScalarFloor((real_T)(b_i + k) / 2.0);
          if (X[maxdimlen - 1] < xln) {
            /*  [X(j,n) or X(j,1) as appropriate]  < xln */
            b_i = maxdimlen + 1;
          } else {
            k = maxdimlen;
          }
          if (*emlrtBreakCheckR2012bFlagVar != 0) {
            emlrtBreakCheckR2012b((emlrtCTX)sp);
          }
        }
        nx = (l + i1) + 1;
        if (b_i == 1) {
          /*  Note: X(nX + 1) == xln   is guaranteed */
          if (muDoubleScalarIsNaN(d)) {
            i1 = e_y->size[0] * e_y->size[1];
            e_y->size[0] = 1;
            e_y->size[1] = 1;
            emxEnsureCapacity_real_T(sp, e_y, i1, &k_emlrtRTEI);
            e_y->data[0] = rtNaN;
          } else if (d < 1.0) {
            e_y->size[0] = 1;
            e_y->size[1] = 0;
          } else if (muDoubleScalarIsInf(d) && (1.0 == d)) {
            i1 = e_y->size[0] * e_y->size[1];
            e_y->size[0] = 1;
            e_y->size[1] = 1;
            emxEnsureCapacity_real_T(sp, e_y, i1, &k_emlrtRTEI);
            e_y->data[0] = rtNaN;
          } else {
            i1 = e_y->size[0] * e_y->size[1];
            e_y->size[0] = 1;
            maxdimlen = (int32_T)muDoubleScalarFloor(d - 1.0);
            e_y->size[1] = maxdimlen + 1;
            emxEnsureCapacity_real_T(sp, e_y, i1, &k_emlrtRTEI);
            for (i1 = 0; i1 <= maxdimlen; i1++) {
              e_y->data[i1] = (real_T)i1 + 1.0;
            }
          }
          st.site = &b_emlrtRSI;
          p = prod(&st, e_y);
          /*  Build integer multiplying the coefficient */
          i1 = d_y->size[0] * 4225;
          if (nx > i1) {
            emlrtDynamicBoundsCheckR2012b(nx, 1, i1, &d_emlrtBCI, (emlrtCTX)sp);
          }
          xln = ((real_T)O - d) + (real_T)(n * MC);
          if (xln != (int32_T)muDoubleScalarFloor(xln)) {
            emlrtIntegerCheckR2012b(xln, &emlrtDCI, (emlrtCTX)sp);
          }
          i1 = C->size[0] * 19 * 65 * 65;
          if (((int32_T)xln < 1) || ((int32_T)xln > i1)) {
            emlrtDynamicBoundsCheckR2012b((int32_T)xln, 1, i1, &e_emlrtBCI,
                                          (emlrtCTX)sp);
          }
          d_y->data[nx - 1] = C->data[(int32_T)xln - 1] * p;
          /*  y(l,n) = C(O - d, i, n) * p; */
        } else {
          /*  Evaluate this piece of the polynomial (see ppval.m) */
          t = xln - X[b_i - 2];
          /*  Switch to local coordinates */
          /*  Overload variable i, to speed with indexing */
          /*  subtract 1 from i so that 1 <= i <= M-1, and X(i) <= xln < X(i+1),
           */
          /*  subtract another 1 for indexing.  */
          b_i = (b_i - 2) * O + n * MC;
          i1 = d_y->size[0] * 4225;
          if (nx > i1) {
            emlrtDynamicBoundsCheckR2012b(nx, 1, i1, &f_emlrtBCI, (emlrtCTX)sp);
          }
          d_y->data[nx - 1] = 0.0;
          xln = (real_T)O - d;
          i1 = (int32_T)xln;
          emlrtForLoopVectorCheckR2021a(1.0, 1.0, xln, mxDOUBLE_CLASS,
                                        (int32_T)xln, &b_emlrtRTEI,
                                        (emlrtCTX)sp);
          for (o = 0; o < i1; o++) {
            st.site = &c_emlrtRSI;
            xln = (real_T)O - ((real_T)o + 1.0);
            p = (xln - d) + 1.0;
            b_st.site = &l_emlrtRSI;
            if (muDoubleScalarIsNaN(p)) {
              k = e_y->size[0] * e_y->size[1];
              e_y->size[0] = 1;
              e_y->size[1] = 1;
              emxEnsureCapacity_real_T(&b_st, e_y, k, &k_emlrtRTEI);
              e_y->data[0] = rtNaN;
            } else if ((int32_T)xln < p) {
              e_y->size[0] = 1;
              e_y->size[1] = 0;
            } else if (muDoubleScalarIsInf(p) && (p == (int32_T)xln)) {
              k = e_y->size[0] * e_y->size[1];
              e_y->size[0] = 1;
              e_y->size[1] = 1;
              emxEnsureCapacity_real_T(&b_st, e_y, k, &k_emlrtRTEI);
              e_y->data[0] = rtNaN;
            } else if (muDoubleScalarFloor(p) == p) {
              k = e_y->size[0] * e_y->size[1];
              e_y->size[0] = 1;
              maxdimlen = (int32_T)((real_T)(int32_T)xln - p);
              e_y->size[1] = maxdimlen + 1;
              emxEnsureCapacity_real_T(&b_st, e_y, k, &k_emlrtRTEI);
              for (k = 0; k <= maxdimlen; k++) {
                e_y->data[k] = p + (real_T)k;
              }
            } else {
              c_st.site = &m_emlrtRSI;
              ndbl = muDoubleScalarFloor(((real_T)(int32_T)xln - p) + 0.5);
              apnd = p + ndbl;
              cdiff = apnd - (real_T)(int32_T)xln;
              if (muDoubleScalarAbs(cdiff) <
                  4.4408920985006262E-16 *
                      muDoubleScalarMax(muDoubleScalarAbs(p),
                                        muDoubleScalarAbs((int32_T)xln))) {
                ndbl++;
                apnd = (int32_T)xln;
              } else if (cdiff > 0.0) {
                apnd = p + (ndbl - 1.0);
              } else {
                ndbl++;
              }
              if (ndbl >= 0.0) {
                calclen = (int32_T)ndbl;
              } else {
                calclen = 0;
              }
              d_st.site = &n_emlrtRSI;
              if (ndbl > 2.147483647E+9) {
                emlrtErrorWithMessageIdR2018a(&d_st, &e_emlrtRTEI,
                                              "Coder:MATLAB:pmaxsize",
                                              "Coder:MATLAB:pmaxsize", 0);
              }
              k = e_y->size[0] * e_y->size[1];
              e_y->size[0] = 1;
              e_y->size[1] = calclen;
              emxEnsureCapacity_real_T(&c_st, e_y, k, &l_emlrtRTEI);
              if (calclen > 0) {
                e_y->data[0] = p;
                if (calclen > 1) {
                  e_y->data[calclen - 1] = apnd;
                  maxdimlen = (calclen - 1) / 2;
                  for (k = 0; k <= maxdimlen - 2; k++) {
                    e_y->data[k + 1] = p + ((real_T)k + 1.0);
                    e_y->data[(calclen - k) - 2] = apnd - ((real_T)k + 1.0);
                  }
                  if (maxdimlen << 1 == calclen - 1) {
                    e_y->data[maxdimlen] = (p + apnd) / 2.0;
                  } else {
                    e_y->data[maxdimlen] = p + (real_T)maxdimlen;
                    e_y->data[maxdimlen + 1] = apnd - (real_T)maxdimlen;
                  }
                }
              }
            }
            st.site = &c_emlrtRSI;
            p = prod(&st, e_y);
            /*  Build integer multiplying the coefficient */
            k = d_y->size[0] * 4225;
            if (nx > k) {
              emlrtDynamicBoundsCheckR2012b(nx, 1, k, &emlrtBCI, (emlrtCTX)sp);
            }
            k = C->size[0] * 19 * 65 * 65;
            maxdimlen = (int32_T)(((uint32_T)o + b_i) + 1U);
            if ((maxdimlen < 1) || (maxdimlen > k)) {
              emlrtDynamicBoundsCheckR2012b(maxdimlen, 1, k, &b_emlrtBCI,
                                            (emlrtCTX)sp);
            }
            k = d_y->size[0] * 4225;
            if (nx > k) {
              emlrtDynamicBoundsCheckR2012b(nx, 1, k, &c_emlrtBCI,
                                            (emlrtCTX)sp);
            }
            d_y->data[nx - 1] =
                t * d_y->data[nx - 1] + C->data[maxdimlen - 1] * p;
            /*  y(l,n) = t * y(l,n) + C(o,i+1,n+1) * p; */
            if (*emlrtBreakCheckR2012bFlagVar != 0) {
              emlrtBreakCheckR2012b((emlrtCTX)sp);
            }
          }
        }
      }
      if (*emlrtBreakCheckR2012bFlagVar != 0) {
        emlrtBreakCheckR2012b((emlrtCTX)sp);
      }
    }
    /*  for l */
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b((emlrtCTX)sp);
    }
  }
  emxFree_real_T(&e_y);
  /*  for n */
  /*  Reshape output to be like input.  Also remove leading dimensions if L ==
   */
  /*  1, but leave row vectors as row vectors. */
  st.site = &d_emlrtRSI;
  b_st.site = &o_emlrtRSI;
  assertValidSizeArg(&b_st, szy);
  if (L * 65 * 65 != d_y->size[0] * 4225) {
    emlrtErrorWithMessageIdR2018a(
        &st, &d_emlrtRTEI, "Coder:MATLAB:getReshapeDims_notSameNumel",
        "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
  }
  f_y = *d_y;
  num[0] = L;
  num[1] = 65;
  num[2] = 65;
  num[3] = 1;
  f_y.size = &num[0];
  f_y.numDimensions = 4;
  st.site = &d_emlrtRSI;
  squeeze(&st, &f_y, y);
  emxFree_real_T(&d_y);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtCTX)sp);
}

/* End of code generation (ppc_val.c) */
