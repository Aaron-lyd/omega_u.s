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
#include "assertValidSizeArg.h"
#include "colon.h"
#include "ppc_val2_data.h"
#include "ppc_val2_emxutil.h"
#include "ppc_val2_types.h"
#include "prod.h"
#include "reshapeSizeChecks.h"
#include "rt_nonfinite.h"
#include "squeeze.h"
#include "mwmathutil.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = {
    79,         /* lineNo */
    "ppc_val2", /* fcnName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/neutral-surfaces/lib/ppc/ppc_val2.m" /* pathName */
};

static emlrtRSInfo b_emlrtRSI = {
    114,        /* lineNo */
    "ppc_val2", /* fcnName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/neutral-surfaces/lib/ppc/ppc_val2.m" /* pathName */
};

static emlrtRSInfo c_emlrtRSI = {
    128,        /* lineNo */
    "ppc_val2", /* fcnName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/neutral-surfaces/lib/ppc/ppc_val2.m" /* pathName */
};

static emlrtRSInfo d_emlrtRSI = {
    142,        /* lineNo */
    "ppc_val2", /* fcnName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/neutral-surfaces/lib/ppc/ppc_val2.m" /* pathName */
};

static emlrtRSInfo e_emlrtRSI = {
    143,        /* lineNo */
    "ppc_val2", /* fcnName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/neutral-surfaces/lib/ppc/ppc_val2.m" /* pathName */
};

static emlrtRSInfo f_emlrtRSI = {
    29,                  /* lineNo */
    "reshapeSizeChecks", /* fcnName */
    "/Applications/MATLAB_R2021a.app/toolbox/eml/eml/+coder/+internal/"
    "reshapeSizeChecks.m" /* pathName */
};

static emlrtRSInfo m_emlrtRSI = {
    28,      /* lineNo */
    "colon", /* fcnName */
    "/Applications/MATLAB_R2021a.app/toolbox/eml/lib/matlab/ops/colon.m" /* pathName
                                                                          */
};

static emlrtRSInfo n_emlrtRSI = {
    117,     /* lineNo */
    "colon", /* fcnName */
    "/Applications/MATLAB_R2021a.app/toolbox/eml/lib/matlab/ops/colon.m" /* pathName
                                                                          */
};

static emlrtRSInfo p_emlrtRSI = {
    18,                  /* lineNo */
    "reshapeSizeChecks", /* fcnName */
    "/Applications/MATLAB_R2021a.app/toolbox/eml/eml/+coder/+internal/"
    "reshapeSizeChecks.m" /* pathName */
};

static emlrtMCInfo emlrtMCI = {
    75,         /* lineNo */
    1,          /* colNo */
    "ppc_val2", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/neutral-surfaces/lib/ppc/ppc_val2.m" /* pName */
};

static emlrtMCInfo b_emlrtMCI = {
    77,         /* lineNo */
    1,          /* colNo */
    "ppc_val2", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/neutral-surfaces/lib/ppc/ppc_val2.m" /* pName */
};

static emlrtRTEInfo c_emlrtRTEI = {
    127,        /* lineNo */
    25,         /* colNo */
    "ppc_val2", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/neutral-surfaces/lib/ppc/ppc_val2.m" /* pName */
};

static emlrtBCInfo emlrtBCI = {
    -1,         /* iFirst */
    -1,         /* iLast */
    91,         /* lineNo */
    15,         /* colNo */
    "x",        /* aName */
    "ppc_val2", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/neutral-surfaces/lib/ppc/ppc_val2.m", /* pName */
    0                                                     /* checkKind */
};

static emlrtBCInfo b_emlrtBCI = {
    -1,         /* iFirst */
    -1,         /* iLast */
    125,        /* lineNo */
    17,         /* colNo */
    "y",        /* aName */
    "ppc_val2", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/neutral-surfaces/lib/ppc/ppc_val2.m", /* pName */
    0                                                     /* checkKind */
};

static emlrtBCInfo c_emlrtBCI = {
    -1,         /* iFirst */
    -1,         /* iLast */
    126,        /* lineNo */
    17,         /* colNo */
    "z",        /* aName */
    "ppc_val2", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/neutral-surfaces/lib/ppc/ppc_val2.m", /* pName */
    0                                                     /* checkKind */
};

static emlrtBCInfo d_emlrtBCI = {
    -1,         /* iFirst */
    -1,         /* iLast */
    115,        /* lineNo */
    25,         /* colNo */
    "C",        /* aName */
    "ppc_val2", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/neutral-surfaces/lib/ppc/ppc_val2.m", /* pName */
    0                                                     /* checkKind */
};

static emlrtDCInfo emlrtDCI = {
    115,        /* lineNo */
    25,         /* colNo */
    "ppc_val2", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/neutral-surfaces/lib/ppc/ppc_val2.m", /* pName */
    1                                                     /* checkKind */
};

static emlrtBCInfo e_emlrtBCI = {
    -1,         /* iFirst */
    -1,         /* iLast */
    115,        /* lineNo */
    17,         /* colNo */
    "y",        /* aName */
    "ppc_val2", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/neutral-surfaces/lib/ppc/ppc_val2.m", /* pName */
    0                                                     /* checkKind */
};

static emlrtBCInfo f_emlrtBCI = {
    -1,         /* iFirst */
    -1,         /* iLast */
    116,        /* lineNo */
    25,         /* colNo */
    "D",        /* aName */
    "ppc_val2", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/neutral-surfaces/lib/ppc/ppc_val2.m", /* pName */
    0                                                     /* checkKind */
};

static emlrtDCInfo b_emlrtDCI = {
    116,        /* lineNo */
    25,         /* colNo */
    "ppc_val2", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/neutral-surfaces/lib/ppc/ppc_val2.m", /* pName */
    1                                                     /* checkKind */
};

static emlrtBCInfo g_emlrtBCI = {
    -1,         /* iFirst */
    -1,         /* iLast */
    116,        /* lineNo */
    17,         /* colNo */
    "z",        /* aName */
    "ppc_val2", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/neutral-surfaces/lib/ppc/ppc_val2.m", /* pName */
    0                                                     /* checkKind */
};

static emlrtBCInfo h_emlrtBCI = {
    -1,         /* iFirst */
    -1,         /* iLast */
    129,        /* lineNo */
    33,         /* colNo */
    "y",        /* aName */
    "ppc_val2", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/neutral-surfaces/lib/ppc/ppc_val2.m", /* pName */
    0                                                     /* checkKind */
};

static emlrtBCInfo i_emlrtBCI = {
    -1,         /* iFirst */
    -1,         /* iLast */
    129,        /* lineNo */
    41,         /* colNo */
    "C",        /* aName */
    "ppc_val2", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/neutral-surfaces/lib/ppc/ppc_val2.m", /* pName */
    0                                                     /* checkKind */
};

static emlrtBCInfo j_emlrtBCI = {
    -1,         /* iFirst */
    -1,         /* iLast */
    129,        /* lineNo */
    21,         /* colNo */
    "y",        /* aName */
    "ppc_val2", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/neutral-surfaces/lib/ppc/ppc_val2.m", /* pName */
    0                                                     /* checkKind */
};

static emlrtBCInfo k_emlrtBCI = {
    -1,         /* iFirst */
    -1,         /* iLast */
    130,        /* lineNo */
    33,         /* colNo */
    "z",        /* aName */
    "ppc_val2", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/neutral-surfaces/lib/ppc/ppc_val2.m", /* pName */
    0                                                     /* checkKind */
};

static emlrtBCInfo l_emlrtBCI = {
    -1,         /* iFirst */
    -1,         /* iLast */
    130,        /* lineNo */
    41,         /* colNo */
    "D",        /* aName */
    "ppc_val2", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/neutral-surfaces/lib/ppc/ppc_val2.m", /* pName */
    0                                                     /* checkKind */
};

static emlrtBCInfo m_emlrtBCI = {
    -1,         /* iFirst */
    -1,         /* iLast */
    130,        /* lineNo */
    21,         /* colNo */
    "z",        /* aName */
    "ppc_val2", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/neutral-surfaces/lib/ppc/ppc_val2.m", /* pName */
    0                                                     /* checkKind */
};

static emlrtRTEInfo k_emlrtRTEI = {
    80,         /* lineNo */
    1,          /* colNo */
    "ppc_val2", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/neutral-surfaces/lib/ppc/ppc_val2.m" /* pName */
};

static emlrtRTEInfo l_emlrtRTEI = {
    81,         /* lineNo */
    1,          /* colNo */
    "ppc_val2", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/neutral-surfaces/lib/ppc/ppc_val2.m" /* pName */
};

static emlrtRTEInfo m_emlrtRTEI = {
    28,      /* lineNo */
    9,       /* colNo */
    "colon", /* fName */
    "/Applications/MATLAB_R2021a.app/toolbox/eml/lib/matlab/ops/colon.m" /* pName
                                                                          */
};

static emlrtRTEInfo n_emlrtRTEI = {
    1,          /* lineNo */
    11,         /* colNo */
    "ppc_val2", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/neutral-surfaces/lib/ppc/ppc_val2.m" /* pName */
};

static emlrtRTEInfo o_emlrtRTEI = {
    1,          /* lineNo */
    13,         /* colNo */
    "ppc_val2", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/neutral-surfaces/lib/ppc/ppc_val2.m" /* pName */
};

static emlrtRTEInfo p_emlrtRTEI = {
    128,        /* lineNo */
    30,         /* colNo */
    "ppc_val2", /* fName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/neutral-surfaces/lib/ppc/ppc_val2.m" /* pName */
};

static emlrtRSInfo r_emlrtRSI = {
    75,         /* lineNo */
    "ppc_val2", /* fcnName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/neutral-surfaces/lib/ppc/ppc_val2.m" /* pathName */
};

static emlrtRSInfo s_emlrtRSI = {
    77,         /* lineNo */
    "ppc_val2", /* fcnName */
    "/Users/yandonglang/Library/CloudStorage/OneDrive-UNSW/1. YL_PHD/10. "
    "omega_hel code/neutral-surfaces/lib/ppc/ppc_val2.m" /* pathName */
};

/* Function Declarations */
static void b_error(const emlrtStack *sp, const mxArray *b,
                    emlrtMCInfo *location);

static void error(const emlrtStack *sp, const mxArray *b, const mxArray *c,
                  emlrtMCInfo *location);

/* Function Definitions */
static void b_error(const emlrtStack *sp, const mxArray *b,
                    emlrtMCInfo *location)
{
  const mxArray *pArray;
  pArray = b;
  emlrtCallMATLABR2012b((emlrtCTX)sp, 0, NULL, 1, &pArray,
                        (const char_T *)"error", true, location);
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

void ppc_val2(const emlrtStack *sp, const real_T X[20],
              const emxArray_real_T *C, const emxArray_real_T *D,
              const real_T x_data[], const int32_T x_size[2], real_T d,
              emxArray_real_T *y, emxArray_real_T *z)
{
  static const int32_T iv[2] = {1, 34};
  static const int32_T iv1[2] = {1, 31};
  static const char_T u[34] = {'x', ' ', 'm', 'u', 's', 't', ' ', 'h', 'a',
                               'v', 'e', ' ', '%', 'd', ' ', 'c', 'o', 'l',
                               'u', 'm', 'n', 's', ' ', 'o', 'r', ' ', '1',
                               ' ', 'c', 'o', 'l', 'u', 'm', 'n'};
  static const char_T b_u[31] = {'C', ' ', 'a', 'n', 'd', ' ', 'D', ' ',
                                 'm', 'u', 's', 't', ' ', 'h', 'a', 'v',
                                 'e', ' ', 't', 'h', 'e', ' ', 's', 'a',
                                 'm', 'e', ' ', 's', 'i', 'z', 'e'};
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack st;
  emxArray_real_T h_y;
  emxArray_real_T *b_z;
  emxArray_real_T *f_y;
  emxArray_real_T *g_y;
  const mxArray *b_y;
  const mxArray *d_y;
  const mxArray *e_y;
  const mxArray *m;
  real_T szy[4];
  real_T a;
  real_T p;
  real_T t;
  real_T xln;
  int32_T b_num_tmp[4];
  int32_T num_tmp[4];
  int32_T L;
  int32_T O;
  int32_T b_i;
  int32_T i;
  int32_T i1;
  int32_T i2;
  int32_T l;
  int32_T maxdimlen;
  int32_T n;
  int32_T nx;
  int32_T o;
  int32_T sz_idx_1;
  int32_T xM;
  boolean_T x[4];
  boolean_T c_y;
  boolean_T exitg1;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtCTX)sp);
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
    d_y = NULL;
    m = emlrtCreateDoubleScalar(4225.0);
    emlrtAssign(&d_y, m);
    st.site = &r_emlrtRSI;
    error(&st, b_y, d_y, &emlrtMCI);
  }
  x[0] = ((int8_T)C->size[0] == D->size[0]);
  x[1] = true;
  x[2] = true;
  x[3] = true;
  c_y = true;
  sz_idx_1 = 0;
  exitg1 = false;
  while ((!exitg1) && (sz_idx_1 < 4)) {
    if (!x[sz_idx_1]) {
      c_y = false;
      exitg1 = true;
    } else {
      sz_idx_1++;
    }
  }
  if (!c_y) {
    e_y = NULL;
    m = emlrtCreateCharArray(2, &iv1[0]);
    emlrtInitCharArrayR2013a((emlrtCTX)sp, 31, m, &b_u[0]);
    emlrtAssign(&e_y, m);
    st.site = &s_emlrtRSI;
    b_error(&st, e_y, &b_emlrtMCI);
  }
  st.site = &emlrtRSI;
  nx = x_size[0] * x_size[1];
  b_st.site = &f_emlrtRSI;
  sz_idx_1 = computeDimsData(&b_st, nx, L);
  n = x_size[0];
  if (x_size[1] > x_size[0]) {
    n = x_size[1];
  }
  maxdimlen = muIntScalarMax_sint32(nx, n);
  if (L > maxdimlen) {
    emlrtErrorWithMessageIdR2018a(&st, &b_emlrtRTEI,
                                  "Coder:toolbox:reshape_emptyReshapeLimit",
                                  "Coder:toolbox:reshape_emptyReshapeLimit", 0);
  }
  if (sz_idx_1 > maxdimlen) {
    emlrtErrorWithMessageIdR2018a(&st, &b_emlrtRTEI,
                                  "Coder:toolbox:reshape_emptyReshapeLimit",
                                  "Coder:toolbox:reshape_emptyReshapeLimit", 0);
  }
  i = L * sz_idx_1;
  if (i != nx) {
    emlrtErrorWithMessageIdR2018a(
        &st, &emlrtRTEI, "Coder:MATLAB:getReshapeDims_notSameNumel",
        "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
  }
  emxInit_real_T(&st, &f_y, 2, &n_emlrtRTEI, true);
  i1 = f_y->size[0] * f_y->size[1];
  f_y->size[0] = L;
  f_y->size[1] = 4225;
  emxEnsureCapacity_real_T(sp, f_y, i1, &k_emlrtRTEI);
  maxdimlen = L * 4225;
  for (i1 = 0; i1 < maxdimlen; i1++) {
    f_y->data[i1] = rtNaN;
  }
  emxInit_real_T(sp, &b_z, 2, &o_emlrtRTEI, true);
  i1 = b_z->size[0] * b_z->size[1];
  b_z->size[0] = L;
  b_z->size[1] = 4225;
  emxEnsureCapacity_real_T(sp, b_z, i1, &l_emlrtRTEI);
  for (i1 = 0; i1 < maxdimlen; i1++) {
    b_z->data[i1] = rtNaN;
  }
  /*  Evaluate each piecewise polynomial */
  nx = (xM > 1);
  /*  used for linear indexing */
  /*  used for linear indexing */
  xM = (int8_T)C->size[0] * 19;
  /*  used for linear indexing */
  /*  used for linear indexing */
  emxInit_real_T(sp, &g_y, 2, &p_emlrtRTEI, true);
  for (n = 0; n < 4225; n++) {
    for (l = 0; l < L; l++) {
      i1 = n * L;
      i2 = (l + i1 * nx) + 1;
      if (i2 > i) {
        emlrtDynamicBoundsCheckR2012b(i2, 1, i, &emlrtBCI, (emlrtCTX)sp);
      }
      xln = x_data[i2 - 1];
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
        sz_idx_1 = 20;
        /*  Result will be <= M always */
        while (b_i < sz_idx_1) {
          maxdimlen =
              (int32_T)muDoubleScalarFloor((real_T)(b_i + sz_idx_1) / 2.0);
          if (X[maxdimlen - 1] < xln) {
            /*  [X(j,n) or X(j,1) as appropriate]  < xln */
            b_i = maxdimlen + 1;
          } else {
            sz_idx_1 = maxdimlen;
          }
          if (*emlrtBreakCheckR2012bFlagVar != 0) {
            emlrtBreakCheckR2012b((emlrtCTX)sp);
          }
        }
        sz_idx_1 = (l + i1) + 1;
        if (b_i == 1) {
          /*  Note: X(nX + 1) == xln   is guaranteed */
          if (muDoubleScalarIsNaN(d)) {
            i1 = g_y->size[0] * g_y->size[1];
            g_y->size[0] = 1;
            g_y->size[1] = 1;
            emxEnsureCapacity_real_T(sp, g_y, i1, &m_emlrtRTEI);
            g_y->data[0] = rtNaN;
          } else if (d < 1.0) {
            g_y->size[0] = 1;
            g_y->size[1] = 0;
          } else if (muDoubleScalarIsInf(d) && (1.0 == d)) {
            i1 = g_y->size[0] * g_y->size[1];
            g_y->size[0] = 1;
            g_y->size[1] = 1;
            emxEnsureCapacity_real_T(sp, g_y, i1, &m_emlrtRTEI);
            g_y->data[0] = rtNaN;
          } else {
            i1 = g_y->size[0] * g_y->size[1];
            g_y->size[0] = 1;
            maxdimlen = (int32_T)muDoubleScalarFloor(d - 1.0);
            g_y->size[1] = maxdimlen + 1;
            emxEnsureCapacity_real_T(sp, g_y, i1, &m_emlrtRTEI);
            for (i1 = 0; i1 <= maxdimlen; i1++) {
              g_y->data[i1] = (real_T)i1 + 1.0;
            }
          }
          st.site = &b_emlrtRSI;
          p = prod(&st, g_y);
          /*  Build integer multiplying the coefficient */
          i1 = C->size[0] * 19 * 65 * 65;
          a = ((real_T)O - d) + (real_T)(n * xM);
          xln = (int32_T)muDoubleScalarFloor(a);
          if (a != xln) {
            emlrtIntegerCheckR2012b(a, &emlrtDCI, (emlrtCTX)sp);
          }
          if (((int32_T)a < 1) || ((int32_T)a > i1)) {
            emlrtDynamicBoundsCheckR2012b((int32_T)a, 1, i1, &d_emlrtBCI,
                                          (emlrtCTX)sp);
          }
          i1 = f_y->size[0] * 4225;
          if (sz_idx_1 > i1) {
            emlrtDynamicBoundsCheckR2012b(sz_idx_1, 1, i1, &e_emlrtBCI,
                                          (emlrtCTX)sp);
          }
          f_y->data[sz_idx_1 - 1] = C->data[(int32_T)a - 1] * p;
          /*  y(l,n) = C(O - d, i, n) * p; */
          i1 = D->size[0] * 19 * 65 * 65;
          if (a != xln) {
            emlrtIntegerCheckR2012b(a, &b_emlrtDCI, (emlrtCTX)sp);
          }
          if (((int32_T)a < 1) || ((int32_T)a > i1)) {
            emlrtDynamicBoundsCheckR2012b((int32_T)a, 1, i1, &f_emlrtBCI,
                                          (emlrtCTX)sp);
          }
          i1 = b_z->size[0] * 4225;
          if (sz_idx_1 > i1) {
            emlrtDynamicBoundsCheckR2012b(sz_idx_1, 1, i1, &g_emlrtBCI,
                                          (emlrtCTX)sp);
          }
          b_z->data[sz_idx_1 - 1] = D->data[(int32_T)a - 1] * p;
        } else {
          /*  Evaluate this piece of the polynomial (see ppval.m) */
          t = xln - X[b_i - 2];
          /*  Switch to local coordinates */
          /*  Overload variable i, to speed with indexing */
          /*  subtract 1 from i so that 1 <= i <= M-1, and X(i) <= xln < X(i+1),
           */
          /*  subtract another 1 for indexing.  */
          b_i = (b_i - 2) * O + n * xM;
          i1 = f_y->size[0] * 4225;
          if (sz_idx_1 > i1) {
            emlrtDynamicBoundsCheckR2012b(sz_idx_1, 1, i1, &b_emlrtBCI,
                                          (emlrtCTX)sp);
          }
          f_y->data[sz_idx_1 - 1] = 0.0;
          i1 = b_z->size[0] * 4225;
          if (sz_idx_1 > i1) {
            emlrtDynamicBoundsCheckR2012b(sz_idx_1, 1, i1, &c_emlrtBCI,
                                          (emlrtCTX)sp);
          }
          b_z->data[sz_idx_1 - 1] = 0.0;
          a = (real_T)O - d;
          i1 = (int32_T)a;
          emlrtForLoopVectorCheckR2021a(1.0, 1.0, a, mxDOUBLE_CLASS, (int32_T)a,
                                        &c_emlrtRTEI, (emlrtCTX)sp);
          for (o = 0; o < i1; o++) {
            st.site = &c_emlrtRSI;
            xln = (real_T)O - ((real_T)o + 1.0);
            a = (xln - d) + 1.0;
            b_st.site = &m_emlrtRSI;
            if (muDoubleScalarIsNaN(a)) {
              i2 = g_y->size[0] * g_y->size[1];
              g_y->size[0] = 1;
              g_y->size[1] = 1;
              emxEnsureCapacity_real_T(&b_st, g_y, i2, &m_emlrtRTEI);
              g_y->data[0] = rtNaN;
            } else if ((int32_T)xln < a) {
              g_y->size[0] = 1;
              g_y->size[1] = 0;
            } else if (muDoubleScalarIsInf(a) && (a == (int32_T)xln)) {
              i2 = g_y->size[0] * g_y->size[1];
              g_y->size[0] = 1;
              g_y->size[1] = 1;
              emxEnsureCapacity_real_T(&b_st, g_y, i2, &m_emlrtRTEI);
              g_y->data[0] = rtNaN;
            } else if (muDoubleScalarFloor(a) == a) {
              i2 = g_y->size[0] * g_y->size[1];
              g_y->size[0] = 1;
              maxdimlen = (int32_T)((real_T)(int32_T)xln - a);
              g_y->size[1] = maxdimlen + 1;
              emxEnsureCapacity_real_T(&b_st, g_y, i2, &m_emlrtRTEI);
              for (i2 = 0; i2 <= maxdimlen; i2++) {
                g_y->data[i2] = a + (real_T)i2;
              }
            } else {
              c_st.site = &n_emlrtRSI;
              eml_float_colon(&c_st, a, (int32_T)xln, g_y);
            }
            st.site = &c_emlrtRSI;
            p = prod(&st, g_y);
            /*  Build integer multiplying the coefficient */
            i2 = f_y->size[0] * 4225;
            if (sz_idx_1 > i2) {
              emlrtDynamicBoundsCheckR2012b(sz_idx_1, 1, i2, &h_emlrtBCI,
                                            (emlrtCTX)sp);
            }
            i2 = C->size[0] * 19 * 65 * 65;
            maxdimlen = (int32_T)(((uint32_T)o + b_i) + 1U);
            if ((maxdimlen < 1) || (maxdimlen > i2)) {
              emlrtDynamicBoundsCheckR2012b(maxdimlen, 1, i2, &i_emlrtBCI,
                                            (emlrtCTX)sp);
            }
            i2 = f_y->size[0] * 4225;
            if (sz_idx_1 > i2) {
              emlrtDynamicBoundsCheckR2012b(sz_idx_1, 1, i2, &j_emlrtBCI,
                                            (emlrtCTX)sp);
            }
            f_y->data[sz_idx_1 - 1] =
                t * f_y->data[sz_idx_1 - 1] + C->data[maxdimlen - 1] * p;
            /*  y(l,n) = t * y(l,n) + C(o,i+1,n+1) * p; */
            i2 = b_z->size[0] * 4225;
            if (sz_idx_1 > i2) {
              emlrtDynamicBoundsCheckR2012b(sz_idx_1, 1, i2, &k_emlrtBCI,
                                            (emlrtCTX)sp);
            }
            i2 = D->size[0] * 19 * 65 * 65;
            if (maxdimlen > i2) {
              emlrtDynamicBoundsCheckR2012b(maxdimlen, 1, i2, &l_emlrtBCI,
                                            (emlrtCTX)sp);
            }
            i2 = b_z->size[0] * 4225;
            if (sz_idx_1 > i2) {
              emlrtDynamicBoundsCheckR2012b(sz_idx_1, 1, i2, &m_emlrtBCI,
                                            (emlrtCTX)sp);
            }
            b_z->data[sz_idx_1 - 1] =
                t * b_z->data[sz_idx_1 - 1] + D->data[maxdimlen - 1] * p;
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
  emxFree_real_T(&g_y);
  /*  for n */
  /*  Reshape output to be like input.  Also remove leading dimensions if L ==
   */
  /*  1, but leave row vectors as row vectors. */
  st.site = &d_emlrtRSI;
  b_st.site = &p_emlrtRSI;
  assertValidSizeArg(&b_st, szy);
  i = L * 65 * 65;
  if (i != f_y->size[0] * 4225) {
    emlrtErrorWithMessageIdR2018a(
        &st, &emlrtRTEI, "Coder:MATLAB:getReshapeDims_notSameNumel",
        "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
  }
  h_y = *f_y;
  num_tmp[0] = L;
  num_tmp[1] = 65;
  num_tmp[2] = 65;
  num_tmp[3] = 1;
  h_y.size = &num_tmp[0];
  h_y.numDimensions = 4;
  st.site = &d_emlrtRSI;
  squeeze(&st, &h_y, y);
  st.site = &e_emlrtRSI;
  b_st.site = &p_emlrtRSI;
  assertValidSizeArg(&b_st, szy);
  emxFree_real_T(&f_y);
  if (i != b_z->size[0] * 4225) {
    emlrtErrorWithMessageIdR2018a(
        &st, &emlrtRTEI, "Coder:MATLAB:getReshapeDims_notSameNumel",
        "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
  }
  h_y = *b_z;
  b_num_tmp[0] = L;
  b_num_tmp[1] = 65;
  b_num_tmp[2] = 65;
  b_num_tmp[3] = 1;
  h_y.size = &b_num_tmp[0];
  h_y.numDimensions = 4;
  st.site = &e_emlrtRSI;
  squeeze(&st, &h_y, z);
  emxFree_real_T(&b_z);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtCTX)sp);
}

/* End of code generation (ppc_val2.c) */
