/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * prod.c
 *
 * Code generation for function 'prod'
 *
 */

/* Include files */
#include "prod.h"
#include "eml_int_forloop_overflow_check.h"
#include "ppc_val_types.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo
    g_emlrtRSI =
        {
            11,     /* lineNo */
            "prod", /* fcnName */
            "/Applications/MATLAB_R2021a.app/toolbox/eml/lib/matlab/datafun/"
            "prod.m" /* pathName */
};

static emlrtRSInfo h_emlrtRSI = {
    99,        /* lineNo */
    "sumprod", /* fcnName */
    "/Applications/MATLAB_R2021a.app/toolbox/eml/lib/matlab/datafun/private/"
    "sumprod.m" /* pathName */
};

static emlrtRSInfo i_emlrtRSI = {
    138,                     /* lineNo */
    "combineVectorElements", /* fcnName */
    "/Applications/MATLAB_R2021a.app/toolbox/eml/lib/matlab/datafun/private/"
    "combineVectorElements.m" /* pathName */
};

static emlrtRSInfo j_emlrtRSI = {
    198,                /* lineNo */
    "colMajorFlatIter", /* fcnName */
    "/Applications/MATLAB_R2021a.app/toolbox/eml/lib/matlab/datafun/private/"
    "combineVectorElements.m" /* pathName */
};

static emlrtRSInfo k_emlrtRSI = {
    21,                               /* lineNo */
    "eml_int_forloop_overflow_check", /* fcnName */
    "/Applications/MATLAB_R2021a.app/toolbox/eml/lib/matlab/eml/"
    "eml_int_forloop_overflow_check.m" /* pathName */
};

/* Function Definitions */
real_T prod(const emlrtStack *sp, const emxArray_real_T *x)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack st;
  real_T y;
  int32_T k;
  int32_T vlen;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &g_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  b_st.site = &h_emlrtRSI;
  vlen = x->size[1];
  if (x->size[1] == 0) {
    y = 1.0;
  } else {
    c_st.site = &i_emlrtRSI;
    y = x->data[0];
    d_st.site = &j_emlrtRSI;
    if ((2 <= x->size[1]) && (x->size[1] > 2147483646)) {
      e_st.site = &k_emlrtRSI;
      check_forloop_overflow_error(&e_st);
    }
    for (k = 2; k <= vlen; k++) {
      y *= x->data[k - 1];
    }
  }
  return y;
}

/* End of code generation (prod.c) */
