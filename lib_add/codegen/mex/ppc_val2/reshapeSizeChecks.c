/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * reshapeSizeChecks.c
 *
 * Code generation for function 'reshapeSizeChecks'
 *
 */

/* Include files */
#include "reshapeSizeChecks.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo g_emlrtRSI = {
    109,               /* lineNo */
    "computeDimsData", /* fcnName */
    "/Applications/MATLAB_R2021a.app/toolbox/eml/eml/+coder/+internal/"
    "reshapeSizeChecks.m" /* pathName */
};

static emlrtRTEInfo d_emlrtRTEI = {
    58,                   /* lineNo */
    23,                   /* colNo */
    "assertValidSizeArg", /* fName */
    "/Applications/MATLAB_R2021a.app/toolbox/eml/eml/+coder/+internal/"
    "assertValidSizeArg.m" /* pName */
};

static emlrtRTEInfo e_emlrtRTEI = {
    116,               /* lineNo */
    9,                 /* colNo */
    "computeDimsData", /* fName */
    "/Applications/MATLAB_R2021a.app/toolbox/eml/eml/+coder/+internal/"
    "reshapeSizeChecks.m" /* pName */
};

static emlrtRSInfo t_emlrtRSI = {
    18,            /* lineNo */
    "indexDivide", /* fcnName */
    "/Applications/MATLAB_R2021a.app/toolbox/eml/eml/+coder/+internal/"
    "indexDivide.m" /* pathName */
};

/* Function Declarations */
static int32_T div_s32(const emlrtStack *sp, int32_T numerator,
                       int32_T denominator);

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

int32_T computeDimsData(const emlrtStack *sp, int32_T nx, real_T varargin_1)
{
  emlrtStack st;
  int32_T calclen;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &g_emlrtRSI;
  if (varargin_1 != varargin_1) {
    emlrtErrorWithMessageIdR2018a(
        &st, &d_emlrtRTEI, "Coder:MATLAB:NonIntegerInput",
        "Coder:MATLAB:NonIntegerInput", 4, 12, MIN_int32_T, 12, MAX_int32_T);
  }
  if ((int32_T)varargin_1 > 0) {
    st.site = &t_emlrtRSI;
    calclen = div_s32(&st, nx, (int32_T)varargin_1);
    if (calclen > nx) {
      emlrtErrorWithMessageIdR2018a(sp, &e_emlrtRTEI,
                                    "Coder:builtins:AssertionFailed",
                                    "Coder:builtins:AssertionFailed", 0);
    }
  } else {
    calclen = 0;
  }
  return calclen;
}

/* End of code generation (reshapeSizeChecks.c) */
