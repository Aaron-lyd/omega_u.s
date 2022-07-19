/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * assertValidSizeArg.c
 *
 * Code generation for function 'assertValidSizeArg'
 *
 */

/* Include files */
#include "assertValidSizeArg.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"

/* Variable Definitions */
static emlrtRTEInfo i_emlrtRTEI = {
    64,                   /* lineNo */
    15,                   /* colNo */
    "assertValidSizeArg", /* fName */
    "/Applications/MATLAB_R2021a.app/toolbox/eml/eml/+coder/+internal/"
    "assertValidSizeArg.m" /* pName */
};

static emlrtRTEInfo j_emlrtRTEI = {
    49,                   /* lineNo */
    19,                   /* colNo */
    "assertValidSizeArg", /* fName */
    "/Applications/MATLAB_R2021a.app/toolbox/eml/eml/+coder/+internal/"
    "assertValidSizeArg.m" /* pName */
};

/* Function Definitions */
void assertValidSizeArg(const emlrtStack *sp, const real_T varargin_1[4])
{
  real_T n;
  int32_T k;
  boolean_T exitg1;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 4)) {
    if ((varargin_1[k] != varargin_1[k]) ||
        muDoubleScalarIsInf(varargin_1[k])) {
      emlrtErrorWithMessageIdR2018a(
          sp, &j_emlrtRTEI,
          "Coder:toolbox:eml_assert_valid_size_arg_invalidSizeVector",
          "Coder:toolbox:eml_assert_valid_size_arg_invalidSizeVector", 4, 12,
          MIN_int32_T, 12, MAX_int32_T);
    } else {
      k++;
    }
  }
  if (varargin_1[0] <= 0.0) {
    n = 0.0;
  } else {
    n = varargin_1[0];
  }
  if (varargin_1[1] <= 0.0) {
    n = 0.0;
  } else {
    n *= varargin_1[1];
  }
  if (varargin_1[2] <= 0.0) {
    n = 0.0;
  } else {
    n *= varargin_1[2];
  }
  if (varargin_1[3] <= 0.0) {
    n = 0.0;
  } else {
    n *= varargin_1[3];
  }
  if (!(n <= 2.147483647E+9)) {
    emlrtErrorWithMessageIdR2018a(sp, &i_emlrtRTEI, "Coder:MATLAB:pmaxsize",
                                  "Coder:MATLAB:pmaxsize", 0);
  }
}

/* End of code generation (assertValidSizeArg.c) */
