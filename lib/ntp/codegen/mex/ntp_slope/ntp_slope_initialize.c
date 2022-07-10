/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * ntp_slope_initialize.c
 *
 * Code generation for function 'ntp_slope_initialize'
 *
 */

/* Include files */
#include "ntp_slope_initialize.h"
#include "_coder_ntp_slope_mex.h"
#include "ntp_slope_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void ntp_slope_initialize(void)
{
  mex_InitInfAndNan();
  mexFunctionCreateRootTLS();
  emlrtClearAllocCountR2012b(emlrtRootTLSGlobal, false, 0U, NULL);
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (ntp_slope_initialize.c) */
