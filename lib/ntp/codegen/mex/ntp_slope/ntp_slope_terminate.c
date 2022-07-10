/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * ntp_slope_terminate.c
 *
 * Code generation for function 'ntp_slope_terminate'
 *
 */

/* Include files */
#include "ntp_slope_terminate.h"
#include "_coder_ntp_slope_mex.h"
#include "ntp_slope_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void ntp_slope_atexit(void)
{
  mexFunctionCreateRootTLS();
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void ntp_slope_terminate(void)
{
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (ntp_slope_terminate.c) */
