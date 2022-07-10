/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * tbs_vertsolve_terminate.c
 *
 * Code generation for function 'tbs_vertsolve_terminate'
 *
 */

/* Include files */
#include "tbs_vertsolve_terminate.h"
#include "_coder_tbs_vertsolve_mex.h"
#include "rt_nonfinite.h"
#include "tbs_vertsolve_data.h"

/* Function Definitions */
void tbs_vertsolve_atexit(void)
{
  mexFunctionCreateRootTLS();
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void tbs_vertsolve_terminate(void)
{
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (tbs_vertsolve_terminate.c) */
