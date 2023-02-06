/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * obs_vertsolve_terminate.c
 *
 * Code generation for function 'obs_vertsolve_terminate'
 *
 */

/* Include files */
#include "obs_vertsolve_terminate.h"
#include "_coder_obs_vertsolve_mex.h"
#include "obs_vertsolve_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void obs_vertsolve_atexit(void)
{
  mexFunctionCreateRootTLS();
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void obs_vertsolve_terminate(void)
{
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (obs_vertsolve_terminate.c) */
