/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * omega_vertsolve_terminate.c
 *
 * Code generation for function 'omega_vertsolve_terminate'
 *
 */

/* Include files */
#include "omega_vertsolve_terminate.h"
#include "_coder_omega_vertsolve_mex.h"
#include "omega_vertsolve_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void omega_vertsolve_atexit(void)
{
  mexFunctionCreateRootTLS();
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void omega_vertsolve_terminate(void)
{
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (omega_vertsolve_terminate.c) */
