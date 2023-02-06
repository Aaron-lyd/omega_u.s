/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * obs_vertsolve_initialize.c
 *
 * Code generation for function 'obs_vertsolve_initialize'
 *
 */

/* Include files */
#include "obs_vertsolve_initialize.h"
#include "_coder_obs_vertsolve_mex.h"
#include "obs_vertsolve_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void obs_vertsolve_initialize(void)
{
  mex_InitInfAndNan();
  mexFunctionCreateRootTLS();
  emlrtClearAllocCountR2012b(emlrtRootTLSGlobal, false, 0U, NULL);
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (obs_vertsolve_initialize.c) */
