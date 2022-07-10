/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * tbs_vertsolve_initialize.c
 *
 * Code generation for function 'tbs_vertsolve_initialize'
 *
 */

/* Include files */
#include "tbs_vertsolve_initialize.h"
#include "_coder_tbs_vertsolve_mex.h"
#include "rt_nonfinite.h"
#include "tbs_vertsolve_data.h"

/* Function Definitions */
void tbs_vertsolve_initialize(void)
{
  mex_InitInfAndNan();
  mexFunctionCreateRootTLS();
  emlrtClearAllocCountR2012b(emlrtRootTLSGlobal, false, 0U, NULL);
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (tbs_vertsolve_initialize.c) */
