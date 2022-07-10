/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * omega_vertsolve_initialize.c
 *
 * Code generation for function 'omega_vertsolve_initialize'
 *
 */

/* Include files */
#include "omega_vertsolve_initialize.h"
#include "_coder_omega_vertsolve_mex.h"
#include "omega_vertsolve_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void omega_vertsolve_initialize(void)
{
  mex_InitInfAndNan();
  mexFunctionCreateRootTLS();
  emlrtClearAllocCountR2012b(emlrtRootTLSGlobal, false, 0U, NULL);
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (omega_vertsolve_initialize.c) */
