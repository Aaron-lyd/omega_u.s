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

/* Function Declarations */
static void omega_vertsolve_once(void);

/* Function Definitions */
static void omega_vertsolve_once(void)
{
  mex_InitInfAndNan();
}

void omega_vertsolve_initialize(void)
{
  mexFunctionCreateRootTLS();
  emlrtClearAllocCountR2012b(emlrtRootTLSGlobal, false, 0U, NULL);
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  if (emlrtFirstTimeR2012b(emlrtRootTLSGlobal)) {
    omega_vertsolve_once();
  }
}

/* End of code generation (omega_vertsolve_initialize.c) */
