/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_omega_vertsolve_mex.c
 *
 * Code generation for function '_coder_omega_vertsolve_mex'
 *
 */

/* Include files */
#include "_coder_omega_vertsolve_mex.h"
#include "_coder_omega_vertsolve_api.h"
#include "omega_vertsolve_data.h"
#include "omega_vertsolve_initialize.h"
#include "omega_vertsolve_terminate.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  mexAtExit(&omega_vertsolve_atexit);
  /* Module initialization. */
  omega_vertsolve_initialize();
  /* Dispatch the entry-point. */
  unsafe_omega_vertsolve_mexFunction(nlhs, plhs, nrhs, prhs);
  /* Module termination. */
  omega_vertsolve_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2022a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL, "windows-1252", true);
  return emlrtRootTLSGlobal;
}

void unsafe_omega_vertsolve_mexFunction(int32_T nlhs, mxArray *plhs[3],
                                        int32_T nrhs, const mxArray *prhs[9])
{
  const mxArray *outputs[3];
  int32_T i;
  /* Check for proper number of arguments. */
  if (nrhs != 9) {
    emlrtErrMsgIdAndTxt(emlrtRootTLSGlobal, "EMLRT:runTime:WrongNumberOfInputs",
                        5, 12, 9, 4, 15, "omega_vertsolve");
  }
  if (nlhs > 3) {
    emlrtErrMsgIdAndTxt(emlrtRootTLSGlobal,
                        "EMLRT:runTime:TooManyOutputArguments", 3, 4, 15,
                        "omega_vertsolve");
  }
  /* Call the function. */
  omega_vertsolve_api(prhs, nlhs, outputs);
  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    i = 1;
  } else {
    i = nlhs;
  }
  emlrtReturnArrays(i, &plhs[0], &outputs[0]);
}

/* End of code generation (_coder_omega_vertsolve_mex.c) */
