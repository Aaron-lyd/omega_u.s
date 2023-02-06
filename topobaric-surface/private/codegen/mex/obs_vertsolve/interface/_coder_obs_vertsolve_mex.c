/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_obs_vertsolve_mex.c
 *
 * Code generation for function '_coder_obs_vertsolve_mex'
 *
 */

/* Include files */
#include "_coder_obs_vertsolve_mex.h"
#include "_coder_obs_vertsolve_api.h"
#include "obs_vertsolve_data.h"
#include "obs_vertsolve_initialize.h"
#include "obs_vertsolve_terminate.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  mexAtExit(&obs_vertsolve_atexit);
  /* Module initialization. */
  obs_vertsolve_initialize();
  /* Dispatch the entry-point. */
  unsafe_obs_vertsolve_mexFunction(nlhs, plhs, nrhs, prhs);
  /* Module termination. */
  obs_vertsolve_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2021a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL);
  return emlrtRootTLSGlobal;
}

void unsafe_obs_vertsolve_mexFunction(int32_T nlhs, mxArray *plhs[3],
                                      int32_T nrhs, const mxArray *prhs[12])
{
  const mxArray *outputs[3];
  int32_T b_nlhs;
  /* Check for proper number of arguments. */
  if (nrhs != 12) {
    emlrtErrMsgIdAndTxt(emlrtRootTLSGlobal, "EMLRT:runTime:WrongNumberOfInputs",
                        5, 12, 12, 4, 13, "obs_vertsolve");
  }
  if (nlhs > 3) {
    emlrtErrMsgIdAndTxt(emlrtRootTLSGlobal,
                        "EMLRT:runTime:TooManyOutputArguments", 3, 4, 13,
                        "obs_vertsolve");
  }
  /* Call the function. */
  obs_vertsolve_api(prhs, nlhs, outputs);
  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }
  emlrtReturnArrays(b_nlhs, &plhs[0], &outputs[0]);
}

/* End of code generation (_coder_obs_vertsolve_mex.c) */
