/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_ntp_slope_mex.c
 *
 * Code generation for function '_coder_ntp_slope_mex'
 *
 */

/* Include files */
#include "_coder_ntp_slope_mex.h"
#include "_coder_ntp_slope_api.h"
#include "ntp_slope_data.h"
#include "ntp_slope_initialize.h"
#include "ntp_slope_terminate.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  mexAtExit(&ntp_slope_atexit);
  /* Module initialization. */
  ntp_slope_initialize();
  /* Dispatch the entry-point. */
  unsafe_ntp_slope_mexFunction(nlhs, plhs, nrhs, prhs);
  /* Module termination. */
  ntp_slope_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2021a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL);
  return emlrtRootTLSGlobal;
}

void unsafe_ntp_slope_mexFunction(int32_T nlhs, mxArray *plhs[2], int32_T nrhs,
                                  const mxArray *prhs[7])
{
  const mxArray *outputs[2];
  int32_T b_nlhs;
  /* Check for proper number of arguments. */
  if (nrhs != 7) {
    emlrtErrMsgIdAndTxt(emlrtRootTLSGlobal, "EMLRT:runTime:WrongNumberOfInputs",
                        5, 12, 7, 4, 9, "ntp_slope");
  }
  if (nlhs > 2) {
    emlrtErrMsgIdAndTxt(emlrtRootTLSGlobal,
                        "EMLRT:runTime:TooManyOutputArguments", 3, 4, 9,
                        "ntp_slope");
  }
  /* Call the function. */
  ntp_slope_api(prhs, nlhs, outputs);
  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }
  emlrtReturnArrays(b_nlhs, &plhs[0], &outputs[0]);
}

/* End of code generation (_coder_ntp_slope_mex.c) */
