/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_bfs_conncomp1_wet_mex.c
 *
 * Code generation for function '_coder_bfs_conncomp1_wet_mex'
 *
 */

/* Include files */
#include "_coder_bfs_conncomp1_wet_mex.h"
#include "_coder_bfs_conncomp1_wet_api.h"
#include "bfs_conncomp1_wet_data.h"
#include "bfs_conncomp1_wet_initialize.h"
#include "bfs_conncomp1_wet_terminate.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  mexAtExit(&bfs_conncomp1_wet_atexit);
  /* Module initialization. */
  bfs_conncomp1_wet_initialize();
  /* Dispatch the entry-point. */
  unsafe_bfs_conncomp1_wet_mexFunction(nlhs, plhs, nrhs, prhs);
  /* Module termination. */
  bfs_conncomp1_wet_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2021a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL);
  return emlrtRootTLSGlobal;
}

void unsafe_bfs_conncomp1_wet_mexFunction(int32_T nlhs, mxArray *plhs[6],
                                          int32_T nrhs, const mxArray *prhs[12])
{
  const mxArray *outputs[6];
  int32_T b_nlhs;
  /* Check for proper number of arguments. */
  if (nrhs != 12) {
    emlrtErrMsgIdAndTxt(emlrtRootTLSGlobal, "EMLRT:runTime:WrongNumberOfInputs",
                        5, 12, 12, 4, 17, "bfs_conncomp1_wet");
  }
  if (nlhs > 6) {
    emlrtErrMsgIdAndTxt(emlrtRootTLSGlobal,
                        "EMLRT:runTime:TooManyOutputArguments", 3, 4, 17,
                        "bfs_conncomp1_wet");
  }
  /* Call the function. */
  bfs_conncomp1_wet_api(prhs, nlhs, outputs);
  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }
  emlrtReturnArrays(b_nlhs, &plhs[0], &outputs[0]);
}

/* End of code generation (_coder_bfs_conncomp1_wet_mex.c) */
