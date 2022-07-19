/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_ppc_val_mex.c
 *
 * Code generation for function '_coder_ppc_val_mex'
 *
 */

/* Include files */
#include "_coder_ppc_val_mex.h"
#include "_coder_ppc_val_api.h"
#include "ppc_val_data.h"
#include "ppc_val_initialize.h"
#include "ppc_val_terminate.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  mexAtExit(&ppc_val_atexit);
  /* Module initialization. */
  ppc_val_initialize();
  /* Dispatch the entry-point. */
  ppc_val_mexFunction(nlhs, plhs, nrhs, prhs);
  /* Module termination. */
  ppc_val_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2021a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL);
  return emlrtRootTLSGlobal;
}

void ppc_val_mexFunction(int32_T nlhs, mxArray *plhs[1], int32_T nrhs,
                         const mxArray *prhs[4])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  const mxArray *outputs;
  st.tls = emlrtRootTLSGlobal;
  /* Check for proper number of arguments. */
  if (nrhs != 4) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 4, 4,
                        7, "ppc_val");
  }
  if (nlhs > 1) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 7,
                        "ppc_val");
  }
  /* Call the function. */
  ppc_val_api(prhs, &outputs);
  /* Copy over outputs to the caller. */
  emlrtReturnArrays(1, &plhs[0], &outputs);
}

/* End of code generation (_coder_ppc_val_mex.c) */
