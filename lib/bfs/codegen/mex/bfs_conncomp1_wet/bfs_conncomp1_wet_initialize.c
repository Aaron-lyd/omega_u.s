/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * bfs_conncomp1_wet_initialize.c
 *
 * Code generation for function 'bfs_conncomp1_wet_initialize'
 *
 */

/* Include files */
#include "bfs_conncomp1_wet_initialize.h"
#include "_coder_bfs_conncomp1_wet_mex.h"
#include "bfs_conncomp1_wet_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void bfs_conncomp1_wet_initialize(void)
{
  mex_InitInfAndNan();
  mexFunctionCreateRootTLS();
  emlrtClearAllocCountR2012b(emlrtRootTLSGlobal, false, 0U, NULL);
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (bfs_conncomp1_wet_initialize.c) */
