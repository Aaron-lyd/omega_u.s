/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * bfs_conncomp1_wet_terminate.c
 *
 * Code generation for function 'bfs_conncomp1_wet_terminate'
 *
 */

/* Include files */
#include "bfs_conncomp1_wet_terminate.h"
#include "_coder_bfs_conncomp1_wet_mex.h"
#include "bfs_conncomp1_wet_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void bfs_conncomp1_wet_atexit(void)
{
  mexFunctionCreateRootTLS();
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void bfs_conncomp1_wet_terminate(void)
{
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (bfs_conncomp1_wet_terminate.c) */
