/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * ppc_val2.h
 *
 * Code generation for function 'ppc_val2'
 *
 */

#pragma once

/* Include files */
#include "ppc_val2_types.h"
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void ppc_val2(const emlrtStack *sp, const real_T X[20],
              const emxArray_real_T *C, const emxArray_real_T *D,
              const emxArray_real_T *x, real_T d, emxArray_real_T *y,
              emxArray_real_T *z);

/* End of code generation (ppc_val2.h) */
