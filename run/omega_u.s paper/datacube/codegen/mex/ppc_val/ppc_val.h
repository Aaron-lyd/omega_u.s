/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * ppc_val.h
 *
 * Code generation for function 'ppc_val'
 *
 */

#pragma once

/* Include files */
#include "ppc_val_types.h"
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void ppc_val(const emlrtStack *sp, const real_T X[20], const emxArray_real_T *C,
             const real_T x_data[], const int32_T x_size[2], real_T d,
             emxArray_real_T *y);

/* End of code generation (ppc_val.h) */
