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
#include "rtwtypes.h"
#include "tbs_vertsolve_types.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void ppc_val2(const emxArray_real_T *X, const real_T C_data[],
              const int32_T C_size[2], const real_T D_data[], real_T x,
              real_T *y, real_T *z);

/* End of code generation (ppc_val2.h) */
