/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * ntp_slope.h
 *
 * Code generation for function 'ntp_slope'
 *
 */

#pragma once

/* Include files */
#include "ntp_slope_types.h"
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void ntp_slope(const emxArray_real_T *Sppc, const emxArray_real_T *Tppc,
               const emxArray_real_T *Z, const emxArray_real_T *z, real_T tolz,
               emxArray_real_T *dx, emxArray_real_T *dy, emxArray_real_T *sx,
               emxArray_real_T *sy);

/* End of code generation (ntp_slope.h) */
