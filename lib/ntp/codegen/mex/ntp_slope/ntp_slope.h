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
               const emxArray_real_T *Z, const real_T z_data[],
               const int32_T z_size[2], real_T tolz, real_T dx_data[],
               int32_T dx_size[2], real_T dy_data[], int32_T dy_size[2],
               real_T sx_data[], int32_T sx_size[2], real_T sy_data[],
               int32_T sy_size[2]);

/* End of code generation (ntp_slope.h) */
