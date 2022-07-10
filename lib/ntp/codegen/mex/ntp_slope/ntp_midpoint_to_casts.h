/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * ntp_midpoint_to_casts.h
 *
 * Code generation for function 'ntp_midpoint_to_casts'
 *
 */

#pragma once

/* Include files */
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
real_T
ntp_midpoint_to_casts(const real_T Sppc_A_data[], const int32_T Sppc_A_size[2],
                      const real_T Tppc_A_data[], const real_T P_A_data[],
                      int32_T P_A_size, real_T k_A, const real_T Sppc_B_data[],
                      const int32_T Sppc_B_size[2], const real_T Tppc_B_data[],
                      const real_T P_B_data[], int32_T P_B_size, real_T k_B,
                      real_T p_A, real_T p_B, real_T tolp);

/* End of code generation (ntp_midpoint_to_casts.h) */
