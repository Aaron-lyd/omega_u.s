/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * obs_vertsolve.h
 *
 * Code generation for function 'obs_vertsolve'
 *
 */

#pragma once

/* Include files */
#include "obs_vertsolve_types.h"
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void obs_vertsolve(const emxArray_real_T *Sppc, const emxArray_real_T *Tppc,
                   const emxArray_real_T *P, const emxArray_real_T *BotK,
                   emxArray_real_T *s, emxArray_real_T *t, emxArray_real_T *p,
                   const real_T dfnb_data[], const int32_T dfnb_size[2],
                   const real_T dfnc_data[], const int32_T dfnc_size[2],
                   real_T s0, real_T t0, real_T tolp);

/* End of code generation (obs_vertsolve.h) */
