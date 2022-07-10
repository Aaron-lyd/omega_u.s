/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * tbs_vertsolve.h
 *
 * Code generation for function 'tbs_vertsolve'
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
void tbs_vertsolve(const emxArray_real_T *Sppc, const emxArray_real_T *Tppc,
                   const emxArray_real_T *P, const emxArray_real_T *BotK,
                   emxArray_real_T *s, emxArray_real_T *t, emxArray_real_T *p,
                   const emxArray_real_T *branchmap,
                   const emxArray_real_T *d_fn, real_T s0, real_T t0,
                   real_T tolp, real_T DP);

/* End of code generation (tbs_vertsolve.h) */
