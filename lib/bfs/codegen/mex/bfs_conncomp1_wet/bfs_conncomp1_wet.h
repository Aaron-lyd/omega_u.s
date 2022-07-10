/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * bfs_conncomp1_wet.h
 *
 * Code generation for function 'bfs_conncomp1_wet'
 *
 */

#pragma once

/* Include files */
#include "bfs_conncomp1_wet_types.h"
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void bfs_conncomp1_wet(const emxArray_real_T *Sppc, const emxArray_real_T *Tppc,
                       const emxArray_real_T *P, emxArray_real_T *s,
                       emxArray_real_T *t, emxArray_real_T *p,
                       const emxArray_real_T *ML, real_T TOL_P,
                       const emxArray_real_T *A, const emxArray_real_T *BotK,
                       real_T r, emxArray_real_T *qu, real_T *freshly_wet,
                       real_T *qt);

/* End of code generation (bfs_conncomp1_wet.h) */
