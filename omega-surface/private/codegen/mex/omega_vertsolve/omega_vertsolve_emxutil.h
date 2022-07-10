/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * omega_vertsolve_emxutil.h
 *
 * Code generation for function 'omega_vertsolve_emxutil'
 *
 */

#pragma once

/* Include files */
#include "omega_vertsolve_types.h"
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void emxEnsureCapacity_real_T(emxArray_real_T *emxArray, int32_T oldNumel);

void emxFree_real_T(emxArray_real_T **pEmxArray);

void emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions,
                    boolean_T doPush);

/* End of code generation (omega_vertsolve_emxutil.h) */
