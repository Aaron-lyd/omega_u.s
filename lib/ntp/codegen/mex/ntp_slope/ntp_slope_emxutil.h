/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * ntp_slope_emxutil.h
 *
 * Code generation for function 'ntp_slope_emxutil'
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
void emxEnsureCapacity_boolean_T(emxArray_boolean_T *emxArray,
                                 int32_T oldNumel);

void emxEnsureCapacity_real_T(emxArray_real_T *emxArray, int32_T oldNumel);

void emxFree_boolean_T(emxArray_boolean_T **pEmxArray);

void emxFree_real_T(emxArray_real_T **pEmxArray);

void emxInit_boolean_T(emxArray_boolean_T **pEmxArray, int32_T numDimensions,
                       boolean_T doPush);

void emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions,
                    boolean_T doPush);

/* End of code generation (ntp_slope_emxutil.h) */
