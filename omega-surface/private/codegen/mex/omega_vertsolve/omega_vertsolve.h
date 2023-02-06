/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * omega_vertsolve.h
 *
 * Code generation for function 'omega_vertsolve'
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
void omega_vertsolve(const emxArray_real_T *Sppc, const emxArray_real_T *Tppc,
                     const emxArray_real_T *P, const emxArray_real_T *BotK,
                     emxArray_real_T *s, emxArray_real_T *t, emxArray_real_T *p,
                     real_T tolp, const emxArray_real_T *phi);

/* End of code generation (omega_vertsolve.h) */
