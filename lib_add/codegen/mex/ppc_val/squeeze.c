/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * squeeze.c
 *
 * Code generation for function 'squeeze'
 *
 */

/* Include files */
#include "squeeze.h"
#include "ppc_val_data.h"
#include "ppc_val_emxutil.h"
#include "ppc_val_types.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"

/* Variable Definitions */
static emlrtRSInfo
    p_emlrtRSI =
        {
            32,        /* lineNo */
            "squeeze", /* fcnName */
            "/Applications/MATLAB_R2021a.app/toolbox/eml/lib/matlab/elmat/"
            "squeeze.m" /* pathName */
};

static emlrtRTEInfo
    o_emlrtRTEI =
        {
            32,        /* lineNo */
            1,         /* colNo */
            "squeeze", /* fName */
            "/Applications/MATLAB_R2021a.app/toolbox/eml/lib/matlab/elmat/"
            "squeeze.m" /* pName */
};

/* Function Definitions */
void squeeze(const emlrtStack *sp, const emxArray_real_T *a, emxArray_real_T *b)
{
  emlrtStack st;
  int32_T k;
  int32_T nx;
  int8_T sqsz[4];
  st.prev = sp;
  st.tls = sp->tls;
  sqsz[0] = 1;
  sqsz[1] = 1;
  sqsz[2] = 1;
  k = 4;
  while ((k > 2) && (a->size[k - 1] == 1)) {
    k--;
  }
  if ((k == 2) && (a->size[0] == 1)) {
    sqsz[1] = (int8_T)a->size[1];
  } else {
    k = 0;
    if (a->size[0] != 1) {
      k = 1;
      sqsz[0] = (int8_T)a->size[0];
    }
    if (a->size[1] != 1) {
      k++;
      sqsz[k - 1] = (int8_T)a->size[1];
    }
    if (a->size[2] != 1) {
      sqsz[k] = (int8_T)a->size[2];
    }
  }
  st.site = &p_emlrtRSI;
  nx = a->size[0] * a->size[1] * a->size[2];
  k = a->size[0];
  if (a->size[1] > a->size[0]) {
    k = a->size[1];
  }
  if (a->size[2] > k) {
    k = a->size[2];
  }
  if (1 > k) {
    k = 1;
  }
  k = muIntScalarMax_sint32(nx, k);
  if (sqsz[0] > k) {
    emlrtErrorWithMessageIdR2018a(&st, &c_emlrtRTEI,
                                  "Coder:toolbox:reshape_emptyReshapeLimit",
                                  "Coder:toolbox:reshape_emptyReshapeLimit", 0);
  }
  if (sqsz[1] > k) {
    emlrtErrorWithMessageIdR2018a(&st, &c_emlrtRTEI,
                                  "Coder:toolbox:reshape_emptyReshapeLimit",
                                  "Coder:toolbox:reshape_emptyReshapeLimit", 0);
  }
  if (sqsz[2] > k) {
    emlrtErrorWithMessageIdR2018a(&st, &c_emlrtRTEI,
                                  "Coder:toolbox:reshape_emptyReshapeLimit",
                                  "Coder:toolbox:reshape_emptyReshapeLimit", 0);
  }
  k = sqsz[0] * sqsz[1] * sqsz[2];
  if (k != nx) {
    emlrtErrorWithMessageIdR2018a(
        &st, &d_emlrtRTEI, "Coder:MATLAB:getReshapeDims_notSameNumel",
        "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
  }
  nx = b->size[0] * b->size[1] * b->size[2] * b->size[3];
  b->size[0] = sqsz[0];
  b->size[1] = sqsz[1];
  b->size[2] = sqsz[2];
  b->size[3] = 1;
  emxEnsureCapacity_real_T(sp, b, nx, &o_emlrtRTEI);
  for (nx = 0; nx < k; nx++) {
    b->data[nx] = a->data[nx];
  }
}

/* End of code generation (squeeze.c) */
