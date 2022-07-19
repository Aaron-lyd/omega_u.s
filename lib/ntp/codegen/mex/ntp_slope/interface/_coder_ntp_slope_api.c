/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_ntp_slope_api.c
 *
 * Code generation for function '_coder_ntp_slope_api'
 *
 */

/* Include files */
#include "_coder_ntp_slope_api.h"
#include "ntp_slope.h"
#include "ntp_slope_data.h"
#include "ntp_slope_emxutil.h"
#include "ntp_slope_types.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static void b_emlrt_marshallIn(const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_real_T *y);

static void c_emlrt_marshallIn(const mxArray *Z, const char_T *identifier,
                               emxArray_real_T *y);

static void d_emlrt_marshallIn(const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_real_T *y);

static void e_emlrt_marshallIn(const mxArray *z, const char_T *identifier,
                               emxArray_real_T *y);

static void emlrt_marshallIn(const mxArray *Sppc, const char_T *identifier,
                             emxArray_real_T *y);

static const mxArray *emlrt_marshallOut(const emxArray_real_T *u);

static void f_emlrt_marshallIn(const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_real_T *y);

static real_T g_emlrt_marshallIn(const mxArray *tolz, const char_T *identifier);

static real_T h_emlrt_marshallIn(const mxArray *u,
                                 const emlrtMsgIdentifier *parentId);

static void i_emlrt_marshallIn(const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_real_T *ret);

static void j_emlrt_marshallIn(const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_real_T *ret);

static void k_emlrt_marshallIn(const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_real_T *ret);

static real_T l_emlrt_marshallIn(const mxArray *src,
                                 const emlrtMsgIdentifier *msgId);

/* Function Definitions */
static void b_emlrt_marshallIn(const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_real_T *y)
{
  i_emlrt_marshallIn(emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void c_emlrt_marshallIn(const mxArray *Z, const char_T *identifier,
                               emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  d_emlrt_marshallIn(emlrtAlias(Z), &thisId, y);
  emlrtDestroyArray(&Z);
}

static void d_emlrt_marshallIn(const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_real_T *y)
{
  j_emlrt_marshallIn(emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void e_emlrt_marshallIn(const mxArray *z, const char_T *identifier,
                               emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  f_emlrt_marshallIn(emlrtAlias(z), &thisId, y);
  emlrtDestroyArray(&z);
}

static void emlrt_marshallIn(const mxArray *Sppc, const char_T *identifier,
                             emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  b_emlrt_marshallIn(emlrtAlias(Sppc), &thisId, y);
  emlrtDestroyArray(&Sppc);
}

static const mxArray *emlrt_marshallOut(const emxArray_real_T *u)
{
  static const int32_T iv[2] = {0, 0};
  const mxArray *m;
  const mxArray *y;
  y = NULL;
  m = emlrtCreateNumericArray(2, (const void *)&iv[0], mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, &u->data[0]);
  emlrtSetDimensions((mxArray *)m, &u->size[0], 2);
  emlrtAssign(&y, m);
  return y;
}

static void f_emlrt_marshallIn(const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_real_T *y)
{
  k_emlrt_marshallIn(emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static real_T g_emlrt_marshallIn(const mxArray *tolz, const char_T *identifier)
{
  emlrtMsgIdentifier thisId;
  real_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = h_emlrt_marshallIn(emlrtAlias(tolz), &thisId);
  emlrtDestroyArray(&tolz);
  return y;
}

static real_T h_emlrt_marshallIn(const mxArray *u,
                                 const emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = l_emlrt_marshallIn(emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static void i_emlrt_marshallIn(const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_real_T *ret)
{
  static const int32_T dims[4] = {8, 19, 257, 257};
  int32_T iv[4];
  int32_T i;
  const boolean_T bv[4] = {true, true, true, true};
  emlrtCheckVsBuiltInR2012b(emlrtRootTLSGlobal, msgId, src,
                            (const char_T *)"double", false, 4U,
                            (void *)&dims[0], &bv[0], &iv[0]);
  ret->allocatedSize = iv[0] * iv[1] * iv[2] * iv[3];
  i = ret->size[0] * ret->size[1] * ret->size[2] * ret->size[3];
  ret->size[0] = iv[0];
  ret->size[1] = iv[1];
  ret->size[2] = iv[2];
  ret->size[3] = iv[3];
  emxEnsureCapacity_real_T(ret, i);
  ret->data = (real_T *)emlrtMxGetData(src);
  ret->canFreeData = false;
  emlrtDestroyArray(&src);
}

static void j_emlrt_marshallIn(const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_real_T *ret)
{
  static const int32_T dims[3] = {20, 257, 257};
  int32_T iv[3];
  int32_T i;
  const boolean_T bv[3] = {true, true, true};
  emlrtCheckVsBuiltInR2012b(emlrtRootTLSGlobal, msgId, src,
                            (const char_T *)"double", false, 3U,
                            (void *)&dims[0], &bv[0], &iv[0]);
  ret->allocatedSize = iv[0] * iv[1] * iv[2];
  i = ret->size[0] * ret->size[1] * ret->size[2];
  ret->size[0] = iv[0];
  ret->size[1] = iv[1];
  ret->size[2] = iv[2];
  emxEnsureCapacity_real_T(ret, i);
  ret->data = (real_T *)emlrtMxGetData(src);
  ret->canFreeData = false;
  emlrtDestroyArray(&src);
}

static void k_emlrt_marshallIn(const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_real_T *ret)
{
  static const int32_T dims[2] = {257, 257};
  int32_T iv[2];
  int32_T i;
  const boolean_T bv[2] = {true, true};
  emlrtCheckVsBuiltInR2012b(emlrtRootTLSGlobal, msgId, src,
                            (const char_T *)"double", false, 2U,
                            (void *)&dims[0], &bv[0], &iv[0]);
  ret->allocatedSize = iv[0] * iv[1];
  i = ret->size[0] * ret->size[1];
  ret->size[0] = iv[0];
  ret->size[1] = iv[1];
  emxEnsureCapacity_real_T(ret, i);
  ret->data = (real_T *)emlrtMxGetData(src);
  ret->canFreeData = false;
  emlrtDestroyArray(&src);
}

static real_T l_emlrt_marshallIn(const mxArray *src,
                                 const emlrtMsgIdentifier *msgId)
{
  static const int32_T dims = 0;
  real_T ret;
  emlrtCheckBuiltInR2012b(emlrtRootTLSGlobal, msgId, src,
                          (const char_T *)"double", false, 0U, (void *)&dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

void ntp_slope_api(const mxArray *const prhs[7], int32_T nlhs,
                   const mxArray *plhs[2])
{
  emxArray_real_T *Sppc;
  emxArray_real_T *Tppc;
  emxArray_real_T *Z;
  emxArray_real_T *dx;
  emxArray_real_T *dy;
  emxArray_real_T *sx;
  emxArray_real_T *sy;
  emxArray_real_T *z;
  const mxArray *prhs_copy_idx_5;
  const mxArray *prhs_copy_idx_6;
  real_T tolz;
  emlrtHeapReferenceStackEnterFcnR2012b(emlrtRootTLSGlobal);
  emxInit_real_T(&Sppc, 4, true);
  emxInit_real_T(&Tppc, 4, true);
  emxInit_real_T(&Z, 3, true);
  emxInit_real_T(&z, 2, true);
  emxInit_real_T(&dx, 2, true);
  emxInit_real_T(&dy, 2, true);
  emxInit_real_T(&sx, 2, true);
  emxInit_real_T(&sy, 2, true);
  prhs_copy_idx_5 = emlrtProtectR2012b(prhs[5], 5, false, -1);
  prhs_copy_idx_6 = emlrtProtectR2012b(prhs[6], 6, false, -1);
  /* Marshall function inputs */
  Sppc->canFreeData = false;
  emlrt_marshallIn(emlrtAlias(prhs[0]), "Sppc", Sppc);
  Tppc->canFreeData = false;
  emlrt_marshallIn(emlrtAlias(prhs[1]), "Tppc", Tppc);
  Z->canFreeData = false;
  c_emlrt_marshallIn(emlrtAlias(prhs[2]), "Z", Z);
  z->canFreeData = false;
  e_emlrt_marshallIn(emlrtAlias(prhs[3]), "z", z);
  tolz = g_emlrt_marshallIn(emlrtAliasP(prhs[4]), "tolz");
  dx->canFreeData = false;
  e_emlrt_marshallIn(emlrtAlias(prhs_copy_idx_5), "dx", dx);
  dy->canFreeData = false;
  e_emlrt_marshallIn(emlrtAlias(prhs_copy_idx_6), "dy", dy);
  /* Invoke the target function */
  ntp_slope(Sppc, Tppc, Z, z, tolz, dx, dy, sx, sy);
  /* Marshall function outputs */
  sx->canFreeData = false;
  plhs[0] = emlrt_marshallOut(sx);
  emxFree_real_T(&sx);
  emxFree_real_T(&dy);
  emxFree_real_T(&dx);
  emxFree_real_T(&z);
  emxFree_real_T(&Z);
  emxFree_real_T(&Tppc);
  emxFree_real_T(&Sppc);
  if (nlhs > 1) {
    sy->canFreeData = false;
    plhs[1] = emlrt_marshallOut(sy);
  }
  emxFree_real_T(&sy);
  emlrtHeapReferenceStackLeaveFcnR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (_coder_ntp_slope_api.c) */
