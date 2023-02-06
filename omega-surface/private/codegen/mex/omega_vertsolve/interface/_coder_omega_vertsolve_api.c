/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_omega_vertsolve_api.c
 *
 * Code generation for function '_coder_omega_vertsolve_api'
 *
 */

/* Include files */
#include "_coder_omega_vertsolve_api.h"
#include "omega_vertsolve.h"
#include "omega_vertsolve_data.h"
#include "omega_vertsolve_emxutil.h"
#include "omega_vertsolve_types.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static void b_emlrt_marshallIn(const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_real_T *y);

static void c_emlrt_marshallIn(const mxArray *P, const char_T *identifier,
                               emxArray_real_T *y);

static void d_emlrt_marshallIn(const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_real_T *y);

static void e_emlrt_marshallIn(const mxArray *BotK, const char_T *identifier,
                               emxArray_real_T *y);

static void emlrt_marshallIn(const mxArray *Sppc, const char_T *identifier,
                             emxArray_real_T *y);

static void emlrt_marshallOut(const emxArray_real_T *u, const mxArray *y);

static void f_emlrt_marshallIn(const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_real_T *y);

static real_T g_emlrt_marshallIn(const mxArray *tolp, const char_T *identifier);

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

static void c_emlrt_marshallIn(const mxArray *P, const char_T *identifier,
                               emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  d_emlrt_marshallIn(emlrtAlias(P), &thisId, y);
  emlrtDestroyArray(&P);
}

static void d_emlrt_marshallIn(const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_real_T *y)
{
  j_emlrt_marshallIn(emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void e_emlrt_marshallIn(const mxArray *BotK, const char_T *identifier,
                               emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  f_emlrt_marshallIn(emlrtAlias(BotK), &thisId, y);
  emlrtDestroyArray(&BotK);
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

static void emlrt_marshallOut(const emxArray_real_T *u, const mxArray *y)
{
  const real_T *u_data;
  u_data = u->data;
  emlrtMxSetData((mxArray *)y, (void *)&u_data[0]);
  emlrtSetDimensions((mxArray *)y, &u->size[0], 2);
}

static void f_emlrt_marshallIn(const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_real_T *y)
{
  k_emlrt_marshallIn(emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static real_T g_emlrt_marshallIn(const mxArray *tolp, const char_T *identifier)
{
  emlrtMsgIdentifier thisId;
  real_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = h_emlrt_marshallIn(emlrtAlias(tolp), &thisId);
  emlrtDestroyArray(&tolp);
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
  static const int32_T dims[4] = {-1, -1, -1, -1};
  int32_T iv[4];
  int32_T i;
  const boolean_T bv[4] = {true, true, true, true};
  emlrtCheckVsBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", false, 4U,
                            (const void *)&dims[0], &bv[0], &iv[0]);
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
  static const int32_T dims[3] = {-1, -1, -1};
  int32_T iv[3];
  int32_T i;
  const boolean_T bv[3] = {true, true, true};
  emlrtCheckVsBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", false, 3U,
                            (const void *)&dims[0], &bv[0], &iv[0]);
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
  static const int32_T dims[2] = {4096, 4096};
  int32_T iv[2];
  int32_T i;
  const boolean_T bv[2] = {true, true};
  emlrtCheckVsBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", false, 2U,
                            (const void *)&dims[0], &bv[0], &iv[0]);
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
  emlrtCheckBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", false, 0U,
                          (const void *)&dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

void omega_vertsolve_api(const mxArray *const prhs[9], int32_T nlhs,
                         const mxArray *plhs[3])
{
  emxArray_real_T *BotK;
  emxArray_real_T *P;
  emxArray_real_T *Sppc;
  emxArray_real_T *Tppc;
  emxArray_real_T *p;
  emxArray_real_T *phi;
  emxArray_real_T *s;
  emxArray_real_T *t;
  const mxArray *prhs_copy_idx_4;
  const mxArray *prhs_copy_idx_5;
  const mxArray *prhs_copy_idx_6;
  real_T tolp;
  emlrtHeapReferenceStackEnterFcnR2012b(emlrtRootTLSGlobal);
  prhs_copy_idx_4 = emlrtProtectR2012b(prhs[4], 4, true, -1);
  prhs_copy_idx_5 = emlrtProtectR2012b(prhs[5], 5, true, -1);
  prhs_copy_idx_6 = emlrtProtectR2012b(prhs[6], 6, true, -1);
  /* Marshall function inputs */
  emxInit_real_T(&Sppc, 4);
  Sppc->canFreeData = false;
  emlrt_marshallIn(emlrtAlias(prhs[0]), "Sppc", Sppc);
  emxInit_real_T(&Tppc, 4);
  Tppc->canFreeData = false;
  emlrt_marshallIn(emlrtAlias(prhs[1]), "Tppc", Tppc);
  emxInit_real_T(&P, 3);
  P->canFreeData = false;
  c_emlrt_marshallIn(emlrtAlias(prhs[2]), "P", P);
  emxInit_real_T(&BotK, 2);
  BotK->canFreeData = false;
  e_emlrt_marshallIn(emlrtAlias(prhs[3]), "BotK", BotK);
  emxInit_real_T(&s, 2);
  s->canFreeData = false;
  e_emlrt_marshallIn(emlrtAlias(prhs_copy_idx_4), "s", s);
  emxInit_real_T(&t, 2);
  t->canFreeData = false;
  e_emlrt_marshallIn(emlrtAlias(prhs_copy_idx_5), "t", t);
  emxInit_real_T(&p, 2);
  p->canFreeData = false;
  e_emlrt_marshallIn(emlrtAlias(prhs_copy_idx_6), "p", p);
  tolp = g_emlrt_marshallIn(emlrtAliasP(prhs[7]), "tolp");
  emxInit_real_T(&phi, 2);
  phi->canFreeData = false;
  e_emlrt_marshallIn(emlrtAlias(prhs[8]), "phi", phi);
  /* Invoke the target function */
  omega_vertsolve(Sppc, Tppc, P, BotK, s, t, p, tolp, phi);
  emxFree_real_T(&phi);
  emxFree_real_T(&BotK);
  emxFree_real_T(&P);
  emxFree_real_T(&Tppc);
  emxFree_real_T(&Sppc);
  /* Marshall function outputs */
  p->canFreeData = false;
  emlrt_marshallOut(p, prhs_copy_idx_6);
  emxFree_real_T(&p);
  plhs[0] = prhs_copy_idx_6;
  if (nlhs > 1) {
    s->canFreeData = false;
    emlrt_marshallOut(s, prhs_copy_idx_4);
    plhs[1] = prhs_copy_idx_4;
  }
  emxFree_real_T(&s);
  if (nlhs > 2) {
    t->canFreeData = false;
    emlrt_marshallOut(t, prhs_copy_idx_5);
    plhs[2] = prhs_copy_idx_5;
  }
  emxFree_real_T(&t);
  emlrtHeapReferenceStackLeaveFcnR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (_coder_omega_vertsolve_api.c) */
