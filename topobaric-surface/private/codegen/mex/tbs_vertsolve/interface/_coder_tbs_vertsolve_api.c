/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_tbs_vertsolve_api.c
 *
 * Code generation for function '_coder_tbs_vertsolve_api'
 *
 */

/* Include files */
#include "_coder_tbs_vertsolve_api.h"
#include "rt_nonfinite.h"
#include "tbs_vertsolve.h"
#include "tbs_vertsolve_data.h"
#include "tbs_vertsolve_emxutil.h"
#include "tbs_vertsolve_types.h"

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

static void g_emlrt_marshallIn(const mxArray *d_fn, const char_T *identifier,
                               emxArray_real_T *y);

static void h_emlrt_marshallIn(const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_real_T *y);

static real_T i_emlrt_marshallIn(const mxArray *s0, const char_T *identifier);

static real_T j_emlrt_marshallIn(const mxArray *u,
                                 const emlrtMsgIdentifier *parentId);

static void k_emlrt_marshallIn(const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_real_T *ret);

static void l_emlrt_marshallIn(const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_real_T *ret);

static void m_emlrt_marshallIn(const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_real_T *ret);

static void n_emlrt_marshallIn(const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_real_T *ret);

static real_T o_emlrt_marshallIn(const mxArray *src,
                                 const emlrtMsgIdentifier *msgId);

/* Function Definitions */
static void b_emlrt_marshallIn(const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_real_T *y)
{
  k_emlrt_marshallIn(emlrtAlias(u), parentId, y);
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
  l_emlrt_marshallIn(emlrtAlias(u), parentId, y);
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
  emlrtMxSetData((mxArray *)y, &u->data[0]);
  emlrtSetDimensions((mxArray *)y, &u->size[0], 2);
}

static void f_emlrt_marshallIn(const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_real_T *y)
{
  m_emlrt_marshallIn(emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void g_emlrt_marshallIn(const mxArray *d_fn, const char_T *identifier,
                               emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  h_emlrt_marshallIn(emlrtAlias(d_fn), &thisId, y);
  emlrtDestroyArray(&d_fn);
}

static void h_emlrt_marshallIn(const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_real_T *y)
{
  n_emlrt_marshallIn(emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static real_T i_emlrt_marshallIn(const mxArray *s0, const char_T *identifier)
{
  emlrtMsgIdentifier thisId;
  real_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = j_emlrt_marshallIn(emlrtAlias(s0), &thisId);
  emlrtDestroyArray(&s0);
  return y;
}

static real_T j_emlrt_marshallIn(const mxArray *u,
                                 const emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = o_emlrt_marshallIn(emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static void k_emlrt_marshallIn(const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_real_T *ret)
{
  static const int32_T dims[4] = {-1, -1, -1, -1};
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

static void l_emlrt_marshallIn(const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_real_T *ret)
{
  static const int32_T dims[3] = {-1, -1, -1};
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

static void m_emlrt_marshallIn(const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_real_T *ret)
{
  static const int32_T dims[2] = {4096, 4096};
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

static void n_emlrt_marshallIn(const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_real_T *ret)
{
  static const int32_T dims[2] = {5, 16777216};
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

static real_T o_emlrt_marshallIn(const mxArray *src,
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

void tbs_vertsolve_api(const mxArray *const prhs[13], int32_T nlhs,
                       const mxArray *plhs[3])
{
  emxArray_real_T *BotK;
  emxArray_real_T *P;
  emxArray_real_T *Sppc;
  emxArray_real_T *Tppc;
  emxArray_real_T *branchmap;
  emxArray_real_T *d_fn;
  emxArray_real_T *p;
  emxArray_real_T *s;
  emxArray_real_T *t;
  const mxArray *prhs_copy_idx_4;
  const mxArray *prhs_copy_idx_5;
  const mxArray *prhs_copy_idx_6;
  real_T DP;
  real_T s0;
  real_T t0;
  real_T tolp;
  emlrtHeapReferenceStackEnterFcnR2012b(emlrtRootTLSGlobal);
  emxInit_real_T(&Sppc, 4, true);
  emxInit_real_T(&Tppc, 4, true);
  emxInit_real_T(&P, 3, true);
  emxInit_real_T(&BotK, 2, true);
  emxInit_real_T(&s, 2, true);
  emxInit_real_T(&t, 2, true);
  emxInit_real_T(&p, 2, true);
  emxInit_real_T(&branchmap, 2, true);
  emxInit_real_T(&d_fn, 2, true);
  prhs_copy_idx_4 = emlrtProtectR2012b(prhs[4], 4, true, -1);
  prhs_copy_idx_5 = emlrtProtectR2012b(prhs[5], 5, true, -1);
  prhs_copy_idx_6 = emlrtProtectR2012b(prhs[6], 6, true, -1);
  /* Marshall function inputs */
  Sppc->canFreeData = false;
  emlrt_marshallIn(emlrtAlias(prhs[0]), "Sppc", Sppc);
  Tppc->canFreeData = false;
  emlrt_marshallIn(emlrtAlias(prhs[1]), "Tppc", Tppc);
  P->canFreeData = false;
  c_emlrt_marshallIn(emlrtAlias(prhs[2]), "P", P);
  BotK->canFreeData = false;
  e_emlrt_marshallIn(emlrtAlias(prhs[3]), "BotK", BotK);
  s->canFreeData = false;
  e_emlrt_marshallIn(emlrtAlias(prhs_copy_idx_4), "s", s);
  t->canFreeData = false;
  e_emlrt_marshallIn(emlrtAlias(prhs_copy_idx_5), "t", t);
  p->canFreeData = false;
  e_emlrt_marshallIn(emlrtAlias(prhs_copy_idx_6), "p", p);
  branchmap->canFreeData = false;
  e_emlrt_marshallIn(emlrtAlias(prhs[7]), "branchmap", branchmap);
  d_fn->canFreeData = false;
  g_emlrt_marshallIn(emlrtAlias(prhs[8]), "d_fn", d_fn);
  s0 = i_emlrt_marshallIn(emlrtAliasP(prhs[9]), "s0");
  t0 = i_emlrt_marshallIn(emlrtAliasP(prhs[10]), "t0");
  tolp = i_emlrt_marshallIn(emlrtAliasP(prhs[11]), "tolp");
  DP = i_emlrt_marshallIn(emlrtAliasP(prhs[12]), "DP");
  /* Invoke the target function */
  tbs_vertsolve(Sppc, Tppc, P, BotK, s, t, p, branchmap, d_fn, s0, t0, tolp,
                DP);
  /* Marshall function outputs */
  p->canFreeData = false;
  emlrt_marshallOut(p, prhs_copy_idx_6);
  plhs[0] = prhs_copy_idx_6;
  emxFree_real_T(&d_fn);
  emxFree_real_T(&branchmap);
  emxFree_real_T(&p);
  emxFree_real_T(&BotK);
  emxFree_real_T(&P);
  emxFree_real_T(&Tppc);
  emxFree_real_T(&Sppc);
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

/* End of code generation (_coder_tbs_vertsolve_api.c) */
