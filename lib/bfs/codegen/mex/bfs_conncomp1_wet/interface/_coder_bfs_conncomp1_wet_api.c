/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_bfs_conncomp1_wet_api.c
 *
 * Code generation for function '_coder_bfs_conncomp1_wet_api'
 *
 */

/* Include files */
#include "_coder_bfs_conncomp1_wet_api.h"
#include "bfs_conncomp1_wet.h"
#include "bfs_conncomp1_wet_data.h"
#include "bfs_conncomp1_wet_emxutil.h"
#include "bfs_conncomp1_wet_types.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static void b_emlrt_marshallIn(const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_real_T *y);

static const mxArray *b_emlrt_marshallOut(const real_T u);

static void c_emlrt_marshallIn(const mxArray *P, const char_T *identifier,
                               emxArray_real_T *y);

static void d_emlrt_marshallIn(const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_real_T *y);

static void e_emlrt_marshallIn(const mxArray *s, const char_T *identifier,
                               emxArray_real_T *y);

static void emlrt_marshallIn(const mxArray *Sppc, const char_T *identifier,
                             emxArray_real_T *y);

static void emlrt_marshallOut(const emxArray_real_T *u, const mxArray *y);

static void f_emlrt_marshallIn(const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_real_T *y);

static real_T g_emlrt_marshallIn(const mxArray *TOL_P,
                                 const char_T *identifier);

static real_T h_emlrt_marshallIn(const mxArray *u,
                                 const emlrtMsgIdentifier *parentId);

static void i_emlrt_marshallIn(const mxArray *qu, const char_T *identifier,
                               emxArray_real_T *y);

static void j_emlrt_marshallIn(const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_real_T *y);

static void k_emlrt_marshallIn(const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_real_T *ret);

static void l_emlrt_marshallIn(const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_real_T *ret);

static void m_emlrt_marshallIn(const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_real_T *ret);

static real_T n_emlrt_marshallIn(const mxArray *src,
                                 const emlrtMsgIdentifier *msgId);

static void o_emlrt_marshallIn(const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_real_T *ret);

/* Function Definitions */
static void b_emlrt_marshallIn(const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_real_T *y)
{
  k_emlrt_marshallIn(emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static const mxArray *b_emlrt_marshallOut(const real_T u)
{
  const mxArray *m;
  const mxArray *y;
  y = NULL;
  m = emlrtCreateDoubleScalar(u);
  emlrtAssign(&y, m);
  return y;
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

static void e_emlrt_marshallIn(const mxArray *s, const char_T *identifier,
                               emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  f_emlrt_marshallIn(emlrtAlias(s), &thisId, y);
  emlrtDestroyArray(&s);
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

static real_T g_emlrt_marshallIn(const mxArray *TOL_P, const char_T *identifier)
{
  emlrtMsgIdentifier thisId;
  real_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = h_emlrt_marshallIn(emlrtAlias(TOL_P), &thisId);
  emlrtDestroyArray(&TOL_P);
  return y;
}

static real_T h_emlrt_marshallIn(const mxArray *u,
                                 const emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = n_emlrt_marshallIn(emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static void i_emlrt_marshallIn(const mxArray *qu, const char_T *identifier,
                               emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  j_emlrt_marshallIn(emlrtAlias(qu), &thisId, y);
  emlrtDestroyArray(&qu);
}

static void j_emlrt_marshallIn(const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_real_T *y)
{
  o_emlrt_marshallIn(emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
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

static real_T n_emlrt_marshallIn(const mxArray *src,
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

static void o_emlrt_marshallIn(const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_real_T *ret)
{
  static const int32_T dims[2] = {16777216, 1};
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

void bfs_conncomp1_wet_api(const mxArray *const prhs[12], int32_T nlhs,
                           const mxArray *plhs[6])
{
  emxArray_real_T *A;
  emxArray_real_T *BotK;
  emxArray_real_T *ML;
  emxArray_real_T *P;
  emxArray_real_T *Sppc;
  emxArray_real_T *Tppc;
  emxArray_real_T *p;
  emxArray_real_T *qu;
  emxArray_real_T *s;
  emxArray_real_T *t;
  const mxArray *prhs_copy_idx_11;
  const mxArray *prhs_copy_idx_3;
  const mxArray *prhs_copy_idx_4;
  const mxArray *prhs_copy_idx_5;
  real_T TOL_P;
  real_T freshly_wet;
  real_T qt;
  real_T r;
  emlrtHeapReferenceStackEnterFcnR2012b(emlrtRootTLSGlobal);
  emxInit_real_T(&Sppc, 4, true);
  emxInit_real_T(&Tppc, 4, true);
  emxInit_real_T(&P, 3, true);
  emxInit_real_T(&s, 2, true);
  emxInit_real_T(&t, 2, true);
  emxInit_real_T(&p, 2, true);
  emxInit_real_T(&ML, 2, true);
  emxInit_real_T(&A, 3, true);
  emxInit_real_T(&BotK, 2, true);
  emxInit_real_T(&qu, 2, true);
  prhs_copy_idx_3 = emlrtProtectR2012b(prhs[3], 3, true, -1);
  prhs_copy_idx_4 = emlrtProtectR2012b(prhs[4], 4, true, -1);
  prhs_copy_idx_5 = emlrtProtectR2012b(prhs[5], 5, true, -1);
  prhs_copy_idx_11 = emlrtProtectR2012b(prhs[11], 11, true, -1);
  /* Marshall function inputs */
  Sppc->canFreeData = false;
  emlrt_marshallIn(emlrtAlias(prhs[0]), "Sppc", Sppc);
  Tppc->canFreeData = false;
  emlrt_marshallIn(emlrtAlias(prhs[1]), "Tppc", Tppc);
  P->canFreeData = false;
  c_emlrt_marshallIn(emlrtAlias(prhs[2]), "P", P);
  s->canFreeData = false;
  e_emlrt_marshallIn(emlrtAlias(prhs_copy_idx_3), "s", s);
  t->canFreeData = false;
  e_emlrt_marshallIn(emlrtAlias(prhs_copy_idx_4), "t", t);
  p->canFreeData = false;
  e_emlrt_marshallIn(emlrtAlias(prhs_copy_idx_5), "p", p);
  ML->canFreeData = false;
  e_emlrt_marshallIn(emlrtAlias(prhs[6]), "ML", ML);
  TOL_P = g_emlrt_marshallIn(emlrtAliasP(prhs[7]), "TOL_P");
  A->canFreeData = false;
  c_emlrt_marshallIn(emlrtAlias(prhs[8]), "A", A);
  BotK->canFreeData = false;
  e_emlrt_marshallIn(emlrtAlias(prhs[9]), "BotK", BotK);
  r = g_emlrt_marshallIn(emlrtAliasP(prhs[10]), "r");
  qu->canFreeData = false;
  i_emlrt_marshallIn(emlrtAlias(prhs_copy_idx_11), "qu", qu);
  /* Invoke the target function */
  bfs_conncomp1_wet(Sppc, Tppc, P, s, t, p, ML, TOL_P, A, BotK, r, qu,
                    &freshly_wet, &qt);
  /* Marshall function outputs */
  s->canFreeData = false;
  emlrt_marshallOut(s, prhs_copy_idx_3);
  plhs[0] = prhs_copy_idx_3;
  emxFree_real_T(&BotK);
  emxFree_real_T(&A);
  emxFree_real_T(&ML);
  emxFree_real_T(&s);
  emxFree_real_T(&P);
  emxFree_real_T(&Tppc);
  emxFree_real_T(&Sppc);
  if (nlhs > 1) {
    t->canFreeData = false;
    emlrt_marshallOut(t, prhs_copy_idx_4);
    plhs[1] = prhs_copy_idx_4;
  }
  emxFree_real_T(&t);
  if (nlhs > 2) {
    p->canFreeData = false;
    emlrt_marshallOut(p, prhs_copy_idx_5);
    plhs[2] = prhs_copy_idx_5;
  }
  emxFree_real_T(&p);
  if (nlhs > 3) {
    plhs[3] = b_emlrt_marshallOut(freshly_wet);
  }
  if (nlhs > 4) {
    qu->canFreeData = false;
    emlrt_marshallOut(qu, prhs_copy_idx_11);
    plhs[4] = prhs_copy_idx_11;
  }
  emxFree_real_T(&qu);
  if (nlhs > 5) {
    plhs[5] = b_emlrt_marshallOut(qt);
  }
  emlrtHeapReferenceStackLeaveFcnR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (_coder_bfs_conncomp1_wet_api.c) */
