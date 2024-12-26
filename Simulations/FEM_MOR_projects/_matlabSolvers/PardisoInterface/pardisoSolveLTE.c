#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"
#include "mkl_types.h"
#include "buildSymmetricCSR_Structure.c"
/* PARDISO prototype. */
#if defined(_WIN32) || defined(_WIN64)
#define pardiso_ PARDISO
#else
#define PARDISO pardiso_
#endif
#if defined(MKL_ILP64)
#define MKL_INT long long
#else
#define MKL_INT int
#endif
#define DEBUG 0

extern MKL_INT PARDISO
(void *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *,
double *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *,
MKL_INT *, double *, double *, MKL_INT *);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* === Input variables ================================================== */
  mxArray *pt_input;
  double *pt, *iparmD;
  
  mxArray *iparm_input;
  MKL_INT *iparm;
  
  mxArray *b_input;
  mxArray *x_output;
  mxArray *A_input;
  
    /* Pardiso control parameters. */
  MKL_INT mtype, error, maxfct, mnum, msglvl, phase;
  
    /* Local variables           */
  MKL_INT full;			/* TRUE if input matrix full, FALSE if sparse */
  MKL_INT n_col;			/* number of columns of A */
  MKL_INT n_row;			/* number of rows of A */
  mwIndex *A_ja = NULL;			/* row indices of input matrix A */
  mwIndex *A_ia = NULL;			/* column pointers of input matrix A */
  double *A_valuesR = NULL;		/* real numerical values of input matrix A */
  double *A_valuesI = NULL;		/* imaginary numerical values of input matrix A */
  
  double *bR = NULL;		/* real numerical values of rhs b */
  double *xR = NULL;		/* real numerical values of solution x */
  double *bI = NULL;		/* imaginary numerical values of rhs b */
  double *xI = NULL;		/* imaginary numerical values of solution x */
  MKL_Complex16 *bC = NULL;		/* complex numerical values of rhs b */
  MKL_Complex16 *xC = NULL;		/* complex numerical values of solution x */
  
  MKL_INT i, j;		/* loop counter */
  MKL_INT nnz, Pnrhs = 1;		/* nonzeros in A */
  MKL_INT *perm = NULL;			/* nonzeros in A */
  
  MKL_INT nnz_sym;
  double *A_values_sym = NULL;	/* real symmetric numerical values of input matrix A */
  MKL_Complex16 *A_values_sym_c = NULL;	/* complex symmetric numerical values of input matrix A */
  MKL_INT *A_ia_sym = NULL;	/* real symmetric IA of input matrix A */
  MKL_INT *A_ja_sym = NULL;	/* real symmetric JA of input matrix A */
  
  bool releaseMemory = false;
  MKL_INT numZerosDiag;
  
    /* === Check inputs ===================================================== */
  if (nlhs != 0 || nrhs != 7)
    mexErrMsgTxt("pardiso: incorrect number of input and/or output arguments.");
  
    /* === pardiso init ====================================================== */
  pt_input = (mxArray *) prhs[0];
  pt = (double *) mxGetPr(pt_input);
  
  mtype = (MKL_INT)(mxGetScalar(prhs[1]));
  
  A_input = (mxArray *) prhs[2];
    /* get size of input matrix A */
  n_col = (MKL_INT)mxGetN(A_input);
  n_row = (MKL_INT)mxGetM(A_input);
  nnz = (MKL_INT)mxGetNzmax(A_input);
  
  A_valuesR = (double *) mxGetPr(A_input);
  if (mtype == 6)
    A_valuesI = (double *) mxGetPi(A_input);
  A_ja	 = mxGetIr(A_input);
  A_ia	 = mxGetJc(A_input);
  
  iparm_input = (mxArray *) prhs[3];
  iparmD = (double *) mxGetData(iparm_input);
  iparm = (MKL_INT *) calloc(64, sizeof(MKL_INT));
  
  for (i = 0; i < 64; i++)
    iparm[i] = (MKL_INT) iparmD[i];
  
  b_input = (mxArray *) prhs[4];
  bR = (double *) mxGetPr(b_input);
  if (mtype == 6) {
    bI = (double *) mxGetPi(b_input);
    if (bI == NULL) {
      /* no imaginary part */
      bI = (double *) calloc(n_col, sizeof(double));
      for (i = 0; i < n_col; i++)
        bI[i] = 0.0;
    }
  }
   
  x_output = (mxArray *) prhs[5];
  xR = (double *) mxGetPr(x_output);
  if (mtype == 6)
    xI = (double *) mxGetPi(x_output);
  
  releaseMemory = (bool)(mxGetScalar(prhs[6]));
  
    /* === If A is full, convert to a sparse matrix ========================= */
  if (mxGetNumberOfDimensions(A_input) != 2)
    mexErrMsgTxt ("pardiso: input matrix must be 2-dimensional.");
  
  full = mxIsSparse(A_input);
  
  if (!full)
    mexErrMsgTxt ("pardiso: input matrix must be MATLAB sparse.");
  
  numZerosDiag = numbersOfZerosOnDiag(A_ja, A_ia, nnz, n_row);
  /* Pardiso explicitly wants values on main diagonal, even if zero */
  nnz_sym =  (nnz - (n_row - numZerosDiag)) / 2 + n_row;
  
  /* build symmetric structure */
  A_ja_sym = (MKL_INT *) calloc(nnz_sym, sizeof(MKL_INT));
  A_ia_sym = (MKL_INT *) calloc(n_row + 1, sizeof(MKL_INT));
  
  maxfct = 1;
  mnum = 1;
  Pnrhs = 1;
  msglvl = 1;
  error = 0;
  
  if (mtype == -2 || mtype == 2) {
    A_values_sym = (double *) calloc(nnz_sym, sizeof(double));
    buildSymmetricCSR_Structure(A_values_sym, A_ja_sym, A_ia_sym, n_col, A_valuesR, A_ja, A_ia);
    #if DEBUG
    for (i = 0; i <= n_col; i++)
      printf("A_ia_sym[%d] = %d\n", i, A_ia_sym[i]);
    for (i = 0; i < nnz_sym; i++)
      printf("i = %d, A_ja_sym[i] = %d\n", i, A_ja_sym[i]);
    for (i = 0; i < n_col; i++)
      for (j = A_ia_sym[i]-1; j < A_ia_sym[i+1]-1; j++)
        printf(" sym i=%d ja=%d  A.real=%e \n", i+1, A_ja_sym[j], A_values_sym[j]);
    #endif
    phase = 33;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
    &n_col, A_values_sym, A_ia_sym, A_ja_sym, perm, &Pnrhs,
    iparm, &msglvl, bR, xR, &error);
  }
  else if (mtype == 6) {
    A_values_sym_c = calloc(nnz_sym, sizeof(MKL_Complex16));
    xC = calloc(n_col, sizeof(MKL_Complex16));
    bC = calloc(n_col, sizeof(MKL_Complex16));
    for (i = 0; i < n_col; i++) {
      bC[i].real = bR[i];
      bC[i].imag = bI[i];
    }
    buildSymmetricCSR_StructureComplex(A_values_sym_c, A_ja_sym, A_ia_sym, n_col, A_valuesR, A_valuesI, A_ja, A_ia);
    #if DEBUG
    for (i = 0; i <= n_col; i++)
      printf("A_ia_sym[%d] = %d\n", i, A_ia_sym[i]);
    for (i = 0; i < nnz_sym; i++)
      printf("i = %d, A_ja_sym[i] = %d\n", i, A_ja_sym[i]);
    for (i = 0; i < n_col; i++)
      for (j = A_ia_sym[i]-1; j < A_ia_sym[i+1]-1; j++)
        printf(" sym i=%d ja=%d  A.real=%e A.imag=%e\n", i+1, A_ja_sym[j], A_values_sym_c[j].real, A_values_sym_c[j].imag);
    for (i = 0; i < n_col; i++)
      printf("bC[%d] = %e + %ej\n", i, bC[i].real, bC[i].imag);
    #endif
    phase = 33;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
    &n_col, (double*)A_values_sym_c, A_ia_sym, A_ja_sym, perm, &Pnrhs,
    iparm, &msglvl, (double*)bC, (double*)xC, &error);
    for (i = 0; i < n_col; i++) {
      xR[i] = xC[i].real;
      xI[i] = xC[i].imag;
    }
  }
  else
    mexErrMsgTxt("pardiso: unknown mtype");
  
  printf(">>> PARDISO >>> Iterative refinements  %d \n", iparm[6]);
  printf(">>> PARDISO >>> Perturbed pivots  %d \n", iparm[13]);
  
  if (releaseMemory) {
    /* -------------------------------------------------------------------- */
    /* .. Termination and release of memory. */
    /* -------------------------------------------------------------------- */
    phase = -1; /* Release internal memory. */
    if (mtype == -2 || mtype == 2)
      PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
      &n_col, A_values_sym, A_ia_sym, A_ja_sym, perm, &Pnrhs,
      iparm, &msglvl, bR, xR, &error);
    else if (mtype == 6)
      PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
      &n_col, (double*)A_values_sym_c, A_ia_sym, A_ja_sym, perm, &Pnrhs,
      iparm, &msglvl, (double*)bC, (double*)xC, &error);
    else
      mexErrMsgTxt("pardiso: unknown mtype");
  }
  
  for (i = 0 ; i < 64 ; i++)
    iparmD[i] = (MKL_INT)iparm [i];
  
  if (A_values_sym != NULL) {
    free(A_values_sym);
    A_values_sym = NULL;
  }
  if (A_values_sym_c != NULL) {
    free(A_values_sym_c);
    A_values_sym_c = NULL;
  }
  if (A_ia_sym != NULL) {
    free(A_ia_sym);
    A_ia_sym = NULL;
  }
  if (xC != NULL) {
    free(xC);
    xC = NULL;
  }
  if (bC != NULL) {
    free(bC);
    bC = NULL;
  }
  if (A_ja_sym != NULL) {
    free(A_ja_sym);
    A_ja_sym = NULL;
  }
  
  free(iparm);
}


