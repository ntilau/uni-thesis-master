#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"
#include "mkl_types.h"
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
#include "buildSymmetricCSR_Structure.c"

extern MKL_INT PARDISO
(void *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *,
double *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *,
MKL_INT *, double *, double *, MKL_INT *);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Internal solver memory pointer pt, */
  /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
  /* or void *pt[64] should be OK on both architectures */
  double *pt;
  
  /* Pardiso control parameters. */
  MKL_INT mtype, phase, error, maxfct, mnum, msglvl;
  MKL_INT nrhsI = 1; /* Number of right hand sides. */
  
    /* Local variables           */
  MKL_INT full;			/* TRUE if input matrix full, FALSE if sparse */
  MKL_INT n_col;			/* number of columns of A */
  MKL_INT n_row;			/* number of rows of A */
  mwIndex *A_ja;			/* row indices of input matrix A */
  mwIndex *A_ia;			/* column pointers of input matrix A */
  double *A_valuesR;		/* real numerical values of input matrix A */
  double *A_valuesI;		/* imaginary numerical values of input matrix A */
  MKL_INT nnz;			/* nonzeros in A */
  
  MKL_INT i, j;		/* loop counter */
  
  MKL_INT idum; /* Integer dummy. */
  double ddum; /* Double dummy */
  
  /* Pardiso control parameters. */
  MKL_INT *iparm;
  /* Auxiliary variables. */
  double *iparm_out;
  
  double *pt_out;
  mxArray *A_input;
  
  MKL_INT nnz_sym;  /* nonzeros in A */
  double *A_values_sym = NULL;	/* real symmetric numerical values of input matrix A */
  MKL_INT *A_ia_sym = NULL;	/* real symmetric IA of input matrix A */
  MKL_INT *A_ja_sym = NULL;	/* real symmetric JA of input matrix A */
  MKL_Complex16 *A_values_sym_c = NULL;	/* complex symmetric numerical values of input matrix A */

  MKL_INT numZerosDiag;

  pt = (double *)calloc(64, sizeof(double));
  iparm = (MKL_INT *)calloc(64, sizeof(MKL_INT));
  mtype = (MKL_INT)(mxGetScalar(prhs[0]));
  A_input = (mxArray *)prhs[1];
    /* get size of input matrix A */
  n_col = (MKL_INT)mxGetN(A_input);
  n_row = (MKL_INT)mxGetM(A_input);
  nnz = (MKL_INT)mxGetNzmax(A_input);
  
  A_valuesR = (double *) mxGetPr(A_input);
  if (mtype == 6)
    A_valuesI = (double *) mxGetPi(A_input);
  
  if (n_col != n_row)
    mexErrMsgTxt ("pardiso: matrix must be square.");
  A_ja = mxGetIr(A_input);
  A_ia = mxGetJc(A_input);
  
  /* === If A is full, convert to a sparse matrix ========================= */
  if (mxGetNumberOfDimensions(A_input) != 2)
    mexErrMsgTxt("pardiso: input matrix must be 2-dimensional.");
  full = mxIsSparse (A_input);
  if (!full)
    mexErrMsgTxt("pardiso: input matrix must be MATLAB sparse.");
  
/* -------------------------------------------------------------------- */
/* .. Setup Pardiso control parameters. */
/* -------------------------------------------------------------------- */
  for (i = 0; i < 64; i++) {
    iparm[i] = 0;
  }
  iparm[0] = 0; /* No solver default: 1, use solver defaults: 0 */
  iparm[1] = 2; /* Fill-in reordering from METIS */
/*   /* Numbers of processors, value of OMP_NUM_THREADS */
  iparm[2] = 2;
  iparm[3] = 0; /* No iterative-direct algorithm */
  /*iparm[4] = 0; /* No user fill-in reducing permutation */
/*   iparm[5] = 0; /* Write solution into x */
/*   iparm[6] = 0; /* Not in use */
/*   iparm[7] = 2; /* Max numbers of iterative refinement steps */
/*   iparm[8] = 0; /* Not in use */
/*   iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
/*   iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
/*   iparm[11] = 0; /* Not in use */
/*   iparm[12] = 0; /* Not in use */
/*   iparm[13] = 0; /* Output: Number of perturbed pivots */
/*   iparm[14] = 0; /* Not in use */
/*   iparm[15] = 0; /* Not in use */
/*   iparm[16] = 0; /* Not in use */
/*   iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
/*   iparm[18] = -1; /* Output: Mflops for LU factorization */
/*   iparm[19] = 0; /* Output: Numbers of CG Iterations */
  maxfct = 1; /* Maximum number of numerical factorizations. */
  mnum = 1; /* Which factorization to use. */
  msglvl = 1; /* Print statistical information in file */
  error = 0; /* Initialize error flag */
/* -------------------------------------------------------------------- */
/* .. Initialize the internal solver memory pointer. This is only */
/* necessary for the FIRST call of the PARDISO solver. */
/* -------------------------------------------------------------------- */
  for (i = 0; i < 64; i++)
    pt[i] = 0;
  
  numZerosDiag = numbersOfZerosOnDiag(A_ja, A_ia, nnz, n_row);
  /* Pardiso explicitly wants values on main diagonal, even if zero */
  nnz_sym =  (nnz - (n_row - numZerosDiag)) / 2 + n_row;
    
  A_ja_sym = (MKL_INT *) calloc(nnz_sym, sizeof(MKL_INT));
  A_ia_sym = (MKL_INT *) calloc(n_row + 1, sizeof(MKL_INT));
  
  if (mtype == -2 || mtype == 2) {
    A_values_sym = (double *) calloc(nnz_sym, sizeof(double));
    buildSymmetricCSR_Structure(A_values_sym, A_ja_sym, A_ia_sym, n_col, A_valuesR, A_ja, A_ia);
    #if DEBUG
    printf("nnz_sym = %d\n", nnz_sym);
    for (i = 0; i <= n_col; i++)
      printf("A_ia_sym[%d] = %d\n", i, A_ia_sym[i]);
    for (i = 0; i < nnz_sym; i++)
      printf("i = %d, A_ja_sym[i] = %d\n", i, A_ja_sym[i]);
    for (i = 0; i < n_col; i++)
      for (j = A_ia_sym[i]-1; j < A_ia_sym[i+1]-1; j++)
        printf(" sym i=%d ja=%d  A.real=%e \n", i+1, A_ja_sym[j], A_values_sym[j]);
    #endif
    phase = 11;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
    &n_col, A_values_sym, A_ia_sym, A_ja_sym, &idum, &nrhsI,
    iparm, &msglvl, &ddum, &ddum, &error);
  }
  else if (mtype == 6) {
    A_values_sym_c = calloc(nnz_sym, sizeof(MKL_Complex16));
    buildSymmetricCSR_StructureComplex(A_values_sym_c, A_ja_sym, A_ia_sym, n_col, A_valuesR, A_valuesI, A_ja, A_ia);
    #if DEBUG
    for (i = 0; i <= n_col; i++)
      printf("A_ia_sym[%d] = %d\n", i, A_ia_sym[i]);
    for (i = 0; i < nnz_sym; i++)
      printf("i = %d, A_ja_sym[i] = %d\n", i, A_ja_sym[i]);
    for (i = 0; i < n_col; i++)
      for (j = A_ia_sym[i]-1; j < A_ia_sym[i+1]-1; j++)
        printf(" sym i=%d ja=%d  A.real=%e A.imag=%e\n", i+1, A_ja_sym[j], A_values_sym_c[j].real, A_values_sym_c[j].imag);
    #endif
    phase = 11;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
    &n_col, (double*)A_values_sym_c, A_ia_sym, A_ja_sym, &idum, &nrhsI,
    iparm, &msglvl, &ddum, &ddum, &error);
  }
  else
    mexErrMsgTxt("pardiso: unknown mtype");
  
  if (error != 0) {
    printf("ERROR during symbolic factorization: %d\n", error);
  }
  
  printf("\nReordering completed ... ");
  printf("\nNumber of nonzeros in factors = %d", iparm[17]);
  printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
  
  if ((iparm[15]+iparm[16]) > iparm[14])
    printf("\n>>> PARDISO >>> Total PARDISO Memory consumption is %d KBytes\n", iparm[15]+iparm[16]);
  else
    printf("\n>>> PARDISO >>> Total PARDISO Memory consumption is %d KBytes\n", iparm[14]);
  
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
  if (A_ja_sym != NULL) {
    free(A_ja_sym);
    A_ja_sym = NULL;
  }
  
  plhs[0] = mxCreateDoubleMatrix(1, 64, mxREAL);
  iparm_out = mxGetPr(plhs[0]);
  
  for (i = 0; i < 64; i++)
    iparm_out[i] = iparm[i];
  
  plhs[1] = mxCreateDoubleMatrix(1, 64, mxREAL);
  pt_out = mxGetPr(plhs[1]);
  
  for (i = 0; i < 64; i++)
    pt_out[i] = pt[i];
  
  free(pt);
  free(iparm);
}



