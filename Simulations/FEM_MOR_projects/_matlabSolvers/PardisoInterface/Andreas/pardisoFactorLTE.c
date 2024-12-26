/* mtype  2 real and symmetric positive definite
 *       -2 real and symmetric indefinite
 *        6 complex and symmetric
 * 
 * Error Information
 *   0   No error.
 *  -1   Input inconsistent.
 *  -2   Not enough memory.
 *  -3   Reordering problem.
 *  -4   Zero pivot, numerical fact. or iterative refinement problem.
 *  -5   Unclassified (internal) error.
 *  -6   Preordering failed (matrix types 11, 13 only).
 *  -7   Diagonal matrix problem.
 *  -8   32-bit integer overflow problem
 * 
 */


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
#define mx_MKL_INT_CLASS mxINT64_CLASS
#else
#define MKL_INT int
#define mx_MKL_INT_CLASS mxINT32_CLASS
#endif

#define DEBUG 0
#include "buildSymmetricCSR_Structure.c"

extern MKL_INT PARDISO
(void *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *,
double *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *,
MKL_INT *, double *, double *, MKL_INT *);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* === Input variables ================================================== */
  mxArray *pt_in;
  double *pt;
  
  mxArray *iparm_in;
  MKL_INT	*iparm;
  
  mxArray *A_values_in, *A_ia_sym_in, *A_ja_sym_in ;
  
    /* Pardiso control parameters. */
  MKL_INT mtype, phase, error, maxfct, mnum, msglvl;
  
    /* Local variables           */
  MKL_INT n_col;			/* number of columns of A */
  MKL_INT nnz_sym;  /* nonzeros in A */
  MKL_INT i, j;		/* loop counter */
  MKL_INT *perm = NULL;			/* nonzeros in A */

  double *A_values_sym = NULL;	/* real symmetric numerical values of input matrix A */
  MKL_INT *A_ia_sym = NULL;	/* real symmetric IA of input matrix A */
  MKL_INT *A_ja_sym = NULL;	/* real symmetric JA of input matrix A */
  MKL_Complex16 *A_values_sym_c = NULL;	/* complex symmetric numerical values of input matrix A */
  
  MKL_INT nrhsI = 1;

  MKL_INT *error_out;
  
    /* === Check inputs ===================================================== */
  if (nrhs != 7)
    mexErrMsgTxt("pardiso: incorrect number of input arguments.\nPlease specify (mtype,iparm,pt,A_val,A_ia,A_ja,ncol)");

  mtype = (MKL_INT)(mxGetScalar(prhs[0]));    

  switch(mtype)
  {
    case -2:
    case 2:
    case 6:
      break;
    default:
      mexErrMsgTxt("pardiso: unknown mtype, only -2, 2 and 6 are implemented");
  }  
  /* pt , mtype , iparm , A_values , A_ia_sym , A_ja_sym , ncol */
    /* === pardiso init ====================================================== */
  iparm_in = (mxArray *) prhs[1];
  iparm=(MKL_INT *) mxGetPr(iparm_in);

  pt_in = (mxArray *) prhs[2];
  pt = (double *) mxGetPr(pt_in);    
  
  maxfct = 1;
  mnum = 1;
  msglvl = 1;
  error = 0;
  
  A_values_in = (mxArray *)prhs[3]; 
  A_ia_sym_in = (mxArray *)prhs[4];
  A_ia_sym = (MKL_INT *) mxGetPr(A_ia_sym_in);
  A_ja_sym_in =  (mxArray *)prhs[5];    
  A_ja_sym = (MKL_INT *) mxGetPr(A_ja_sym_in);
  
  n_col =(MKL_INT) mxGetScalar(prhs[6]);  

  #if DEBUG
  printf("sizeof(MKL_INT) %i\n",sizeof(MKL_INT));
  printf("sizeof(double) %i\n",sizeof(double));
  printf("mxGetElementSize(A_ia_sym_in) %i\n",mxGetElementSize(A_ia_sym_in));
  printf("mxGetElementSize(A_values_in) %i\n",mxGetElementSize(A_values_in));
  printf("n_col %i\n",n_col);
  #endif

  if (mtype == -2 || mtype == 2) {
    A_values_sym = (double *)mxGetPr(A_values_in);
    nnz_sym=mxGetM(A_values_in);

    #if DEBUG
    printf("nnz_sym %i\n",nnz_sym);
    for (i = 0; i <= n_col; i++)
      printf("A_ia_sym[%d] = %d\n", i, A_ia_sym[i]);
    for (i = 0; i < nnz_sym; i++)
      printf("i = %d, A_ja_sym[i] = %d\n", i, A_ja_sym[i]);
    for (i = 0; i < n_col; i++)
      for (j = A_ia_sym[i]-1; j < A_ia_sym[i+1]-1; j++)
        printf(" sym i=%d ja=%d  A.real=%e \n", i+1, A_ja_sym[j], A_values_sym[j]);
    #endif

    phase = 22;

    PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
    &n_col, A_values_sym, A_ia_sym, A_ja_sym, perm, &nrhsI,
    iparm, &msglvl, NULL, NULL, &error);

  }
  else if (mtype == 6) {
    A_values_sym_c = (MKL_Complex16 *)mxGetPr(A_values_in);
    nnz_sym=mxGetM(A_values_in)/2;

    #if DEBUG
    printf("nnz_sym %i\n",nnz_sym);
            
    for (i = 0; i < nnz_sym; i++)
      printf("A_ja_sym[%d] = %d\n", i, A_ja_sym[i]);
    for (i = 0; i < n_col; i++)
      for (j = A_ia_sym[i]-1; j < A_ia_sym[i+1]-1; j++)
        printf(" sym i=%d ja=%d  A.real=%e A.imag=%e\n", i+1, A_ja_sym[j], A_values_sym_c[j].real, A_values_sym_c[j].imag);
    #endif

    phase = 22;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
    &n_col, (double*)A_values_sym_c, A_ia_sym, A_ja_sym, perm, &nrhsI,
    iparm, &msglvl, NULL, NULL, &error);
  }
  
  #if DEBUG
  if ((iparm[15]+iparm[16]) > iparm[14])
    printf(">>> PARDISO >>> Total PARDISO Memory consumption is %d KBytes\n", iparm[15]+iparm[16]);
  else
    printf(">>> PARDISO >>> Total PARDISO Memory consumption is %d KBytes\n", iparm[14] );
  printf(">>> PARDISO >>> Postive Eigenvalues of matrix A:  %d\n", iparm[21]);
  printf(">>> PARDISO >>> Negative Eigenvalues of matrix A:  %d\n", n_col - iparm[21]);
  #endif

  plhs[0] =  mxCreateNumericMatrix(1,1, mx_MKL_INT_CLASS, mxREAL);
  error_out = (MKL_INT*)mxGetPr(plhs[0]);
  *error_out = error;
}

