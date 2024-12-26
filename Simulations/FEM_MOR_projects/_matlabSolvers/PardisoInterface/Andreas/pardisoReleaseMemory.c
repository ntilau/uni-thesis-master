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
#include "buildSymmetricCSR_Structure.c"
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
  MKL_INT mtype, error, maxfct, mnum, msglvl, phase,Pnrhs;
  
    /* Local variables           */
  MKL_INT n_col;			/* number of columns of A */
  MKL_INT i, j;		/* loop counter */
  MKL_INT *perm = NULL;			/* nonzeros in A */

  double *A_values_sym = NULL;	/* real symmetric numerical values of input matrix A */
  MKL_Complex16 *A_values_sym_c = NULL;	/* complex symmetric numerical values of input matrix A */
  MKL_INT *A_ia_sym = NULL;	/* real symmetric IA of input matrix A */
  MKL_INT *A_ja_sym = NULL;	/* real symmetric JA of input matrix A */
  
  MKL_INT *error_out;
  double ddummy;
  int idummy;
  
  bool isreal = true;
  
    /* === Check inputs ===================================================== */
  if (nrhs !=7)
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
      
  maxfct = 1;
  mnum = 1;
  Pnrhs = 1;
  msglvl = 1;
  error = 0;
  
  #if DEBUG
  printf("sizeof(MKL_INT) %i\n",sizeof(MKL_INT));
  printf("sizeof(double) %i\n",sizeof(double));
  printf("mxGetElementSize(A_ia_sym_in) %i\n",mxGetElementSize(A_ia_sym_in));
  printf("mxGetElementSize(A_values_in) %i\n",mxGetElementSize(A_values_in));
  printf("n_col %i\n",n_col);
  #endif
  
  
  if (mtype == -2 || mtype == 2) 
  { 
    A_values_sym = (double *)mxGetPr(A_values_in);
  }
  else if (mtype == 6) {
    A_values_sym_c = (MKL_Complex16 *)mxGetPr(A_values_in);
  }
    
  /* -------------------------------------------------------------------- */
  /* .. Termination and release of memory. */
  /* -------------------------------------------------------------------- */
  phase = -1; /* Release internal memory. */   
  PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
  &n_col, &ddummy, A_ia_sym, A_ja_sym, &idummy, &Pnrhs,
  iparm, &msglvl, &ddummy, &ddummy, &error);
  
  
  plhs[0] = mxCreateDoubleMatrix(1, 1, mxINT32_CLASS);
  error_out = (MKL_INT*)mxGetPr(plhs[0]);
  *error_out = error;  
    
}


