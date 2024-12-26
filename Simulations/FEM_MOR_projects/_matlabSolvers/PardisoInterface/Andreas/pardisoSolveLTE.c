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
  MKL_INT mtype, error, maxfct, mnum, msglvl, phase;
  
    /* Local variables           */
  MKL_INT n_col;			/* number of columns of A */
  MKL_INT nnz_sym;  /* nonzeros in A */
  MKL_INT i, j;		/* loop counter */
  MKL_INT *perm = NULL;			/* nonzeros in A */
  MKL_INT Pnrhs = 1;		/* number if RHS vectors*/

  double *A_values_sym = NULL;	/* real symmetric numerical values of input matrix A */
  MKL_INT *A_ia_sym = NULL;	/* real symmetric IA of input matrix A */
  MKL_INT *A_ja_sym = NULL;	/* real symmetric JA of input matrix A */
  MKL_Complex16 *A_values_sym_c = NULL;	/* complex symmetric numerical values of input matrix A */
  
  mxArray *b_in;
  double *bR = NULL;		/* real numerical values of rhs b */
  double *bI = NULL;		/* imaginary numerical values of rhs b */
  MKL_Complex16 *bC = NULL;		/* complex numerical values of rhs b */
  
  mxArray *x_out;
  double *xR = NULL;		/* real numerical values of solution x */
  double *xI = NULL;		/* imaginary numerical values of solution x */
  MKL_Complex16 *xC = NULL;		/* complex numerical values of solution x */
      
  bool releaseMemory = false;

  MKL_INT *error_out;
  /*bool isreal = true;*/
  double ddummy;
  int idummy;

  
    /* === Check inputs ===================================================== */
  if (nrhs !=8 && nrhs!=9)
    mexErrMsgTxt("pardiso: incorrect number of input arguments.\nPlease specify (mtype,iparm,pt,A_val,A_ia,A_ja,ncol,B,[,ReleaseMemory])");

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
         
  b_in = (mxArray *) prhs[7];
  bR = (double *) mxGetPr(b_in);
  bI = (double *) mxGetPi(b_in);
  Pnrhs = (MKL_INT) ( mxGetM(b_in)>mxGetN(b_in) ? mxGetN(b_in) : mxGetM(b_in) );
     
  if(nrhs==9)releaseMemory = (bool)(mxGetScalar(prhs[8]));
     
  if (mxIsSparse (b_in))
    mexErrMsgTxt("pardiso: right hand side matrix B must be full, not sparse.");

  #if DEBUG
  printf("releaseMemory: %d\n",releaseMemory);
  printf("b is complex: %i\n",mxIsComplex(b_in));
  printf("sizeof(MKL_INT) %i\n",sizeof(MKL_INT));
  printf("sizeof(MKL_Complex16) %i\n",sizeof(MKL_Complex16));
  printf("sizeof(double) %i\n",sizeof(double));  
  printf("mxGetElementSize(A_ia_sym_in) %i\n",mxGetElementSize(A_ia_sym_in));
  printf("mxGetElementSize(A_values_in) %i\n",mxGetElementSize(A_values_in));
  printf("n_col %i\n",n_col);
  #endif
  
  if( mtype == 6 || mxIsComplex(b_in) )
  {
    #if DEBUG
      printf("allocating complex %ix%i solution vector\n",n_col, Pnrhs);
    #endif
    x_out = mxCreateDoubleMatrix(n_col, Pnrhs, mxCOMPLEX);
    xR = (double *) mxGetPr(x_out);
    xI = (double *) mxGetPi(x_out);  
  } else
  {
    #if DEBUG
      printf("allocating real %ix%i solution vector\n",n_col, Pnrhs);
    #endif
    x_out = mxCreateDoubleMatrix(n_col, Pnrhs, mxREAL);
    xR = (double *) mxGetPr(x_out);
  }    
  
  if (mtype == -2 || mtype == 2) 
  { 
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
    phase = 33;
    
    if(mxIsComplex(b_in))
    {
      PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
      &n_col, A_values_sym, A_ia_sym, A_ja_sym, perm, &Pnrhs,
      iparm, &msglvl, bR, xR, &error);
      PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
      &n_col, A_values_sym, A_ia_sym, A_ja_sym, perm, &Pnrhs,
      iparm, &msglvl, bI, xI, &error);
    } else      
      PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
      &n_col, A_values_sym, A_ia_sym, A_ja_sym, perm, &Pnrhs,
      iparm, &msglvl, bR, xR, &error);
  }
  else if (mtype == 6) {
    A_values_sym_c = (MKL_Complex16 *)mxGetPr(A_values_in);
    nnz_sym=mxGetM(A_values_in)/2;

    xC = mxCalloc(n_col*Pnrhs, sizeof(MKL_Complex16));   
    bC = mxCalloc(n_col*Pnrhs, sizeof(MKL_Complex16));
    
    for (i = 0; i < n_col*Pnrhs; i++) {
      bC[i].real = bR[i];
    }
    if(mxIsComplex(b_in))
    {
      for (i = 0; i < n_col*Pnrhs; i++) {
        bC[i].imag = bI[i];             
      }
    }
    else
    {
      #if DEBUG
        printf("zeroing imag part of b\n");
      #endif
      for (i = 0; i < n_col*Pnrhs; i++) {
        bC[i].imag = 0.0;
      }       
    }         
    #if DEBUG
    printf("sizeof(bC[0]) %d\n",sizeof(bC[0]));
    printf("sizeof(bC[0].real) %d\n",sizeof(bC[0].real));
    printf("sizeof(bC[0].imag) %d\n",sizeof(bC[0].imag));
    printf("sizeof(bR[0]) %d\n",sizeof(bR[0]));    
    printf("sizeof(bI[0]) %d\n",sizeof(bI[0]));
    printf("nnz_sym %i\n",nnz_sym);            
    for (i = 0; i <= n_col; i++)
      printf("A_ia_sym[%d] = %d\n", i, A_ia_sym[i]);
    for (i = 0; i < nnz_sym; i++)
      printf("i = %d, A_ja_sym[i] = %d\n", i, A_ja_sym[i]);
    for (i = 0; i < n_col; i++)
      for (j = A_ia_sym[i]-1; j < A_ia_sym[i+1]-1; j++)
        printf(" sym i=%d ja=%d  A.real=%e A.imag=%e\n", i+1, A_ja_sym[j], A_values_sym_c[j].real, A_values_sym_c[j].imag);
    for(j = 0; j < Pnrhs; j++)
      for (i = 0; i < n_col; i++)      
        printf("bC[%d,%d] = %e + %ej         bR = %e    bI = %e\n", i,j, bC[j*n_col+i].real, bC[j*n_col+i].imag,bR[j*n_col+i],bI[j*n_col+i]);
    #endif
    phase = 33;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
    &n_col, (double*)A_values_sym_c, A_ia_sym, A_ja_sym, perm, &Pnrhs,
    iparm, &msglvl, (double*)bC, (double*)xC, &error);
    for (i = 0; i < n_col*Pnrhs; i++) {
      xR[i] = xC[i].real;
      xI[i] = xC[i].imag;
    }
    mxFree(bC);
    mxFree(xC);
  }
    
  #if DEBUG
  printf(">>> PARDISO >>> Iterative refinements  %d \n", iparm[6]);
  printf(">>> PARDISO >>> Perturbed pivots  %d \n", iparm[13]);
  #endif
  
  if (releaseMemory) {
    /* -------------------------------------------------------------------- */
    /* .. Termination and release of memory. */
    /* -------------------------------------------------------------------- */
    phase = -1; /* Release internal memory. */
    /*if (isreal)
      PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
      &n_col, A_values_sym, A_ia_sym, A_ja_sym, perm, &Pnrhs,
      iparm, &msglvl, bR, xR, &error);
    else
      PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
      &n_col, (double*)A_values_sym_c, A_ia_sym, A_ja_sym, perm, &Pnrhs,
      iparm, &msglvl, (double*)bC, (double*)xC, &error);*/
    
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
    &n_col, &ddummy, A_ia_sym, A_ja_sym, &idummy, &Pnrhs,
    iparm, &msglvl, &ddummy, &ddummy, &error);
  }
  

  plhs[0] = x_out;
  plhs[1] = mxCreateDoubleMatrix(1, 1, mxINT32_CLASS);
  error_out = (MKL_INT*) mxGetPr(plhs[1]);
  *error_out = error;  
}


