/* ========================================================================== */
/* === pardiso mexFunction ================================================== */
/* ========================================================================== */

/*
    Usage:

    For compilation please look to the README file.
 
    Computes the solution for symmetric indefinite systems with
    PARDISO solver version 3.3
    
    Example:

    % for solving indefinite complex symmetric systems 
    mtype = -2 ;
    [iparm, pt] = pardisoinit(mtype);
    
    % iparm(3) must always be equal to the environment variable OMP_NUM_THREADS
    % only one processor here
    iparm (3) = 1; 

    % reorder the system
    pardisoreorder(pt, mtype, A, iparm);

    % factorize 
    pardisofactor(pt, mtype, A, iparm);

    % get numbers of equations and allocate memory for solution vector x
    % b is the right hand side
    neqns = size(A);
    x = ones(neqns(1),1);
    % create right hand side b
    b =  A*x;
    pardisosolve(pt, mtype, A, iparm, b, x);

    % check residual
    norm(A*x-b)/norm(b)

    Authors:

	The author of the code itself is Olaf Schenk, University of
	Basel. The solver was developed in collaboration with Klaus
	Gaertner, WIAS, Germany.

    Date:

        December  04, 2008. Update to Version 7.6.0.324 (R2008a)
        December  19, 2007. PARDISO Version 3.3
        September 26, 2005. PARDISO Version 2.1

    Acknowledgements:

	This work was supported by the Swiss Comission for Technology
	and Innovation (CTI) under grant number 7036 ENS-ES.

    Notice:

        Copyright (c)  2009, 2008, 2007, 2006, 2005, 2004  University of Basel.
        All Rights Reserved.

	THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
	EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.

    Availability:

	This file is located at

        http://www.pardiso-project.org
*/

/* ========================================================================== */
/* === Include files and prototypes ========================================= */
/* ========================================================================== */

#include "mex.h"
#include "matrix.h"
#include <stdlib.h>

#define DEBUG 0

/* PARDISO prototype. */

extern  int pardiso_
        (void *, int *, int *, int *, int *, int *,
         double *, int *, int *, int *, int *, int *,
         int *, double *, double *, int *);



static void pardiso_help (void) ;

/* ========================================================================== */
/* === pardiso mexFunction ================================================== */
/* ========================================================================== */

void mexFunction
(
    /* === Parameters ======================================================= */

    int nlhs,			/* number of left-hand sides */
    mxArray *plhs [],		/* left-hand side matrices */
    int nrhs,			/* number of right--hand sides */
    const mxArray *prhs []	/* right-hand side matrices */
)
{
    /* === Input variables ================================================== */
    mxArray *pt_input;
    double    *pt, *iparmD;
    
    mxArray *iparm_input;
    int	    *iparm;

    mxArray *b_input;
    mxArray *A_input ;
    
    /* Pardiso control parameters. */
    int      mtype, phase, error, maxfct, mnum, msglvl;

    /* Local variables           */
    int full ;			/* TRUE if input matrix full, FALSE if sparse */
    int n_col ;			/* number of columns of A */
    int n_row ;			/* number of rows of A */   
    mwIndex *A_ja ;			/* row indices of input matrix A */
    mwIndex *A_ia ;			/* column pointers of input matrix A */
   
    double *A_valuesR ;		/* real numerical values of input matrix A */  
    


    int j, k, i, ptr ;		/* loop counter */
    int nnz;			/* nonzeros in A */
    int *perm;			/* nonzeros in A */

    /* === Check inputs ===================================================== */

    if (nlhs != 0 || nrhs != 4)
    {
    	pardiso_help () ;
    	mexErrMsgTxt (
	    "pardiso: incorrect number of input and/or output arguments.") ;
    }

    /* === pardiso init ====================================================== */
    pt_input  = (mxArray *) prhs [0];
    pt        = (double *) mxGetPr(pt_input);
   
    mtype = (int)(mxGetScalar(prhs[1]));
    
    A_input = (mxArray *) prhs [2] ;
    /* get size of input matrix A */
    n_col = mxGetN (A_input) ;
    n_row = mxGetM (A_input) ;
    nnz = mxGetNzmax(A_input);

    A_valuesR        = (double *) mxGetPr(A_input);
 
   
#if DEBUG
    for (i = 0 ; i < nnz ; i++) 
        printf(" matrix A[%d] real: %e\n", i,  A_valuesR[i]);
#endif
  

    iparm_input  = (mxArray *) prhs [3];

    iparmD        = (double *) mxGetData(iparm_input);
    
    iparm = (int *) mxCalloc (64, sizeof (int)) ;
  
    for (i = 0 ; i < 64 ; i++)
    {
    	iparm[i] = (int) iparmD [i];
    }

	
    /* === If A is full, convert to a sparse matrix ========================= */

    if (mxGetNumberOfDimensions (A_input) != 2)
    {
    	pardiso_help () ;
    	mexErrMsgTxt ("pardiso: input matrix must be 2-dimensional.") ;
    }
   
    full = mxIsSparse (A_input) ;

    if (!full)
    {
	mexErrMsgTxt ("pardiso: input matrix must be MATLAB sparse.") ;
    }


    /* === Allocate workspace for construction of (A, IA, JA)  =========== */

    /* get size of input matrix A */
    n_col = mxGetN (A_input) ;
    n_row = mxGetM (A_input) ;

    if (n_col != n_row)
    {
    	pardiso_help () ;
    	mexErrMsgTxt ("pardiso: matrix must be square.") ;
    }

    A_ja	 = mxGetIr (A_input) ;
    A_ia	 = mxGetJc (A_input) ;
  
    nnz = mxGetNzmax(A_input);

    /* -------------------------------------------------------------------- */
    /* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
    /*     notation.                                                        */
    /* -------------------------------------------------------------------- */

    maxfct = 1;
    mnum = 1;
    nrhs = 1;
    msglvl = 0;
    error = 0;

#if DEBUG
 
    for (i = 0 ; i < n_col ; i++) 
	for (j = A_ia[i] ; j < A_ia[i+1] ; j++) 
	    printf("i=%d ja=%d  A.real=%e \n", i+1,  A_ja[j]+1, A_valuesR[j]);
   
#endif
    
   
    {
	int nnz_sym;
	int k, kk, i, j;
	double  *A_values_sym;	/* real symmetric numerical values of input matrix A */  
	int     *A_ia_sym;	/* real symmetric IA of input matrix A */  
	int     *A_ja_sym;	/* real symmetric JA of input matrix A */  

	nnz_sym =  (nnz - n_row)/2 + n_row;

	/* build symmetric structure */
	A_values_sym = (double *) mxCalloc (nnz_sym, sizeof (double)) ;
	A_ja_sym = (int *) mxCalloc (nnz_sym, sizeof (int)) ;
	A_ia_sym = (int *) mxCalloc (n_row, sizeof (int)) ;
	/* build symmetric CSR structure */
       
	k = 0;
	A_ia_sym[0] = 1;
	for (i = 0 ; i < n_col ; i++) 
	{
	    kk = 0;
	    for (j = A_ia[i] ; j < A_ia[i+1] ; j++) 
	    {
		/* upper part */
		if (A_ja[j] >= i)
		{
			    A_ja_sym[k] = A_ja[j]+1;
			    A_values_sym [k] = A_valuesR[j];
			    k +=1;
			    kk +=1;
			}
		    }
		    A_ia_sym[i+1] =  A_ia_sym[i] + kk;
		}
	#if DEBUG
		for (i = 0 ; i < n_col ; i++) 
		    for (j = A_ia_sym[i]-1 ; j < A_ia_sym[i+1]-1 ; j++) 
			printf(" sym i=%d ja=%d  A.real=%e \n", i+1,  A_ja_sym[j], A_values_sym[j]);
#endif

	phase = 11;
	pardiso_ (pt, &maxfct, &mnum, &mtype, &phase,
		  &n_col, A_values_sym, A_ia_sym, A_ja_sym, perm, &nrhs,
		  iparm, &msglvl, NULL, NULL, &error);


	if ( (iparm[15]+iparm[16]) > iparm[14] ) 
	{
	    printf(">>> PARDISO >>> Total PARDISO Memory consumption is %d KBytes\n", iparm[15]+iparm[16]);
	}
	else
	{
	    printf(">>> PARDISO >>> Total PARDISO Memory consumption is %d KBytes\n", iparm[14] );
	}
        for (i = 0 ; i < 64 ; i++)
        {
        	iparmD [i] = iparm[i];
        }


	if (A_values_sym!= NULL)  
	{
	    mxFree (A_values_sym);
	    A_values_sym = NULL;
	}
	if (A_ia_sym!= NULL)  
	{
	    mxFree (A_ia_sym);
	    A_ia_sym = NULL;
	}
	if (A_ja_sym!= NULL)  
	{
	    mxFree (A_ja_sym);
	    A_ja_sym = NULL;
	}

    }
  
     mxFree (iparm) ;
   
}


/* ========================================================================== */
/* === pardisinit_help ====================================================== */
/* ========================================================================== */

static void pardiso_help (void)
{
    mexPrintf ("All information on PARDISO can be obtained at\n") ;
    mexPrintf ("http://www.pardiso-project.org\n") ;
}
