/* ========================================================================== */
/* === pardisoinit mexFunction ============================================== */
/* ========================================================================== */

/*
    Usage:

    For compilation please look to the README file.
    
    Return the vector iparm and the internal address pointer pt for
    PARDISO solver version 3.3
    
    Example:

    % for solving indefinite symmetric systems 
    mtype = -2 ;
    [iparm, pt] = pardisoinit(mtype);

    % iparm(3) must always be equal to the environment variable OMP_NUM_THREADS
    iparm (3) = OMP_NUM_THREADS;


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

	Copyright (c), 2009, 2008, 2007, 2006, 2005, 2004  University of Basel.
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

static void pardisoinit_help (void) ;

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
    /* === Local variables ================================================== */
 
    /* PardisoInit control parameters. */
    int      *iparm;
    double   *iparm_out;
   
    double      *pt;
    double   *pt_out;

    int       mtype, i;

    /* === Check inputs ===================================================== */

    if (nlhs != 2 || nrhs != 1)
    {
    	pardisoinit_help () ;
    	mexErrMsgTxt (
    	"pardisoinit: incorrect number of input and/or output arguments.") ;
    }
    
    /* === pardiso init ====================================================== */

    pt	  = (double *) mxCalloc (64, sizeof (double)) ;
    iparm = (int *) mxCalloc (64, sizeof (int)) ;

    mtype = (int)(mxGetScalar(prhs[0]));
 
    pardisoinit_(pt,  &mtype, iparm);
    iparm[0] = 1;
    iparm[10] = 1;
    iparm[12] = 2;

    /* === Return the handle pointer  ==================================== */

    plhs [0]  = mxCreateDoubleMatrix (1, 64, mxREAL) ;
    iparm_out = mxGetPr (plhs [0]);

    for (i = 0 ; i < 64 ; i++)
    {
    	iparm_out [i] = iparm [i];
    }

    plhs [1]  = mxCreateDoubleMatrix (1, 64, mxREAL) ;
    pt_out = mxGetPr (plhs [1]);

    for (i = 0 ; i < 64 ; i++)
    {
    	pt_out [i] = pt [i];
    }

    mxFree (pt) ;
    mxFree (iparm) ;
}


/* ========================================================================== */
/* === pardisinit_help ====================================================== */
/* ========================================================================== */

static void pardisoinit_help (void)
{
    mexPrintf ("All information on PARDISO can be obtained at\n") ;
    mexPrintf ("http://www.pardiso-project.org\n") ;
}

