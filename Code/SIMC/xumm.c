/*********************************************************************
 * xumm.c
 * This file shows how to get partial entries of the multiplication of two matrices.
 * For more information, contact me at xum@lamda.nju.edu.cn
 * by: Miao Xu
 ********************************************************************/
#include <matrix.h>
#include <mex.h>   
/* Definitions to keep compatibility with earlier versions of ML */
#ifndef MWSIZE_MAX
typedef int mwSize;
typedef int mwIndex;
typedef int mwSignedIndex;

#if (defined(_LP64) || defined(_WIN64)) && !defined(MX_COMPAT_32)
/* Currently 2^48 based on hardware limitations */
# define MWSIZE_MAX    281474976710655UL
# define MWINDEX_MAX   281474976710655UL
# define MWSINDEX_MAX  281474976710655L
# define MWSINDEX_MIN -281474976710655L
#else
# define MWSIZE_MAX    2147483647UL
# define MWINDEX_MAX   2147483647UL
# define MWSINDEX_MAX  2147483647L
# define MWSINDEX_MIN -2147483647L
#endif
#define MWSIZE_MIN    0UL
#define MWINDEX_MIN   0UL
#endif
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*
     *About Input£∫
     *A:n*r
     *B:r*m
     *row£∫1*|Omega|
     *column£∫1*|Omega|
     *output:A*B' entries decided by row and column: 1*|Omega|
     */
    mxArray *a_in_m, *b_in_m, *c_in_m,*d_in_m, *e_out_m;// ‰»Î «a,b,c,d, ‰≥ˆe
    const mwSize *dims,*dimb,*dimc;
    double *a, *b, *c,*d,*e;
    int dimx, dimy;//, numdims;
    int dimbx,dimby;//,numdimbs;
    int dimcx,dimcy,numdimcs;
    int i,j,k,l;
    int sum;
    a_in_m = mxDuplicateArray(prhs[0]);
    b_in_m = mxDuplicateArray(prhs[1]);
    c_in_m = mxDuplicateArray(prhs[2]);
    d_in_m = mxDuplicateArray(prhs[3]);
    dims = mxGetDimensions(prhs[0]);
	dimy = (int)dims[0]; 
	dimx = (int)dims[1];
    dimb = mxGetDimensions(prhs[1]);
	dimby = (int)dimb[0];
	dimbx = (int)dimb[1];
    dimc = mxGetDimensions(prhs[2]);
	dimcy = (int)dimc[0];
	dimcx = (int)dimc[1];
    //numdims = mxGetNumberOfDimensions(prhs[0]);
    //numdimbs = mxGetNumberOfDimensions(prhs[1]);    
    e_out_m = plhs[0] = mxCreateDoubleMatrix(dimcy,dimcx,mxREAL);
    a = mxGetPr(a_in_m);
    b = mxGetPr(b_in_m);
    c = mxGetPr(c_in_m);
    d = mxGetPr(d_in_m);
    e = mxGetPr(e_out_m);
    for(l=0;l<dimcx;l++) 
    {
           i=c[l]-1;
           j=d[l]-1;
            e[l]=0;
            for(k=0;k<dimx;k++)
            {   
                e[l]+=a[k*dimy+i]*b[j*dimby+k];
            }        
    }
    return;
}
