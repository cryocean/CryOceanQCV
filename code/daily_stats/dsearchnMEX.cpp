/* C++ function [d] = distance2Coast(x,y,xc,yc,maxd)

Computes the distance from the coast for distances smaller than maxd

INPUT:   x : column-vector (n x 1) containing the longitudes of observations
         y : column-vector (n x 1) containing the latitudes of observations
         xc : column-vector (nc x 1) containing the longitudes of the coastline
         yc : column-vector (nc x 1) containing the latitudes of the coastline
              yc should be sorted in ascending order.
         maxd : distance in km (int). The search is restricted to points closer
                than maxd from the coast. Using smaller values for maxd significantly
                speeds up the computation. maxd must be a multiple of 100.

OUTPUT:  d : column-vector (n x 1) containing the distances to coast for each 
             (x,y) point.

Note: xc and yc must be the same length and should not contain NaN values.

Author : Francisco Mir Calafat (francisco.calafat@noc.ac.uk)
 
*/

#include "mex.h"
#include <cmath>
#include <stdio.h>
#include <cstring>
#include <algorithm>
#include <matrix.h>
using namespace std;




/* construct the gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    
    
    /* Declare variables */
    const double *x, *y, *xc, *yc;
    mwIndex n, nc;
    mwSize dims[2];
    double *d;
    int *k;
    
    /* the number of inputs should be 4 */
    if (nrhs < 4)
	{
		mexErrMsgTxt("Not enough inputs");
	}
    
    /* the inputs x and y must be column vectors */
    if(mxGetM(prhs[0]) < mxGetN(prhs[0]) || mxGetM(prhs[1]) < mxGetN(prhs[1]))
    {
        mexErrMsgTxt("x and y must be column vectors");   
    }
    
    /* the inputs xc and yc must be column vectors */
    if(mxGetM(prhs[2]) < mxGetN(prhs[2]) || mxGetM(prhs[3]) < mxGetN(prhs[3]))
    {
        mexErrMsgTxt("xc and yc must be column vectors");   
    }
    
    
    /* create pointers to the input vectors and matrices */
    x = mxGetPr(prhs[0]);
    y = mxGetPr(prhs[1]);
    xc = mxGetPr(prhs[2]);
    yc = mxGetPr(prhs[3]);
    n = mxGetM(prhs[0]);
    nc = mxGetM(prhs[2]);
    dims[0] = nc;
    dims[1] = 1;
    
    
    /* set the output pointer to the output matrix */
    plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(2, dims, mxUINT32_CLASS, mxREAL);
    
    /*  create pointer to a copy of the output matrices */
    d = mxGetPr(plhs[0]);
    k = (int*) mxGetData(plhs[1]);
    
    double dtmp;
    double dtmp0;
    int ktmp;
    for(int i=0;i<nc;++i){
        dtmp0 = 1e20;
        for(int j=0;j<n;++j){
            if (abs(y[j]-yc[i]) > 5 || abs(x[j]-xc[i]) > 20)
                    continue;
           dtmp = (x[j]-xc[i])*(x[j]-xc[i])+(y[j]-yc[i])*(y[j]-yc[i]);
           if (dtmp < dtmp0){
               dtmp0 = dtmp;
               ktmp = j;}
        }
        d[i] = sqrt(dtmp0);
        k[i] = ktmp+1;
    }

}