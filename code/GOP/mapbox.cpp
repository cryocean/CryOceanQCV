/* C++ function [h] = mapbox2(lon,lat,z,X,Y,dx)

INPUT:   lon : vector containing the longitudes of the sla
         lat : vector containing the latitudes of the sla
         z : vector containing the field (sla, swh,...)
         X : vector containing the longitudes of the grid
         Y : vector containing the latitudes of the grid
         dx : half length

OUTPUT:  h : 1D matrix containing the sla (1 x lonlat)
 
*/

#include "mex.h"
#include "matrix.h"
#include <cmath>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <cstring>
using namespace std;

/* >>> construct the gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    
    
    /* Declare variables */
    const double *X, *Y, *dx;
    mwIndex nw, ngrid, nr;
    double *h;
    mxArray *lon, *lat;
    double *lontmp, *lattmp, *slatmp;
    mwSize dims[2];
    
    /* create pointers to the input vectors and matrices */
    X = mxGetPr(prhs[3]);
    Y = mxGetPr(prhs[4]);
    dx = mxGetPr(prhs[5]);
    nw = mxGetM(prhs[2]);
    ngrid = mxGetM(prhs[3]);
    dims[0] = nw;
    dims[1] = ngrid;
    
    /* set the output pointer to the output matrix */
    plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    
    /*  create pointer to a copy of the output matrices */
    h = mxGetPr(plhs[0]);
    
// starting computation
    int ct[ngrid];
    double tmp[ngrid];
    for(int i=0; i<nw; ++i){
        cout<<i<<endl;
        memset(ct,0,sizeof(ct));
        memset(tmp,0,sizeof(tmp));
        lontmp = mxGetPr(mxGetCell(prhs[0], i));
        lattmp = mxGetPr(mxGetCell(prhs[1], i));
        slatmp = mxGetPr(mxGetCell(prhs[2], i));
        nr = mxGetN(mxGetCell(prhs[0], i));
        for(int k=0; k<nr; ++k){
            for(int j=0; j<ngrid; ++j){
                if (lontmp[k] >= X[j]-*dx && lontmp[k] < X[j]+*dx &&
                        lattmp[k] >= Y[j]-*dx && lattmp[k] < Y[j]+*dx){
                    tmp[j] += slatmp[k];
                    ct[j]++;
                }
            }
        }
        for(int j=0; j<ngrid; ++j){
            if(ct[j] == 0){
                continue;}
            h[i+j*nw] = tmp[j]/ct[j];
        }
    }
    
}
