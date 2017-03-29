/* C++ function [d] = bias_C2(ssh1,ssh2,x1,x2,y1,y2,t1,t2,off,dh,DT,DX)

Computes the difference in SSH from tracks between two data sets, possibly
from different satellite missions

INPUT:   ssh1 : column-vector (n1 x 1) containing the ssh1 of data set 1
         ssh2 : column-vector (n2 x 1) containing the ssh1 of data set 2
         x1, x2 : column-vectors (n1 x 1) and (n2 x 1) constaining longitudes
         y1, y2 : column-vectors (n1 x 1) and (n2 x 1) constaining latitudes
         t1, t2 : column-vectors (n1 x 1) and (n2 x 1) constaining datenums
         off : column-vector (n1 x 1) containing the ref_frame_offset (RADS)
         dh : column-vector (n2 x 1) containing the difference in elevation
              for WGS84 minus TOPEX ellipsoid 
         dt : maximum time difference in days between two points (double)
         dx : maximum distance in km between two points (double)

OUTPUT:  d : column-vector (n2 x 1) containing the difference in SSH

Author : Francisco Mir Calafat (francisco.calafat@noc.ac.uk)
 
*/

#include "mex.h"
#include <cmath>
#include <vector>
#include <stdio.h>
#include <cstring>
#include <algorithm>
#include <matrix.h>
using namespace std;

/* inline function to compute distance on Earth */
const double pi = acos(-1.0);
const double deg2rad = pi/180.0;
const double R = 6371000;
inline double distan(double x1,double y1,double x2,double y2){
    return R*acos(sin(y1*deg2rad)*sin(y2*deg2rad)+
            cos(y1*deg2rad)*cos(y2*deg2rad)*cos((x2-x1)*deg2rad));
}

/* construct the gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    
    
    /* Declare variables */
    const double *ssh1, *ssh2, *x1, *x2, *y1, *y2, *t1, *t2, *off, *dh;
    mwIndex n1, n2;
    mwSize dims[2];
    double *d;
    
    /* the number of inputs should be 5 */
    if (nrhs < 12)
	{
		mexErrMsgTxt("Not enough inputs");
	}
    
    /* the inputs ssh1 and ssh2 must be column vectors */
    if(mxGetM(prhs[0]) < mxGetN(prhs[0]) || mxGetM(prhs[1]) < mxGetN(prhs[1]))
    {
        mexErrMsgTxt("ssh1 and ssh2 must be column vectors");   
    }
    
    /* the inputs x1 and x2 must be column vectors */
    if(mxGetM(prhs[2]) < mxGetN(prhs[2]) || mxGetM(prhs[3]) < mxGetN(prhs[3]))
    {
        mexErrMsgTxt("x1 and x2 must be column vectors");   
    }
    
    if(mxGetM(prhs[0]) !=mxGetM(prhs[8]))
    {
        mexErrMsgTxt("ssh1 must have the same dimensions as off");   
    }
    
    if(mxGetM(prhs[1]) != mxGetM(prhs[9]))
    {
        mexErrMsgTxt("ssh2 must have the same dimensions as dh");   
    }
    
    
    
    
    /* create pointers to the input vectors and matrices */
    ssh1 = mxGetPr(prhs[0]);
    ssh2 = mxGetPr(prhs[1]);
    x1 = mxGetPr(prhs[2]);
    x2 = mxGetPr(prhs[3]);
    y1 = mxGetPr(prhs[4]);
    y2 = mxGetPr(prhs[5]);
    t1 = mxGetPr(prhs[6]);
    t2 = mxGetPr(prhs[7]);
    off = mxGetPr(prhs[8]);
    dh = mxGetPr(prhs[9]);
    const double DT = (double) mxGetScalar(prhs[10]);
    const double DX = (double) mxGetScalar(prhs[11]);
    n1 = mxGetM(prhs[0]);
    n2 = mxGetM(prhs[1]);
    dims[0] = n2;
    dims[1] = 1;
    
    
    /* set the output pointer to the output matrix */
    plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    
    /*  create pointer to a copy of the output matrices */
    d = mxGetPr(plhs[0]);
    
    // start computation
    double dt, dx, dssh;
    int ct;
    for(int i=0; i<n2; ++i){
        dssh = 0;
        ct = 0;
        for(int j=0; j<n1; ++j){
            dx = distan(x2[i],y2[i],x1[j],y1[j])/1000.0;
            dt = abs(t2[i]-t1[j]);
            if (dx < DX & dt < DT){
                ++ct;
                dssh += ssh1[j] - ssh2[i] - off[j] - dh[i];
            }
        }
        if (ct != 0)
            d[i] = dssh/ct;
        else
            d[i] = 0;
    }
    
} //close gateway function
        
        
        