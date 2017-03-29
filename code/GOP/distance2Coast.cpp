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


/* structure that will contain the indices and number of points for each area */
struct ind{
    int nk;
    vector<int> k;
};

/* construct the gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    
    
    /* Declare variables */
    const double *x, *y, *xc, *yc;
    mwIndex n, nc;
    mwSize dims[2];
    double *d;
    
    /* the number of inputs should be 5 */
    if (nrhs < 5)
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
    const int maxd = (int) mxGetScalar(prhs[4]);
    n = mxGetM(prhs[0]);
    nc = mxGetM(prhs[2]);
    dims[0] = n;
    dims[1] = 1;
    
    /* the input maxd must be a multiple of 100 */
    if(maxd%100 != 0)
    {
        mexErrMsgTxt("maxd must be a multiple of 100");   
    }
    const int MAXD = maxd/100;
    
    /* set the output pointer to the output matrix */
    plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    
    /*  create pointer to a copy of the output matrices */
    d = mxGetPr(plhs[0]);
    
    /* build areas of 1º x 1º based on the coastline data */
    const int zoney1 = (int)floor(yc[0])+90;
    const int zoney2 = (int)floor(yc[nc-1])+90;
    const int nzonesy = zoney2-zoney1+1;
   // const int nzonesy = 170;
    const int nzonesx = 361;
    const int offset = 90-zoney1;
    const int offsetx = 180;
    int k0[nzonesy];
    int kf[nzonesy];
    memset(k0,0,sizeof(k0));
    memset(kf,-1,sizeof(kf));
    int dz;
    kf[nzonesy-1] = nc-1;
    for(int i=1; i<nc; ++i){
        dz = floor(yc[i]) - floor(yc[i-1]);
        if (dz != 0){
            kf[(int)(floor(yc[i-1])+offset)] = i-1;
            k0[(int)(floor(yc[i])+offset)] = i;
        }
    }
    
    /* store indices for each area in structure */
    ind zoneInd[nzonesy][nzonesx]; // structure containing indices
    int np,zonex;
    for(int i=0; i<nzonesy; ++i){
        np = kf[i]-k0[i]+1;
        for(int j=k0[i]; j<np+k0[i]; ++j){
            zonex = (int)(floor(xc[j])+offsetx);
            zoneInd[i][zonex].k.push_back(j);
        }
        for(int j=0; j<nzonesx; ++j)
            zoneInd[i][j].nk = zoneInd[i][j].k.size();
    }
    
    
    /* start computation */
    double dtmp;
    double tmp;
    int zoney, iy, fy,ix,fx;
    int indy,lati,nx;
    int DX[nzonesy];
    for(int i=0; i<nzonesy; ++i){
        lati = abs(i-offset);
        if(lati <= 60)
            DX[i] = 2*MAXD;
        else if(lati <=70)
            DX[i] = 3*MAXD;
        else if(lati <= 80)
            DX[i] = 6*MAXD;
        else
            DX[i] = 20*MAXD;
    }
    
    int rs;
    for(int i=0; i<n; ++i){
        dtmp = 1e+30;
        zoney = (int)(floor(y[i])+offset);
        if(zoney > -2 || zoney < nzonesy+1){
            iy = max(0,zoney-MAXD);
            fy = min(nzonesy-1,zoney+MAXD);
            iy = min(iy,nzonesy-1);
            fy = max(0,fy);
            zonex = (int)(floor(x[i])+offsetx);
            ix = zonex-DX[max(abs(iy),abs(fy))];
            fx = zonex+DX[max(abs(iy),abs(fy))];
            for(int j=iy; j<fy+1; ++j){
                for(int r=ix; r<fx+1; ++r){
                    rs = r%361;
                    if(rs < 0)
                        rs = r + 361;
                    for(int jj=0; jj<zoneInd[j][rs].nk; ++jj){
                        indy = zoneInd[j][rs].k[jj];
                        tmp = distan(x[i],y[i],xc[indy],yc[indy]);
                        if (tmp < dtmp)
                            dtmp = tmp;
                    }
                }
            }
        }
        d[i] = dtmp;
    }

}
