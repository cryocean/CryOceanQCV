/* C++ function [Te] = boxes(lon,lat,X,Y,depth,T,levs,DZ,nZ)

INPUT:   lon : vector containing the longitudes of the profiles
         lat : vector containing the latitudes of the profiles
         X : vector containing the longitudes of the grid
         Y : vector containing the latitudes of the grid
         depth : 2D matrix containing the depths of the profiles
         T : 2D matrix containing the temperature of the profiles
         levs : vector containing the standard levels
         DZ : vector containing factors for level decision
         nZ : number of levels

OUTPUT:  Te : 2D matrix containing the temperatures (level x lonlat)
 
 NOTE: this C++ function is called by the MATLAB script "read_profiles2.m"
*/


#include "mex.h"
#include "matrix.h"
#include <cmath>
#include <stdio.h>
#include <vector>
#include <math.h>
using namespace std;

/* >>> construct the gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]) {
    
    
/* Declare variables */
    const int ngrid = 1426; // 5 grados = 2664, 10 grados = 684; 1 grado = 6026
    const double Dx = 1; // 5 grados = 5.0, 10 grados = 10;
    const double Dy = 1; // 5 grados = 3.0, 10 grados = 5.0;
    const int nlev = 38;
    const int *nx = mxGetDimensions(prhs[0]);
    double *Te,*lon,*lat,*X,*Y,*depth,*T,*levs,*DZ;
    const mwSize ndims[] = {nlev,ngrid};
            
/* create pointers to the input vectors and matrices */
    lon = mxGetPr(prhs[0]);
    lat = mxGetPr(prhs[1]);
    X = mxGetPr(prhs[2]);
    Y = mxGetPr(prhs[3]);
    depth = mxGetPr(prhs[4]);
    T = mxGetPr(prhs[5]);
    levs = mxGetPr(prhs[6]);
    DZ = mxGetPr(prhs[7]);
    const double* nZ = (double*) mxGetData(prhs[8]);
    int ny = (int) *nZ;
    
/* set the output pointer to the output matrix */
  plhs[0] = mxCreateNumericArray(2,ndims,mxDOUBLE_CLASS,mxREAL);

/*  create pointer to a copy of the output matrices */
  Te = mxGetPr(plhs[0]);


/* >>> starting computation */
  for(int i=0; i<ngrid; ++i){
      
//----------- find profiles within the corresponding box ------------      
     vector<int> kd;
     for(int ix=0; ix<*nx; ++ix){
         if(*(lon+ix) > *(X+i)-Dx && *(lon+ix) <= *(X+i)+Dx &&
         *(lat+ix) > *(Y+i)-Dy && *(lat+ix) <= *(Y+i)+Dy){
             kd.push_back(ix);}
             }
//------------------------------------------------------ end finding
     
     int skd = kd.size();
     if(skd != 0){ // check that there is at least one profile
    // if so select level
         for(int nz=0; nz<nlev; ++nz){
             double xx = 0;
             int ct = 0;
             // iterate over all the profiles (maybe there is more than one)
                 for(int j=0; j<skd; ++j){ 
                                  
             //----------- find the corresponding level -----------    
                 vector<int> kz;
                 for(int Ix=0; Ix<ny; ++Ix){
                     if(*(depth+kd[j]*ny+Ix)>levs[nz]-DZ[nz] &&
                     *(depth+kd[j]*ny+Ix)<=levs[nz]+DZ[nz+1] &&
                     *(T+kd[j]*ny+Ix) < 99990)
                     kz.push_back(Ix);
                 }
               int skz = kz.size();        
            //---------------------------------------------- end finding 
                 if(skz != 0){
                     for(int jj=0; jj<skz; ++jj){
                     xx = xx+*(T+kd[j]*ny+kz[jj]);
                     ct++;}
                 }
                 }
             *(Te+i*nlev+nz) = xx/ct;
         }
     } // end if
  } // end first for
  
} // close gateway function
