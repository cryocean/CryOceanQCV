/*
 function [in,bnd] = inPolygonMex(X,Xv)

 Points inside polygonal region  
 
 USAGE : [in , bnd] = inPolygonMex(X,Xv)
 
 INPUT :      X : points to be tested test as an 2 x N array [x1;y1 x2;y2 etc].
              Xv : vertices of the polygons as an 2 x NXv array [xv1;yv1 xv2;yv2 etc]. The vertices must be
                   specified in consecutive order. If testing for more than 1 polygon, polygons must be
                   separated by NaN ([xv1;yv1 xv2;yv2 ... xvn;yvn NaN;NaN xvn+2;yvn+2 xvn+3;yvn+3... NaN;NaN]).
                   The last polygon must also end with NaN;NaN. 
 
 
 OUTPUT :     in : logical 1 x N vector whose elements are assigned a value of 1 if the associated X point
                   is inside the polygonal region or on the polygon boundary and 0 otherwise.
              bnd : logical 1 x N vector whose elements are assigned a value of 1 if the associated X point
                    is on the polygon boundary and 0 otherwise.
 
 
 To build the MEX-file: 
 mex  inPolygonMex.c

 Mex file based on the inpoly function of Darren Engwirda. Originally created by Sebastien Paris, modified
 by Francisco Mir Calafat (corrected a segmentation violation due to dynamic memory allocation/deallocation
 using standard C malloc and free functions)
 
 */

#include <math.h>
#include "mex.h"

#ifndef max
    #define max(a,b) (a >= b ? a : b)
    #define min(a,b) (a <= b ? a : b)
#endif

/*-------------------------------------------------------------------------------------------------------------- */
/* Function prototypes */

void qsindex( double * , int * , int , int  );

/*-------------------------------------------------------------------------------------------------------------- */

/* construct gateway function */
void mexFunction(int nlhs, mxArray* plhs[] ,  int nrhs, const mxArray* prhs[])
{
	/* variables declaration */
    double *X , *Xv, *edge;
	double *Xtemp , *Xvtemp;
    bool *in , *bnd;
    bool *cn , *on;
    double *x, *y;
    int *index;
    double tol, tol1 , eps = 2.2204e-016, temp , minix=10e20, maxix=-10e20;
    double miniy=10e20, maxiy=-10e20, norminf=-10e20;
    double x1, x2 , y1, y2, xmin, xmax , ymin, ymax, XX, YY , x2x1 , y2y1;
    double lim , ub;
    int i , j, l , i2 , ind , k ;
    int n1 , n2 , lower , upper, start;
    int N, N1, NXv, Nedge;

    /* the number of inputs should be 2 */
    if (nrhs < 2)
	{
		mexPrintf(
			"\n"
			"Not enough inputs.\n"
			);
		return;
	}
    
	/* create pointers to input matrices */  
    X = mxGetPr(prhs[0]);  /* (2 x N) vector */
	Xv = mxGetPr(prhs[1]); /* (2 x NXv) vector */
    N = mxGetN(prhs[0]);
	N1 = N - 1;
	NXv = mxGetN(prhs[1]);
	
    if(mxGetM(prhs[0]) != 2)
    {
        mxErrMsgTxt("X must be (2 x N)");   
    }
	
    if(mxGetM(prhs[1]) != 2)
    {
        mxErrMsgTxt("Xv must be (2 x NXv)");   
    }
    
	Xtemp = mxMalloc(2*N*sizeof(double));
	for(i = 0 ; i < 2*N ; i++)
	{
		Xtemp[i] = X[i];
	}

	Xvtemp = mxMalloc(2*NXv*sizeof(double));
	for(i = 0 ; i < 2*NXv ; i++)
	{
		Xvtemp[i] = Xv[i];
	}
	
    /* build edge */
    Nedge = NXv;
	edge = mxMalloc(2*Nedge*sizeof(double));
    for (i = 0 ; i < Nedge - 1 ; i++)
    {
        i2 = 2*i;
        edge[i2] = i + 1;
        edge[i2 + 1] = i + 2;
        }
     edge[NXv*2-2] = NXv;
     edge[NXv*2-1] = 1.0;

    /* set output pointers to output matrices */
    plhs[0]   = mxCreateLogicalMatrix(1 , N);
	plhs[1]   = mxCreateLogicalMatrix(1 , N);
	
    /* create pointer to a copy of the output matrices */
    in        = (bool *)mxGetPr(plhs[0]);
    bnd       = (bool *)mxGetPr(plhs[1]);
      
    /* Temporal vectors */
    cn = mxMalloc(N*sizeof(bool));
    on = mxMalloc(N*sizeof(bool));
    x = mxMalloc(N*sizeof(double));
    y = mxMalloc(N*sizeof(double));
    index = mxMalloc(N*sizeof(int));
    
    for (i = 0; i < 2*N ; i=i+2)
    {        
        if(X[i] < minix)
        {
            minix = X[i];
        }
        if(X[i] > maxix)
        {
            maxix = X[i];
        }

        if(X[i+1] < miniy)
        {
            miniy = X[i+1];
        }
        if(X[i+1] > maxiy)
        {
            maxiy = X[i+1];
        }
    }
    for (i = 0 ; i < 2*NXv ; i = i + 2)
    {
        if(fabs(Xv[i]) > norminf)
        {
            norminf = fabs(Xv[i]);
        }
        if(fabs(Xv[i + 1]) > norminf)
        {
            norminf = fabs(Xv[i + 1]);
        }
    }
    tol  = norminf*eps;
    lim  = norminf + tol;
	tol1 = tol + 1.0;
    
    if ((maxix - minix) > (maxiy - miniy))
    {
        for (i = 0; i < 2*N ; i=i+2)
        {
            temp             = Xtemp[i];
            Xtemp[i]         = Xtemp[i + 1];
            Xtemp[i + 1]     = temp;
        }
        for (i = 0 ; i < 2*NXv ; i=i+2)
        {
            temp             = Xvtemp[i];
            Xvtemp[i]      = Xvtemp[i + 1];
            Xvtemp[i + 1]  = temp;
        }
    }
    for (i = 0 ; i< N ; i++)
    {
        cn[i]    = 0;
        on[i]    = 0;
        y[i]     = Xtemp[2*i + 1];
        index[i] = i;
    }
    
    qsindex( y , index , 0 , N1 );
    for (i = 0 ; i < N ; i++)
    {
        x[i]     = Xtemp[2*index[i]];
    }
    
    for (k = 0 ; k < Nedge ; k++)
    {
        n1   = 2*(((int) edge[2*k]) - 1);
        n2   = 2*(((int) edge[2*k + 1]) - 1);
        x1   = Xvtemp[n1];
        y1   = Xvtemp[n1 + 1];
        x2   = Xvtemp[n2];
        y2   = Xvtemp[n2 + 1];
		x2x1 = x2 - x1;
		y2y1 = y2 - y1;

        if (x1 > x2)
        {
            xmin = x2;
            xmax = x1;
        }
        else
        {
            xmin = x1;
            xmax = x2;
        }
        if (y1 > y2)
        {
            ymin = y2;
            ymax = y1;
        }
        else
        {
            ymin = y1;
            ymax = y2;
        }
        if (y[0] == ymin)
        {
            start = 0;
        }
        else if (y[N1] <= ymin)
        {
            start = N1;
        }
        else
        {
            lower = 0;
            upper = N1;
            start = ((lower+upper)/2);
            for (l = 0 ; l < N ; l++)
            {
                if (y[start] < ymin)
                {
                    lower = start;
                    start = ((lower+upper)/2);
                }
                else if (y[start-1]<ymin)
                {
                    break;
                }
                else
                {
                    upper = start;
                    start = ((lower+upper)/2);
                }
            }
            start--;
        }

		start = max(0 , start);
        for (j = start ; j < N ; j++)
        {
            YY = y[j];
            if (YY <= ymax)
            {
                if (YY >= ymin)
                {
                    XX = x[j];
                    
                    if (XX >= xmin)
                    {
                        if (XX <= xmax)
                        {
                            on[j] = on[j] || ( fabs( (y2 - YY)*(x1 - XX) - (y1 - YY)*(x2 - XX) ) < tol );
                            if (YY < ymax)
                            {
                                ub = ( x2x1*(y1 - YY) - y2y1*(x1 - XX) )/( (XX - lim)*y2y1 );
                                if ( (ub > -tol) && (ub < tol1 ) )
                                {
                                    cn[j] = !cn[j];
                                }
                            }
                        }
                    }
                    else if (YY < ymax)   
                    {
                        cn[j] = !cn[j];   
                    }    
                }
            }
            else   
            {
                break;
            }   
        }    
    }
    for(i = 0 ; i < N ; i++)
    {
        ind      = index[i];
        in[ind]  = (cn[i] || on[i]);
        bnd[ind] = on[i];   
    }
    /* deallocate dynamic memory */
    mxFree(cn);
    mxFree(on);
    mxFree(x);
    mxFree(y);
    mxFree(index);
    mxFree(Xtemp);
    mxFree(Xvtemp);
    mxFree(edge);
}

/* ----------------------------------------------------------------------------- */
void qsindex (double  *a, int *index , int lo, int hi)
{
/*  
  lo is the lower index, hi is the upper index
  of the region of array a that is to be sorted
*/
    int i=lo, j=hi , ind;
    double x=a[(lo+hi)/2] , h;

    /*  partition */
    do
    {    
        while (a[i]<x) i++; 
        while (a[j]>x) j--;
        if (i<=j)
        {
            h        = a[i]; 
			a[i]     = a[j]; 
			a[j]     = h;
			ind      = index[i];
			index[i] = index[j];
			index[j] = ind;
            i++; 
			j--;
        }
    }
	while (i<=j);

    /*  recursion */
    if (lo<j) qsindex(a , index , lo , j);
    if (i<hi) qsindex(a , index , i , hi);
}
/* ----------------------------------------------------------------------------- */
