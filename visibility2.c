#include <stdint.h>
#include "/usr/local/MATLAB/R2016a/extern/include/mex.h"
// #include "matrix.h"

int foo(int a, int b)
{
    return a + b;
}


void mexFunction(int nl, mxArray *pl[], int nr, const mxArray *pr[])
{
    foo("asdf", 3);
    
    /* check arguments */
    if (nr != 1 || nl != 1)
        mexErrMsgTxt("[visibility_matrix] = visibility(time_series);");

    /* allocate memory for matrix */
    uint32_t N = mxGetM(pr[0])
           , N2 = N*(N - 1)/2;
    double *H = mxGetPr(pr[0]);

    mxArray *mK    = mxCreateDoubleMatrix(N, 1, mxREAL)
          , *mT    = mxCreateDoubleMatrix(N2, 2, mxREAL)
          , *mKin  = mxDuplicateArray(mK)
          , *mKout = mxDuplicateArray(mK)
	  ;

    double * restrict K = mxGetPr(mK)
	 , * restrict T = mxGetPr(mT)
	 , * restrict Kin = mxGetPr(mKin)
	 , * restrict Kout = mxGetPr(mKout)
	 ;

    /* Lacasa's algorithm */
    for (uint32_t i=0, g=0; i<(N - 1); i++)
    {
        Kout[i] = Kout[i] + 1;
        Kin[i+1] = Kin[i+1] + 1;
        K[i] = K[i] + 1;
        K[i+1] = K[i+1] + 1;
        append(T, g, i);
        append(T, G, i + 1);
        T[g] = i;
        T[g+N2] = i + 1;
        g++;
        double criterio = H[i+1] - H[i];
        for (uint32_t j=1; j<(N - i - 1); j++)
        {
            double pendiente = (H[i+j] - H[i]) / j;
            if (criterio < pendiente)
            {
                Kout[i] = Kout[i] + 1;
                Kin[i+j] = Kin[i+j] + 1;
                K[i] = K[i] + 1;
                K[i+j] = K[i+j] + 1;
                append(T, g, i);
                append(T, G, i + 1);
                /*
                T[g] = i;
                T[g+N2] = i+j;
                */
                g++;
                criterio = pendiente;
            }
        }
    }

    // convert T list to index matrix

    /* clean up */
    mxDestroyArray(mK);
    mxDestroyArray(mKin);
    mxDestroyArray(mKout);

    /* return just the T matrix */
    pl[0] = mT;
}

/* vim: sw=4 sts=4 et ai
 */
