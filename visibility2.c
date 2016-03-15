#include <stdint.h>
#include "mex.h"
// #include "matrix.h"

void mexFunction(int nl, mxArray *pl[], int nr, const mxArray *pr[])
{
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

    double *K = mxGetPr(mK)
	 , *T = mxGetPr(mT)
	 , *Kin = mxGetPr(mKin)
	 , *Kout = mxGetPr(mKout)
	 ;

    /* Lacasa's algorithm */
    for (uint32_t i=0, g=0; i<(N - 1); i++)
    {
        Kout[i] = Kout[i] + 1;
        Kin[i+1] = Kin[i+1] + 1;
        K[i] = K[i] + 1;
        K[i+1] = K[i+1] + 1;
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
                T[g] = i;
                T[g+N2] = i+j;
                g++;
                criterio = pendiente;
            }
        }
    }

    /* clean up */
    mxDestroyArray(mK);
    mxDestroyArray(mKin);
    mxDestroyArray(mKout);

    /* return just the T matrix */
    pl[0] = mT;
}

/* vim: sw=4 sts=4 et ai
 */
