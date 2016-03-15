#include <stdint.h>
#include "mex.h"
#include "matrix.h"

void mexFunction(int nl, mxArray *pl[], int nr, const mxArray *pr[])
{
    /* check arguments */
    if (nr != 1 || nl != 1)
        mexErrMsgTxt("[visibility_matrix] = visibility(time_series);");
    
    /* allocate memory for matrix */
    uint32_t len = mxGetM(pr[0]);
    double *y = mxGetPr(pr[0]);
    pl[0] = mxCreateNumericMatrix(len, len, mxUINT8_CLASS, mxREAL);
    uint8_t *vis = mxGetData(pl[0]);
    
    /* zero matrix */
    for (uint32_t i=0; i<(len*len); i++)
        vis[i] = 0;
    
    /* define macro for 2D access */
    #define V(i, j) vis[i + len*j]
    
    /* compute visibility */
    for (uint32_t a=0; a<len; a++)
    {
        double y_a = y[a];
        
        /* neighbors always visible (?) */
        if (a > 0)
            V(a, a - 1) = 1;
        
        if (a < (len - 1))
            V(a, a + 1) = 1;
        
        /* all non-neighbors after i, test */
        for (uint32_t b=a+2; b<len; b++)
        {
            double y_b = y[b];
            V(a, b) = 1;
            for (uint32_t c=a+1; c<b; c++)
                if (y[c] >= (y_b + (y_a - y_b)*((b - c)*1.0)/(b - a)))
                {
                    V(a, b) = 0;
                    break;
                }
        }
    }
}