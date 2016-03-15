#include <stdint.h>
#include "mex.h"
#include "matrix.h"

void mexFunction(int nl, mxArray *pl[], int nr, 
  const mxArray *pr[])
{
    /* check arguments */
    if (nr != 1 || nl != 1)
        mexErrMsgTxt("[series] = conway(length);");
    
    /* allocate memory for array */
    uint64_t len = mxGetScalar(pr[0]);
    mxArray *series = mxCreateNumericMatrix(len + 1, 1, mxUINT64_CLASS, mxREAL);
    uint64_t *a = mxGetData(series);
    
    /* compute Conway series */
    a[1] = a[2] = 1;
    for (uint64_t n=3; n<(len + 1); n++)
        a[n] = a[a[n - 1]] + a[n - a[n - 1]];
    
    /* convert to conventional form, a(n) - n/2, in double precision */
    pl[0] = mxCreateDoubleMatrix(len, 1, mxREAL);
    double *a_ = mxGetPr(pl[0]);
    for (uint64_t n=0; n<len; n++)
        a_[n] = (double) a[n + 1]  - n / 2.0;
    
    /* clean up */
    mxDestroyArray(series);
}