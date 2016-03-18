/* Lacasa's visibility graph algorithm.
 *
 * Compile with 
 *
 * >>> mex CFLAGS='-fPIC -std=c99 -O3 -ffast-math -march=native' visibility4.c
 *
 */

#include <stdint.h>

#include "mex.h"

struct ll { uint32_t count, capacity, *data; } ll_zero = {0};

void ll_append(struct ll *ll, uint32_t el1, uint32_t el2)
{
    if (ll->count == ll->capacity)
    {
        ll->data = mxRealloc(ll->data, sizeof(uint32_t) * 2 * (ll->capacity + 1000000));
        if (ll->data == NULL)
            mexErrMsgTxt("t list realloc failed.");
        ll->capacity += 1000000;
    }
    ll->data[ll->count*2 + 0] = el1;
    ll->data[ll->count*2 + 1] = el2;
    ll->count++;
}

struct lacasa_data
{
    uint32_t M, g;
    struct ll ll;
    mxArray *K, *Kin, *Kout, *VG;
    const mxArray *H;
};

void lacasa_compute(struct lacasa_data *data)
{
    data->M = mxGetNumberOfElements(data->H);
    double * restrict H = mxGetPr(data->H);
   
    data->g = 0;
    data->K = mxCreateNumericMatrix(data->M, 1, mxUINT32_CLASS, mxREAL);
    data->Kin = mxDuplicateArray(data->K);
    data->Kout = mxDuplicateArray(data->K);
    
    uint32_t * restrict K = mxGetData(data->K)
           , * restrict Kin = mxGetData(data->Kin)
           , * restrict Kout = mxGetData(data->Kout)
           ;
    
    for (uint32_t i=1; i<=(data->M - 1); i++)
    {
        Kout[i-1] = Kout[i-1] + 1;
        Kin[i] = Kin[i] + 1;
        K[i-1] = K[i-1] + 1;
        K[i] = K[i] + 1;
        ll_append(&data->ll, i, i + 1);
        data->g++;
        double criterio = H[i] - H[i-1];
        for (uint32_t j=2; j<=(data->M - i); j++)
        {
            double pendiente = (H[i+j-1] - H[i-1]) / j;
            if (criterio < pendiente)
            {
                Kout[i-1] = Kout[i-1] + 1;
                Kin[i+j-1] = Kin[i+j-1] + 1;
                K[i-1] = K[i-1] + 1;
                K[i+j-1] = K[i+j-1] + 1;
                ll_append(&data->ll, i, i + j);
                data->g++;
                criterio = pendiente;
            }
        }
    }
    
    data->VG = mxCreateNumericMatrix(0, 0, mxUINT32_CLASS, mxREAL);
    mxSetData(data->VG, data->ll.data);
    mxSetM(data->VG, 2);
    mxSetN(data->VG, data->ll.count);
}

void mexFunction(int nl, mxArray *pl[], int nr, const mxArray *pr[])
{
    /* check arguments */
    if (nr != 1 || nl != 4)
        mexErrMsgTxt("[IX_visibility_matrix, total-degree, in-degree, out-degree] = visibility(time_series);");
    
    /* setup & compute v g */
    struct lacasa_data ld = {0};
    ld.H = pr[0];    
    lacasa_compute(&ld);

    /* return results */
    pl[0] = ld.VG;
    pl[1] = ld.K;
    pl[2] = ld.Kin;
    pl[3] = ld.Kout;
}

/* vim: sw=4 sts=4 et ai
 */
