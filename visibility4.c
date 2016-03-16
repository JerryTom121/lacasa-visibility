#include <stdint.h>

#include "mex.h"

struct row { uint32_t el1, el2; } ;

struct t_list
{
    uint32_t n_row, capacity;
    struct row *rows;
};

struct t_list *new_t_list()
{
    struct t_list *tl = mxMalloc(sizeof(struct t_list));
    if (tl == NULL)
        mexErrMsgTxt("alloc t list failed.");
    tl->n_row = 0;
    tl->capacity = 0;
    tl->rows = NULL;
    return tl;
}

void free_t_list(struct t_list *t_list)
{
    if (t_list->rows != NULL)
        mxFree(t_list->rows);
    mxFree(t_list);
}

#define REALLOC_INCREMENT 100000

void append_t_list(struct t_list *tl, uint32_t el1, uint32_t el2)
{
    if (tl->n_row == tl->capacity)
    {
        tl->rows = mxRealloc(tl->rows, sizeof(struct row) * (tl->capacity + REALLOC_INCREMENT));
        if (tl->rows == NULL)
            mexErrMsgTxt("t list realloc failed.");
        tl->capacity += REALLOC_INCREMENT;
    }
    tl->rows[tl->n_row].el1 = el1;
    tl->rows[tl->n_row].el2 = el2;
    tl->n_row++;
}

void mexFunction(int nl, mxArray *pl[], int nr, const mxArray *pr[])
{
    /* check arguments */
    if (nr != 1 || nl != 4)
        mexErrMsgTxt("[IX_visibility_matrix, total-degree, in-degree, out-degree] = visibility(time_series);");
    
    /* allocate memory for matrix */
    uint32_t N = mxGetM(pr[0]);
    double * restrict H = mxGetPr(pr[0]);
    
    mxArray *mK    = mxCreateNumericMatrix(N, 1, mxUINT32_CLASS, mxREAL)
    , *mKin  = mxDuplicateArray(mK)
    , *mKout = mxDuplicateArray(mK)
    ;
    
    uint32_t * restrict K    = mxGetData(mK)
    , * restrict Kin  = mxGetData(mKin)
    , * restrict Kout = mxGetData(mKout)
    ;
    
    uint32_t g=0;
    struct t_list *tl = new_t_list();
    /* Lacasa's algorithm */
    for (uint32_t i=1; i<=(N - 1); i++)
    {
        Kout[i-1] = Kout[i-1] + 1;
        Kin[i] = Kin[i] + 1;
        K[i-1] = K[i-1] + 1;
        K[i] = K[i] + 1;
        append_t_list(tl, i, i + 1);
        g++;
        double criterio = H[i] - H[i-1];
        for (uint32_t j=2; j<=(N - i); j++)
        {
            double pendiente = (H[i+j-1] - H[i-1]) / j;
            if (criterio < pendiente)
            {
                Kout[i-1] = Kout[i-1] + 1;
                Kin[i+j-1] = Kin[i+j-1] + 1;
                K[i-1] = K[i-1] + 1;
                K[i+j-1] = K[i+j-1] + 1;
                append_t_list(tl, i, i + j);
                g++;
                criterio = pendiente;
            }
        }
    }

    /* convert t list to 2 col index matrix */
    pl[0] = mxCreateNumericMatrix(tl->n_row, 2, mxUINT32_CLASS, mxREAL);
    uint32_t *Tc1 = mxGetData(pl[0]);
    uint32_t *Tc2 = Tc1 + tl->n_row;

    for (uint32_t i=0; i<tl->n_row; i++)
    {
        struct row *row = tl->rows + i;
        Tc1[i] = row->el1;
        Tc2[i] = row->el2;
    }
    
    /* clean up */
    free_t_list(tl);

    /* return other matrices */
    pl[1] = mK;
    pl[2] = mKin;
    pl[3] = mKout;
}

/* vim: sw=4 sts=4 et ai
 */
