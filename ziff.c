#include "stdlib.h"
#include "stdio.h"
#include "math.h"

#define M  16383
#define NewRandomInteger (++nd,ra[nd&M] = ra[(nd-471)&M]^ra[(nd-1586)&M]^ra[(nd-6988)&M]^ra[(nd-9689)&M])

long     ra[M+1], nd, j;

double rann_()
{
  return (NewRandomInteger/2147483648.0);
}

void randinit_(long seed)
{
    double a, ee = -1 + 1/2147483648.0;
    long i;
    extern long nd, ra[M+1];

    a = seed/2147483648.0;
    for (nd = 0; nd <= M; nd++)
    {
        a *= 16807;
        a += ee * (long)(a);
        if (a >= 1) a += ee;
        ra[nd] = (long) (2147483648.0 * a);
    }
    nd = M;
}
