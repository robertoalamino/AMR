#include "genotype.h"

double total(int *v, int n)
{
    int i;
    double r;

    r = 0;
    for(i=1;i<=n;i++) r += v[i];

    return r;
}