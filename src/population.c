#include "population.h"
#include "random.h"

int population_cull(int *alive, int *occupied, int count,
                    double removal_probability, long *seed)
{
    int survivors = 0;

    for(int i = 1; i <= count; i++)
    {
        int site = occupied[i];

        if(ran2(seed) < removal_probability)
        {
            alive[site] = 0;
        }
        else
        {
            occupied[++survivors] = site;
        }
    }

    for(int i = survivors + 1; i <= count; i++)
        occupied[i] = 0;

    return survivors;
}