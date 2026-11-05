#include <stdio.h>
#include <stdlib.h>

#include "lattice.h"

int lattice_neighbour(int site, int direction, int size)
{
    int sites = size * size;

    switch(direction)
    {
        case 1:
            if(site % size == 0) return site - size + 1;
            return site + 1;

        case 2:
            if((site - 1) % size == 0) return site + size - 1;
            return site - 1;

        case 3:
            if(site - size < 1) return site - size + sites;
            return site - size;

        case 4:
            if(site + size > sites) return site + size - sites;
            return site + size;

        default:
            fprintf(stderr, "Invalid neighbour direction: %d\n",
                    direction);
            exit(EXIT_FAILURE);
    }
}

int lattice_empty_neighbours(int site, int size, const int *occupied, int empty[4])
{
    int count = 0;

    for(int direction = 1; direction <= 4; direction++)
    {
        int neighbour = lattice_neighbour(site, direction, size);

        if(occupied[neighbour] == 0)
            empty[count++] = neighbour;
    }

    return count;
}