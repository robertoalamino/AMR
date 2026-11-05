#include <stdio.h>
#include <stdlib.h>

#include "lattice.h"

static void check(int site, int direction, int size, int expected)
{
    int actual = lattice_neighbour(site, direction, size);

    if(actual != expected)
    {
        fprintf(stderr,
                "FAIL: site=%d direction=%d size=%d: expected %d, got %d\n",
                site, direction, size, expected, actual);
        exit(EXIT_FAILURE);
    }
}

int main(void)
{
    // Interior site in a 3 x 3 lattice.
    check(5, 1, 3, 6);
    check(5, 2, 3, 4);
    check(5, 3, 3, 2);
    check(5, 4, 3, 8);

    // Wrapping across each edge.
    check(6, 1, 3, 4);
    check(4, 2, 3, 6);
    check(2, 3, 3, 8);
    check(8, 4, 3, 2);

    // Corner wrapping.
    check(1, 2, 3, 3);
    check(1, 3, 3, 7);
    check(9, 1, 3, 7);
    check(9, 4, 3, 3);

    // A single-site lattice is its own neighbour.
    for(int direction = 1; direction <= 4; direction++)
        check(1, direction, 1, 1);

    puts("PASS: lattice neighbour tests.");
    return EXIT_SUCCESS;
}