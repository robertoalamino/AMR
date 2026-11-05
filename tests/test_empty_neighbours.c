#include <stdio.h>
#include <stdlib.h>

#include "lattice.h"

static void check(const char *name, int site, const int *occupied,
                  int expected_count, const int *expected)
{
    int empty[4];
    int count = lattice_empty_neighbours(site, 3, occupied, empty);

    if(count != expected_count)
    {
        fprintf(stderr, "FAIL: %s: expected %d neighbours, got %d.\n",
                name, expected_count, count);
        exit(EXIT_FAILURE);
    }

    for(int i = 0; i < count; i++)
    {
        if(empty[i] != expected[i])
        {
            fprintf(stderr,
                    "FAIL: %s: entry %d: expected %d, got %d.\n",
                    name, i, expected[i], empty[i]);
            exit(EXIT_FAILURE);
        }
    }
}

int main(void)
{
    // One-based indexing; initially all nine sites are occupied.
    int occupied[10];
    for(int i = 0; i < 10; i++)
        occupied[i] = 1;

    check("no empty neighbours", 5, occupied, 0, NULL);

    occupied[6] = 0;
    int one[] = {6};
    check("one empty neighbour", 5, occupied, 1, one);

    occupied[2] = 0;
    int two[] = {6, 2};
    check("two empty neighbours", 5, occupied, 2, two);

    occupied[4] = 0;
    occupied[8] = 0;
    int four[] = {6, 4, 2, 8};
    check("four empty neighbours", 5, occupied, 4, four);

    // Corner site 1: left wraps to 3, up wraps to 7.
    for(int i = 0; i < 10; i++)
        occupied[i] = 1;

    occupied[3] = 0;
    occupied[7] = 0;
    int corner[] = {3, 7};
    check("corner wrapping", 1, occupied, 2, corner);

    puts("PASS: empty neighbour enumeration.");
    return EXIT_SUCCESS;
}