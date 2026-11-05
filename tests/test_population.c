#include <stdio.h>
#include <stdlib.h>

#include "population.h"
#include "random.h"

static const double draws[] = {0.1, 0.8, 0.2, 0.9};
static int draw_count;

// Test replacement for the production random generator.
double ran2(long *seed)
{
    (void)seed;

    if(draw_count >= 4)
    {
        fprintf(stderr, "FAIL: too many random draws.\n");
        exit(EXIT_FAILURE);
    }

    return draws[draw_count++];
}

static void require(int condition, const char *message)
{
    if(!condition)
    {
        fprintf(stderr, "FAIL: %s\n", message);
        exit(EXIT_FAILURE);
    }
}

static void check_case(double probability, const int expected[4])
{
    // Deliberately non-consecutive occupied sites.
    int alive[8] = {0};
    int occupied[5] = {0, 2, 4, 6, 7};
    const int sites[4] = {2, 4, 6, 7};
    long seed = -1;
    int expected_count = 0;

    for(int i = 0; i < 4; i++)
        alive[sites[i]] = 1;

    draw_count = 0;

    int count = population_cull(alive, occupied, 4,
                                probability, &seed);

    require(draw_count == 4, "each original cell must be tested once");

    for(int i = 0; i < 4; i++)
    {
        require(alive[sites[i]] == expected[i],
                "incorrect lattice occupation");

        if(expected[i])
        {
            expected_count++;
            require(occupied[expected_count] == sites[i],
                    "incorrect survivor list");
        }
    }

    require(count == expected_count, "incorrect survivor count");

    for(int i = count + 1; i <= 4; i++)
        require(occupied[i] == 0, "unused list entries must be cleared");

    require(alive[1] == 0 && alive[3] == 0 && alive[5] == 0,
            "unlisted sites must remain unchanged");
}

int main(void)
{
    const int all_survive[4] = {1, 1, 1, 1};
    const int none_survive[4] = {0, 0, 0, 0};
    const int mixed[4] = {0, 1, 0, 1};

    check_case(0.0, all_survive);
    check_case(1.0, none_survive);
    check_case(0.5, mixed);

    int alive[1] = {0};
    int occupied[1] = {0};
    long seed = -1;

    draw_count = 0;
    require(population_cull(alive, occupied, 0, 0.5, &seed) == 0,
            "empty population must remain empty");
    require(draw_count == 0, "empty population must consume no draws");

    puts("PASS: culling decisions and population bookkeeping.");
    return EXIT_SUCCESS;
}