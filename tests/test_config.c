#include <stdio.h>

#include "config.h"

int main(int argc, char *argv[])
{
    RunConfig config = {
        .seed = -1,
        .steps = 120,
        .realisations = 2,
        .treatment_start = 40,
        .treatment_end = 60
    };

    parse_config(argc, argv, &config, 1);

    printf("%ld %d %d %d %d\n",
        config.seed, config.steps, config.realisations,
        config.treatment_start, config.treatment_end);
    return 0;
}