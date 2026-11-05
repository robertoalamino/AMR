#ifndef AMR_CONFIG_H
#define AMR_CONFIG_H

typedef struct
{
    long seed;
    int steps;
    int realisations;
    int treatment_start;
    int treatment_end;
} RunConfig;

void parse_config(int argc, char *argv[],
                  RunConfig *config, int has_treatment);

#endif