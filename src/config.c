#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"

static void usage(const char *program, int has_treatment)
{
    fprintf(stderr,
            "Usage: %s [--seed N] [--steps N] [--realisations N]",
            program);

    if(has_treatment)
        fprintf(stderr, " [--start N] [--end N]");

    fprintf(stderr, "\n");
}

static long positive_integer(const char *option, const char *text,
                             long maximum)
{
    char *end;
    long value;

    errno = 0;
    value = strtol(text, &end, 10);

    if(errno == ERANGE || end == text || *end != '\0' ||
       value <= 0 || value > maximum)
    {
        fprintf(stderr, "%s requires an integer between 1 and %ld.\n",
                option, maximum);
        exit(EXIT_FAILURE);
    }

    return value;
}

void parse_config(int argc, char *argv[],
                  RunConfig *config, int has_treatment)
{
    for(int i = 1; i < argc; i += 2)
    {
        const char *option = argv[i];

        if(i + 1 >= argc)
        {
            fprintf(stderr, "Missing value for %s.\n", option);
            usage(argv[0], has_treatment);
            exit(EXIT_FAILURE);
        }

        const char *value = argv[i + 1];

        if(strcmp(option, "--seed") == 0)
            config->seed = -positive_integer(option, value, LONG_MAX);

        else if(strcmp(option, "--steps") == 0)
            config->steps =
                (int)positive_integer(option, value, INT_MAX - 1);

        else if(strcmp(option, "--realisations") == 0)
            config->realisations =
                (int)positive_integer(option, value, INT_MAX);

        else if(has_treatment && strcmp(option, "--start") == 0)
            config->treatment_start =
                (int)positive_integer(option, value, INT_MAX);

        else if(has_treatment && strcmp(option, "--end") == 0)
            config->treatment_end =
                (int)positive_integer(option, value, INT_MAX);

        else
        {
            fprintf(stderr, "Unknown option: %s.\n", option);
            usage(argv[0], has_treatment);
            exit(EXIT_FAILURE);
        }
    }

    if(has_treatment &&
       (config->treatment_start >= config->treatment_end ||
        config->treatment_start > config->steps ||
        config->treatment_end > config->steps + 1))
    {
        fprintf(stderr,
                "Treatment must satisfy "
                "1 <= start < end <= steps + 1, with start <= steps.\n");
        exit(EXIT_FAILURE);
    }
}