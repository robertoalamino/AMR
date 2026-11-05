#ifndef AMR_POPULATION_H
#define AMR_POPULATION_H

// Test each listed cell once for removal.
// Compact survivors into occupied[1..new_count] and return new_count.
// removal_probability must be between 0 and 1.
int population_cull(int *alive, int *occupied, int count,
                    double removal_probability, long *seed);

#endif