#include "random.h"

// Test-only generator: initial occupation requires a draw > 0.5,
// so returning 0.25 leaves every site empty.
double ran2(long *seed)
{
    (void)seed;
    return 0.25;
}