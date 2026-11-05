# Short reversal simulation reference

These outputs were generated after the initial code corrections, before structural refactoring.

## Settings

- Lattice: 50 × 50
- Genotype dimension: 50
- Time steps: 120
- Treatment start: 40
- Treatment end: 60 (exclusive)
- Realisations: 2
- Initial random seed: -1, set once before the realisation loop
- Compiler: GCC 16.2.0, MSYS2
- Optimisation: -O2

## Checks performed

- Both output files contain two rows.
- Each population row contains 121 observations.
- Each death-probability row contains 81 observations.
- Two executions of the same executable produced byte-identical outputs.

These files provide a regression reference for subsequent refactoring. They do not establish scientific validity or reproduction of the published results. Exact agreement across platforms and compiler versions has not been checked.