# Short Antimicrobial Parameter Sweep Reference

Regression reference generated during repository maintenance.

## Settings

- Lattice: 50 × 50
- Genotype dimension: 50
- Time steps: 30
- Realisations per parameter combination: 2
- Parameter grid: 11 × 11, with both parameters spanning -1 to 1
- Initial random seed: -1, set once before the parameter sweep
- Compiler: GCC 16.2.0, MSYS2
- Optimisation: -O2

## Checks Performed

- Output contains 121 rows and three columns per row.
- All reported death probabilities are finite and within [0, 1].
- Two executions of the same executable produced identical outputs.

This short run checks execution and supports regression testing. It does not establish convergence to asymptotic behaviour or reproduction of the published results. Exact agreement across platforms and compiler versions has not been checked.