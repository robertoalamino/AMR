# Experiment settings

All reported comparisons used seed 1 and a continuing random stream.

| Program | Steps | Realisations | Treatment interval | Mutation rate |
|---|---:|---:|---|---:|
| Sweep | 3000 | 100 per pair | Not applicable | 0 |
| Reversal | 15000 | 1000 | [5000, 5500) | 0.001 |
| Reversal | 15000 | 1000 | [5000, 6000) | 0.001 |
| Reversal | 15000 | 1000 | [5000, 7500) | 0.001 |
| Reversal | 15000 | 1000 | [5000, 10000) | 0.001 |

The grid has 121 alpha/beta pairs spanning -1 to 1 in steps of 0.2.
Each setting is run with either `historical` or `corrected` sources.

Short checks use 30 steps and two realisations per pair for the sweep,
and 120 steps, two realisations and [40, 60) treatment for reversal.
These check execution and analysis, not long-time convergence.
