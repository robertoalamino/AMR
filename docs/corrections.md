# Corrections and Scientific Scope

This repository accompanies a reproducibility review of the published model. No journal correction has been or will be submitted or published as part of this review. The historical baseline used for comparison includes earlier maintenance changes and differs from the code as originally published.

## Implementation Changes

- Reproduction now selects among empty neighbours rather than choosing among all four neighbours and failing when the chosen site is occupied.
- Culling tests every cell present at the start of a pass once, correcting skipped entries during population-list compaction.
- Parent selection is without replacement even when reproduction fails. This latent issue does not affect the published setting $r = 1$.
- Early-extinction reversal output is aligned to the observation window. Extinction records and individual sweep outcomes make reporting explicit.

## Archived Analysis and Interpretation

An analysis of archived snapshot observations identified histogram truncation above $q = 0.6$, renormalisation after truncation, and a formula that did not implement the stated $q \geq 0.27$ threshold. Direct archived counts give 69.4%, 70.6%, 68.0% and 65.9% for durations 500, 1000, 2500 and 5000, rather than the approximately 73%, 73%, 72% and 70% reported. These are archived-data recalculations, not corrected simulation results. The [archived reversal snapshot extract](../evidence/amr-archived-evidence/README.md) contains the observations and instructions for reproducing the threshold counts.

Archived final snapshot timing remains partly unresolved. A verified 500-step treatment-end snapshot corresponds to the last treated step. An additional final row is not explained by a header offset. Its origin has not been established. Sensitivity checks for three available durations leave the checked threshold and histogram counts unchanged. The equivalent check is unavailable for the 5000-step archived trajectory.

Histogram percentages describe realisations, not individual cells. Current histograms cover the interval $[0, 1]$ and retain all realisations in their denominator.

## Measured Consequences

The seed-1 sweep retained its broad parameter dependence. Across four reversal durations, recovery remained slower than acquisition at the 50% ensemble-change threshold, while survival and recovery times changed substantially. Individual corrections have not been isolated. See the results reports for numerical details and definitions.

Incomplete recovery within a finite observation window does not establish permanent resistance. Treatment duration and follow-up time vary together. Overshoot relative to recorded onset is distinct from exceeding $q = 0.3$. The recorded onset follows the first treatment step and is not a directly
measured pre-treatment baseline. Robustness across seeds remains untested.
