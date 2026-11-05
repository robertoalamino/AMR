# Reversal Results Across Treatment Durations

The historical and corrected implementations were compared at treatment durations of 500, 1000, 2500 and 5000 steps. Each run used seed 1, 1000 realisations, treatment onset at step 5000 and a final observation at step 15000.

| Treatment steps | Zero final population, historical → corrected | Final $q ≥ 0.27$, historical → corrected | 50% acquisition / recovery, historical | 50% acquisition / recovery, corrected |
|---:|---:|---:|---:|---:|
| 500 | 15.0% → 1.2% | 69.4% → 80.4% | 68 / 1,135 | 82 / 641 |
| 1,000 | 15.1% → 2.3% | 70.6% → 77.8% | 91 / 1,205 | 91 / 735 |
| 2,500 | 17.6% → 1.8% | 68.0% → 78.1% | 98 / 1,375 | 100 / 759 |
| 5,000 | 13.6% → 1.3% | 65.9% → 75.7% | 104 / 1,572 | 125 / 822 |

Crossing times are measured from treatment onset for acquisition and from withdrawal for recovery. Each implementation uses its own decline in ensemble mean $q$, including carried-forward values. The final $q$ threshold percentages in the table use all realisations.

## What Held After Corrections

At the 50% threshold, recovery took 13.2–16.7 times as long as acquisition in the historical runs and 6.6–8.1 times as long in the corrected runs. This supports the qualitative conclusion that resistance reversal was slower than acquisition in these simulations. The crossing times describe ensemble means, they do not estimate how long individual populations
took to recover.

## What Changed After Corrections

The fraction ending with zero recorded population fell by 12.3–15.8 percentage points. The fraction with final $q \geq 0.27$ rose by 7.2–11.0 points across all realisations and by 4.8–7.8 points among final recorded survivors. The direction therefore held when arried-forward values were excluded. The survivor groups differ between implementations, so this is a comparison of their observed survivor distributions.

The change at the low-susceptibility threshold was less clear. Across all realisations, the approximate 95% intervals for the difference in the fraction with final $q \leq 0.1$ included zero at every duration.

## Observation Limits

The corrected ensemble reached 90% recovery after 4973, 5721 and 5545 steps for the 500, 1000 and 2500 steps treatments, respectively. The historical ensemble did not reach that threshold in any of the four runs. Neither implementation reached it after the 5000 steps treatment, which had only 5000 steps of follow-up. Because all runs end at step 15000, final susceptibility across durations reflects both treatment duration and time available for recovery.

These comparisons use one seed per setting. Robustness across seeds and uncertainty in the crossing times were not assessed. The results quantify differences between the historical baseline and corrected implementation. They do not establish exact reproduction of the published figures.

Run details: [500 steps](results-reversal-500.md),
[1,000 steps](results-reversal-1000.md),
[2,500 steps](results-reversal-2500.md), and
[5,000 steps](results-reversal-5000.md).

## Within-realisation Susceptibility Recovery

Each corrected realisation surviving through the final step was compared with its own first recorded treatment-step susceptibility.

| Treatment steps | Final survivors | Final q above own onset | Mean final minus onset q |
|---:|---:|---:|---:|
| 500 | 988 | 627 (63.46%) | −0.005947525 |
| 1000 | 977 | 565 (57.83%) | −0.009283146 |
| 2500 | 982 | 578 (58.86%) | −0.010112316 |
| 5000 | 987 | 456 (46.20%) | −0.017463963 |

No survivor measurements were unavailable. One survivor in the 500-step run had equal onset and final values at output precision. All other survivors not counted as above onset were below it.

Overshoot relative to recorded onset occurred at every duration. Nevertheless, the mean within-realisation change was negative in every run. Decreases in some populations outweighed increases in others.

The baseline is measured after the first treatment step’s life-cycle and death events, rather than immediately before treatment. Differences use the six-decimal output precision. These are descriptive,
survivor-conditioned results from one seed per setting. Treatment duration and post-withdrawal observation time also vary together, therefore the differences across durations cannot be attributed to treatment
duration alone.
