# Archived reversal snapshot evidence

This extract contains 8,000 numeric observations from the archived distribution
worksheet: two snapshots of 1,000 realisations for each of four durations.
The `index` column preserves row order within each snapshot. The source columns
were selected by duration in the combined observations table. No original
workbooks, formulas, workbook metadata, or local paths are included.

These are observations used in the archived analysis, separate from the
historical-baseline and corrected simulation reruns. `archived_final` names
the observation's role in the archived analysis; its exact simulation step
has not been established. For the 5000-step duration the source trajectory
was unavailable for timing verification.

From this directory, check file integrity and reproduce the final counts:

```bash
sha256sum -c checksums.sha256
awk -F '\t' 'NR > 1 && $2 == "archived_final" {
  n[$1]++; ge[$1] += ($4 >= 0.27); above[$1] += ($4 > 0.6)
}
END {
  for (d = 500; d <= 5000; d++) if (d in n)
    printf "%d: %d/%d at or above 0.27 (%.1f%%); %d above 0.6\n",
           d, ge[d], n[d], 100*ge[d]/n[d], above[d]
}' snapshots.tsv
```

| Duration | Final q >= 0.27 | Final q > 0.6 |
|---:|---:|---:|
| 500 | 694/1000 | 2 |
| 1000 | 706/1000 | 2 |
| 2500 | 680/1000 | 1 |
| 5000 | 659/1000 | 1 |

The threshold denominator is all 1,000 realisations, including observations
above the archived histogram's display limit of 0.6.
