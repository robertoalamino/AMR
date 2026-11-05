# Numerical inputs and recorded results

These are data from the completed seed-1 experiments, not new simulations
of the release candidate. Original source revision identifiers are recorded
in the portable provenance files. They identify research revisions and may
not be reachable in the release repository's clean history.

## Contents

- `sweep/`: historical and corrected aggregate outputs (no header), corrected
  individual outcomes and extinction records (no header), the corrected
  summary, the aggregate comparison and a descriptive difference report.
- `reversal/{500,1000,2500,5000}/`: corrected summary trajectories, historical
  mean trajectories, corrected snapshots and histograms for the last treated
  step, first untreated step and final step, corrected extinction records,
  and the recorded crossing-time and final-comparison reports.
- Per-setting `provenance/`: source revision identifiers, arguments, compiler
  details, run dates and original output checksum records without personal
  paths. Reversal provenance also preserves the original crossing-analysis
  script and raw input identities where recorded.
- `checksums.sha256`: checksums for the files included under `results/`.

The corrected mean trajectory is column `historical_style_mean` in each
`corrected-summary.tsv`; it includes carried observations. It was verified
against the separately supplied corrected mean table before import, so that
redundant table is not duplicated here. Snapshot and histogram filenames
express their observation role; the snapshot's step column gives its time.

## Verification and regeneration

From the repository root:

```bash
bash scripts/check-released-results.sh
```

The script verifies included-file checksums; rebuilds all twelve corrected
histograms from individual snapshots; checks their observation steps,
extinction associations and summary means; rebuilds the sweep summary from
12100 individual records; and checks sweep aggregate reconstruction and
historical/corrected comparison values. Regeneration uses the existing
Bash/AWK tools and writes temporary data under `build/`.

Numeric comparisons of regenerated tables use an absolute tolerance of
1e-12 for floating-point text differences. Integer counts and regenerated
histogram category totals are represented numerically by the same comparison.
Comparisons with six-decimal simulation aggregates allow 0.00000051.

All four corrected reversal summaries have 10001 observation rows for steps
5000 through 15000. Every snapshot has 1000 realisations. The sweep has
121 parameter pairs with 100 realisations each. Every corrected sweep
realisation survived in this experiment.

## Scope and limitations

The compact data are sufficient inputs for mean-trajectory plots, corrected
snapshot distributions, the existing sweep comparison and extinction plots.
They do not contain historical individual snapshots or the full reversal
population and susceptibility trajectories. The original complete outputs
must be archived separately before claiming independent reconstruction of
all published result reports, including final-outcome uncertainty estimates.
Reversal summaries cannot be recomputed from snapshots alone.

`original-output.sha256` records original full-run output hashes, including
files not supplied in this compact bundle. It is provenance, not the manifest
for verifying this directory. Some included originals have been normalised
to LF text endings, so their original Windows byte hashes may also differ.
Use `checksums.sha256` for the released copies.

These data are simulation reruns. They are not the private historical
workbook observations used to assess the paper's histogram analysis.
A separate reviewed evidence extract for that claim remains a release gate.
One seed per setting does not establish robustness across seeds. Different
implementations consume random streams differently, so matching seeds do
not create paired realisations.

Figure rendering and paper-style plotting scripts are the next preparation
stage. This data bundle does not yet claim a completed one-command figure
workflow or container verification.
