# Antimicrobial resistance: an agent-based lattice model

C simulations accompanying R. C. Alamino (2020), *An agent-based lattice
model for the emergence of anti-microbial resistance*, Journal of
Theoretical Biology 486, 110080.
https://doi.org/10.1016/j.jtbi.2019.110080

**Release candidate.** The scientific implementations are frozen for release
preparation. Compact result data, checks and paper-style plots are included.
Native Windows and Linux checks and Docker/Linux checks have passed, and
representative rendered PNGs have been visually inspected. Licensing and
remaining archival release work are not yet complete.

## Choose a reproduction route

- **Native:** install the compiler and command-line tools yourself, then run
  the scripts below. Docker is not required.
- **Docker:** build the supplied environment and run the same scripts inside
  it. Docker must be installed and its engine running.

Both routes use the public code and included data; no private workbooks or
external research folders are needed. See the
[setup and reproduction instructions](docs/reproducibility.md).

Viewing the committed figures needs no software installation. Regenerating
figures from included tables does not rerun simulations. Full production
experiments can take hours.

## Build and check

Requirements: GCC with C99 support, Bash, AWK and standard GNU command-line
utilities. Gnuplot with pngcairo is required only for rendering figures.
No Python dependency is introduced. Native setup commands for Linux and Windows
are provided in [Reproducibility](docs/reproducibility.md).

From the repository root:

```bash
bash scripts/check-all.sh
bash scripts/check-historical.sh
bash scripts/check-released-results.sh
```

The first command builds and tests the corrected implementation. It checks
lattice operations, culling, configuration, output errors, simulation
structure and repeatability, extinction handling, summaries, snapshots,
histograms and parent selection. It does not compare corrected trajectories
to the historical references. The second command builds the isolated
historical baseline and compares it with the saved short references, normalising only CRLF/LF
line endings. The third command verifies the compact released result data
and regenerates their snapshot histograms and sweep summary.

## Programs

- `amdspace`: 121-pair antimicrobial parameter sweep.
- `reversaltimes`: environmental adaptation, treatment and withdrawal.

Corrected executables are in `build/`; historical executables are in
`build/historical/`. Each simulation writes fixed filenames in its current
working directory. Use a separate directory for every run.

```bash
mkdir -p build/example-reversal
(cd build/example-reversal && ../reversaltimes.exe \
  --steps 120 --realisations 2 --start 40 --end 60 --seed 1)
```

The reversal writes `population.txt`, `qav.txt` and `extinction.txt`.
The sweep writes `asymdrug.txt`, `extinction.txt` and
`sweep-realisations.txt` in the corrected implementation.

## Documentation

- [Model and output conventions](docs/model.md)
- [Corrections and scientific scope](docs/corrections.md)
- [Reproducibility](docs/reproducibility.md)
- [Experiment settings](experiments/README.md)
- [Released results status](results/README.md)
- [Historical baseline](historical/README.md)
- [Sweep findings](docs/results-sweep.md)
- [Reversal findings](docs/results-reversal-summary.md)

Diagnostic audit plots and paper-style figures are both included. Their
presentation differs without changing the scientific calculations. Reported
findings are based on seed 1; they do not establish robustness across seeds
or exact reconstruction of every published figure.

## Paper-style figures

```bash
bash scripts/render-paper-figures.sh
```

This generates PNG and SVG versions under `build/paper-figures/`, using
verified result tables without rerunning simulations. See
[figure definitions and verification status](docs/paper-style-figures.md).
The plots use paper-style presentation of identified rerun data, not the
original published images. Inspect generated outputs before publication.

## Docker quick start

```bash
docker build -t amr-repro .
docker run --rm amr-repro
docker run --rm amr-repro bash scripts/check-historical.sh
docker run --rm amr-repro bash scripts/check-released-results.sh
```

These run short checks, not full production experiments. See
[Reproducibility](docs/reproducibility.md) for rendering, retrieving outputs
and running recorded experiments through either route.
