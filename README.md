# Antimicrobial resistance: an agent-based lattice model

C simulations accompanying R. C. Alamino (2020), *An agent-based lattice
model for the emergence of anti-microbial resistance*, Journal of
Theoretical Biology 486, 110080.
https://doi.org/10.1016/j.jtbi.2019.110080

**Release candidate, stage 1.** The scientific implementations are frozen
for release preparation. Numerical figure inputs, paper-style plots,
container verification and licensing review are not yet complete.
Do not describe this candidate as a finished reproducibility release.

## Build and check

Requirements: GCC with C99 support, Bash, AWK and standard GNU command-line
utilities. Gnuplot with pngcairo is required only for rendering figures.
No Python dependency is introduced.

From the repository root:

```bash
bash scripts/check-all.sh
bash scripts/check-historical.sh
```

The first command builds and tests the corrected implementation. It checks
lattice operations, culling, configuration, output errors, simulation
structure and repeatability, extinction handling, summaries, snapshots,
histograms and parent selection. It does not compare corrected trajectories
to the historical references. The second command builds the isolated
historical baseline and compares it with the saved short references, normalising only CRLF/LF
line endings.

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

Existing figures are diagnostic audit plots. Paper-style versions will be
added separately without changing the scientific calculations. Reported
findings are based on seed 1; they do not establish robustness across seeds
or exact reconstruction of every published figure.
