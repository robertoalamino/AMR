# Reproducibility

The code and checks run from a source archive without Git metadata.
Use `bash scripts/check-all.sh` for corrected checks and
`bash scripts/check-historical.sh` for historical short-reference checks.
Historical sources and corrected sources are built into separate paths.

## Recorded experiments

```bash
bash scripts/run-experiment.sh corrected reversal 1 1000 500
bash scripts/run-experiment.sh historical sweep 1 100
```

The runner copies the selected local source tree into a unique run directory,
builds that copy and records source hashes, compiler information, commands,
exit status and output hashes. It uses the local source bytes, including
uncommitted edits; commit or archive changes when preparing an official run.
Git commit metadata is supplementary, while `source.sha256` identifies
what was actually built. Run records are local and may contain absolute
paths; publish only reviewed portable provenance records.

Full experiments can take hours. A successful process exit is not scientific
validation. Check dimensions, ranges, parameter association and aggregate
reconstruction before interpreting output.

## Numerical portability

Historical reference comparisons normalise CRLF to LF and then require exact
content. No numerical tolerance is applied. Agreement
on another platform must be checked rather than assumed. Corrected checks
use structural properties, controlled cases and same-executable repeatability.
Cross-platform floating-point and image identity are not promised.

The corrected suite and historical regression checks pass on Linux
with GCC 14.2; the development workflow has also been exercised on Windows
with Git Bash and MSYS2. Container and plotting verification are pending.

## Figure Reproduction

Existing PNGs can be viewed, and compact numerical figure inputs are now
included under `results/`. Run `bash scripts/check-released-results.sh` to
verify them and regenerate histograms and the sweep summary. Full reversal
trajectories are not included, so trajectory summaries and historical
final-outcome reports are not independently reconstructed by this command.
The paper-style generation command is now available:
`bash scripts/render-paper-figures.sh`. See
[figure definitions and rendering status](paper-style-figures.md).
Actual rendering and visual inspection remain required before release. No private files are
required by the code or test workflow. The public archived-data evidence
extract and the full-result archive must be completed before release.

## Optional Docker Environment

Native execution remains supported. Docker is an alternative that packages the compiler, analysis utilities, Gnuplot and fonts.

Docker must be installed and its engine running. From the repository root:

```bash
docker build -t amr-repro .
docker run --rm amr-repro
docker run --rm amr-repro bash scripts/check-historical.sh
docker run --rm amr-repro bash scripts/check-released-results.sh
```

The default command runs the short corrected-implementation checks. These commands do not rerun the full production experiments.

To regenerate figures from the included result tables:

```bash
docker run --name amr-render amr-repro bash scripts/render-paper-figures.sh
mkdir -p build/paper-figures-docker
docker cp amr-render:/amr/build/paper-figures/. ./build/paper-figures-docker/
docker rm amr-render
```

In Windows Git Bash, prefix the `docker cp` command with
`MSYS_NO_PATHCONV=1` to prevent conversion of the container's Linux path.

The image contains a copy of the repository taken during the build.
Rebuild after changing source files or plotting scripts.

The Debian base image is pinned by digest. Additional packages come from
maintained Debian repositories, so their versions may change between builds.
This configuration passed the corrected checks, historical regressions and
released-result checks on Linux/amd64 through Docker Desktop with WSL 2.
Representative PNGs from each paper-style plot layout were visually inspected.

These checks do not establish identical full-production trajectories across
all operating systems and processor architectures.
