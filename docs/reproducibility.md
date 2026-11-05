# Reproducibility

The code and checks run from a source archive without Git metadata.
Use `bash scripts/check-all.sh` for corrected checks and
`bash scripts/check-historical.sh` for historical short-reference checks.
Historical sources and corrected sources are built into separate paths.

## Native setup (no Docker)

Requirements are a C99 compiler (the scripts invoke `gcc`), C development
headers and the maths library, Bash, AWK and GNU utilities including
Coreutils, Diffutils, Findutils, Grep, Sed, Tar and Gzip. Git supplies optional
provenance metadata. Rendering also needs Gnuplot with the `pngcairo` and
`svg` terminals. Paper-style plots request the Liberation Serif font.
Python is not required.

### Debian or Ubuntu Linux

Install dependencies using your system package manager:

```bash
sudo apt-get update
sudo apt-get install --no-install-recommends \
  bash gcc libc6-dev gawk coreutils diffutils findutils grep sed \
  git tar gzip gnuplot-nox fonts-liberation
```

These commands install system packages and may request your administrator
password. Other Linux distributions use their own package names and managers.

### Windows with MSYS2

Install Git for Windows (including Git Bash) and MSYS2. In the MSYS2 UCRT64
terminal, update the installation with `pacman -Syu`; reopen the terminal
and repeat the update if instructed. Install the compiler and Gnuplot:

```bash
pacman -S --needed mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-gnuplot
```

For the default MSYS2 installation location, run the following in Git Bash:

```bash
export PATH="/c/msys64/ucrt64/bin:$PATH"
gcc --version
gnuplot --version
```

Adapt the path if MSYS2 is installed elsewhere. Install Liberation Serif as a
Windows font for the intended typography; missing fonts may be substituted.
The Docker route includes this font. Font substitution does not change the
numerical data, but can affect labels and spacing.

### Build and validate

From the repository root:

```bash
bash scripts/check-all.sh
bash scripts/check-historical.sh
bash scripts/check-released-results.sh
```

These compile the programs and run short tests. They do not rerun the full
production experiments. To build without running the checks, use
`bash scripts/build.sh` and `bash scripts/build-historical.sh`.

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
with Git Bash and MSYS2. Docker/Linux validation also passed with Debian 12,
GCC 12.2.0 and Gnuplot 5.4.4 on amd64. Representative PNG layouts were
visually inspected after UTF-8 and font settings were configured.

## Figure Reproduction

Existing PNGs can be viewed, and compact numerical figure inputs are now
included under `results/`. Run `bash scripts/check-released-results.sh` to
verify them and regenerate histograms and the sweep summary. Full reversal
trajectories are not included, so trajectory summaries and historical
final-outcome reports are not independently reconstructed by this command.
The paper-style generation command is now available:
`bash scripts/render-paper-figures.sh`. See
[figure definitions and rendering status](paper-style-figures.md).
Native Windows rendering and representative Docker/Linux PNG layouts have
been visually inspected. SVG files are generated but have not received the
same visual verification. No private files are required by the code or test workflow. The public archived-data evidence
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

### Recorded experiments in Docker

Use a named container so that its run records can be copied out after it
finishes. For example, the following is a full corrected 500-step treatment
experiment and may take hours:

```bash
docker run --name amr-experiment amr-repro \
  bash scripts/run-experiment.sh corrected reversal 1 1000 500
mkdir -p build/docker-experiments
docker cp amr-experiment:/amr/build/experiments/. ./build/docker-experiments/
docker rm amr-experiment
```

In Windows Git Bash, prefix `docker cp` with `MSYS_NO_PATHCONV=1`.
Only remove the container after successfully copying the run records.
Use a different container name for concurrent runs. To run the historical
sweep example above, replace the runner arguments with
`historical sweep 1 100`. Interrupted simulations have no checkpoint resume;
a new invocation starts a new run. Validate outputs before interpretation.

### Preserving a tested image

A package list records the installed versions but does not lock future
builds. An optional local image archive preserves the installed environment:

```bash
mkdir -p build/docker-validation
docker run --rm amr-repro dpkg-query -W '-f=${binary:Package}\t${Version}\n' \
  > build/docker-validation/packages.tsv
docker image save -o build/docker-validation/amr-repro.tar amr-repro
(cd build/docker-validation && sha256sum amr-repro.tar > amr-repro.tar.sha256)
```

Verify and reload an archive on a compatible Docker host using:

```bash
(cd build/docker-validation && sha256sum -c amr-repro.tar.sha256)
docker image load -i build/docker-validation/amr-repro.tar
```

Local archives under `build/` are not distributed by the Git repository.
An archive is specific to the image's processor architecture. Loading it
avoids rebuilding and package downloads, but still requires Docker. Native
execution remains a separate supported route.
