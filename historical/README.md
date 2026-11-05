# Historical comparison baseline

These source files are copied without modification from the supplied
`pre-scientific-corrections` archive, revision
`fb82dd0ed027256c1dfe512314fa0a31dcdf06da`.

This baseline includes maintenance and refactoring. It is not the untouched
publication implementation, and exact figure provenance is not implied.
It retains known scientific discrepancies for comparison only.

Build with `bash scripts/build-historical.sh`. Executables are isolated in
`build/historical/`; corrected executables remain in `build/`.
Run `bash scripts/check-historical.sh` for both historical regressions.
Reference outputs are in `tests/reference/`.

`source.sha256` records the imported source bytes. From the repository root:

```bash
(cd historical && sha256sum -c source.sha256)
```
