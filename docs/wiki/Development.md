# Development

## Layout

```
brucetyper/
  __init__.py       the version
  assays.py         every primer and expected band size, and nowhere else
  insilico_pcr.py   primer matching and amplicon finding
  genotype.py       sizes to canonical bands to a call
  cli.py            argument parsing and the reports; the brucetyper entry point
bruce_ladder.py     a shim at the root, so python bruce_ladder.py still works
bruce_rates.py      legacy 2023
suis_rates.py       legacy 2023
```

`pip install -e .` for a working install from a clone. See
[Installation](https://github.com/duceppemo/brucetyper/wiki/Installation).

## Tests

```bash
python -m unittest discover -s tests -v
```

The suite runs the CLI as `python -m brucetyper.cli`, so it exercises the same entry point the
installed `brucetyper` command uses, and runs unchanged from a clone or from an unpacked source
distribution. One test covers the `bruce_ladder.py` shim and skips when it is absent, as it is in
the distribution.

The suite builds synthetic genomes whose products reproduce each expected band profile, runs the
CLI over them and checks the calls. It exercises the pipeline end to end; it is **not** a
validation of the assays against real *Brucella* sequence.

## Legacy scripts

`bruce_rates.py` and `suis_rates.py` are the original 2023 scripts, kept as they were apart from
three fixes:

- they no longer hard-code an input path
- `bruce_rates.py` no longer divides by zero on a species that is absent from the report, which
  killed it on its first iteration for any report without an *abortus* genome
- its `wrongPreID` and `notID` counters now match the call against the species list by substring
  rather than by equality, which is what that list holds: `'abort'` has to be found inside
  `abortus`. `wrongPreID` compared the species to the list it had just been removed from, so it
  could never fire at all

A third 2023 script, `suis_ladder.py`, was removed in 1.1.0. It called genotypes from the
consolidated report of a separate in silico PCR tool, a step `brucetyper` now does itself, and it
was never finished: `write_output` was a stub that wrote the literal string `test`, and
`do_amos_pcr` wrote one `-` row per profile that did not match instead of one row per sample. It
remains in the history as `git show e13310a:suis_ladder.py`.

Both `*_rates.py` scripts score calls against species or biovar names embedded in file names, and
take the report to score as their first argument:

```bash
python bruce_rates.py results/report_bruce.tsv
```

## Coverage

Coverage is measured in CI and reported to Codecov. The suite drives the CLI through
`subprocess`, and `brucetyper/cli.py` fans out over a multiprocessing pool, so `.coveragerc` turns
on `parallel` mode and the workflow sets `COVERAGE_PROCESS_START`; without that only the parent
process is measured and the figure reads about 35% instead of the real one.

Every file in the repository is measured. Nothing is excluded from the figure except the suite
itself, so the number is not flattered by an omit list:

| Module | Coverage |
| --- | --- |
| `brucetyper/__init__.py` | 100% |
| `brucetyper/assays.py` | 100% |
| `brucetyper/cli.py` | 100% |
| `brucetyper/genotype.py` | 100% |
| `brucetyper/insilico_pcr.py` | 100% |
| `bruce_ladder.py` (the shim) | 100% |
| `bruce_rates.py` (legacy 2023) | 100% |
| `suis_rates.py` (legacy 2023) | 100% |
| **Total** | **100%** |

This arrived in 1.1.0. Before it the two legacy `*_rates.py` scripts had no tests at all and the
total was 59%, and the pair of bugs in `bruce_rates.py` above is exactly what untested code had
been hiding. `tests/test_rates.py` now runs both scripts the way this page documents and reads
their printed blocks back, including the cases that used to crash.
