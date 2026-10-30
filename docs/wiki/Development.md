# Development

## Tests

```bash
python -m unittest discover -s tests -v
```

The suite builds synthetic genomes whose products reproduce each expected band profile, runs the
CLI over them and checks the calls. It exercises the pipeline end to end; it is **not** a
validation of the assays against real *Brucella* sequence.

## Legacy scripts

`suis_ladder.py`, `bruce_rates.py` and `suis_rates.py` are the original 2023 scripts, kept as they
were except that the two `*_rates.py` scripts no longer hard-code an input path. `suis_ladder.py` called genotypes from the consolidated report of a separate in silico PCR
tool, a step `bruce_ladder.py` now does itself. The two `*_rates.py` scripts score calls against
species or biovar names embedded in file names, and take the report to score as their first
argument:

```bash
python bruce_rates.py results/report_bruce.tsv
```

## Coverage

Coverage is measured in CI and reported to Codecov. The suite drives the CLI through
`subprocess`, and `bruce_ladder.py` fans out over a multiprocessing pool, so `.coveragerc` turns
on `parallel` mode and the workflow sets `COVERAGE_PROCESS_START`; without that only the parent
process is measured and the figure reads about 35% instead of the real one.

The headline number is held down by the legacy scripts, which have no tests:

| Module | Coverage |
| --- | --- |
| `assays.py` | 100% |
| `bruce_ladder.py` | 100% |
| `genotype.py` | 100% |
| `insilico_pcr.py` | 99% |
| `bruce_rates.py`, `suis_rates.py`, `suis_ladder.py` (legacy 2023) | 0% |
| **Total** | **58%** |

The legacy scripts are deliberately left in the measurement rather than excluded. Omitting them
would put the badge near 100% and overstate what is tested.
