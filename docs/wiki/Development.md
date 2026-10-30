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
