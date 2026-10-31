# Installation

## Requirements

- Python 3.7+ (profile order matters, so insertion-ordered dicts are required). CI covers 3.9
  through 3.13; 3.7 and 3.8 are past end of life and are not exercised.
- [Biopython](https://biopython.org/), which `pip` installs for you.

Nothing else. No external PCR tool, no aligner, no database to download.

## pip

```bash
pip install brucetyper
```

That puts a `brucetyper` command on your path:

```bash
brucetyper --version
brucetyper -i genomes/ -o results/
```

## conda

```bash
conda install -c conda-forge -c bioconda brucetyper
```

`mamba` resolves the same thing considerably faster:

```bash
mamba install -c conda-forge -c bioconda brucetyper
```

Or set the channels up once, which is what Bioconda recommends, and drop the `-c` flags from
every later install:

```bash
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install brucetyper
```

The recipe is `noarch`, so the one build serves every platform conda supports.

## pipx

A command line tool in its own isolated environment, which keeps Biopython out of your user
site-packages:

```bash
pipx install brucetyper
brucetyper -i genomes/ -o results/
```

## From a clone

For development, or to run a version that has not been released:

```bash
git clone https://github.com/duceppemo/brucetyper
cd brucetyper
pip install -e .
```

`-e` installs in place, so edits to `brucetyper/assays.py` take effect without reinstalling, which
is the usual reason to work from a clone.

Without installing at all, a clone also runs directly, which is how every example in this wiki
worked before the package existed:

```bash
pip install biopython
python bruce_ladder.py -i genomes/ -o results/
```

`bruce_ladder.py` is a three line shim at the repository root that dispatches to
`brucetyper.cli`. It takes exactly the same arguments as the `brucetyper` command and is kept so
that the older invocation, and the `bruce_rates.py` and `suis_rates.py` workflows, keep working.
It is not part of the installed package.

## Checking the install

```bash
$ brucetyper --version
brucetyper 1.1.0
```

The test suite is shipped inside the source distribution, so a release can be verified without a
clone:

```bash
pip download --no-binary :all: --no-deps brucetyper
tar xzf brucetyper-*.tar.gz
cd brucetyper-*/
python -m unittest discover -s tests
```

One test is skipped there: the shim lives at the repository root and is not in the distribution.

## What is installed

| | |
| --- | --- |
| `brucetyper` | The command, equivalent to `python bruce_ladder.py` from a clone |
| `brucetyper.assays` | Every primer and expected band size. The only place either is defined |
| `brucetyper.insilico_pcr` | Primer matching and amplicon finding |
| `brucetyper.genotype` | Sizes to canonical bands to a call |

`brucetyper.assays` is importable, so the profiles can be read or patched from your own code:

```python
from brucetyper import assays
print(assays.PROFILES['bruce']['RB51'])
```

See [Usage](https://github.com/duceppemo/brucetyper/wiki/Usage) for the options and the reports,
and [Methods](https://github.com/duceppemo/brucetyper/wiki/Methods) for how a call is made.
