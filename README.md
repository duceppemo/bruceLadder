<p align="center">
  <img src="docs/logo.svg" alt="" width="150">
</p>

<h1 align="center">brucetyper</h1>

<p align="center">
  <a href="https://github.com/duceppemo/brucetyper/actions/workflows/tests.yml"><img src="https://github.com/duceppemo/brucetyper/actions/workflows/tests.yml/badge.svg" alt="Tests"></a>
  <a href="https://codecov.io/gh/duceppemo/brucetyper"><img src="https://codecov.io/gh/duceppemo/brucetyper/graph/badge.svg" alt="Coverage"></a>
  <img src="https://img.shields.io/badge/python-3.7%2B-blue" alt="Python 3.7+">
  <a href="LICENSE"><img src="https://img.shields.io/github/license/duceppemo/brucetyper" alt="License"></a>
  <a href="https://github.com/duceppemo/brucetyper/wiki"><img src="https://img.shields.io/badge/docs-wiki-informational" alt="Documentation"></a>
  <a href="https://doi.org/10.5281/zenodo.23815138"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.23815138.svg" alt="DOI"></a>
</p>

In silico **Bruce-ladder**, **Suis-ladder**, **AMOS PCR** and an **RB51 junction PCR** for
*Brucella*, run end to end from assembled genomes.

For every FASTA file the primer pairs of each requested assay are located directly in the
nucleotide sequence, the resulting amplicon sizes are read the way a gel would be read, and a
genotype is called. No external PCR tool, no alignment step, no database.

| Assay | Flag | Separates |
| --- | --- | --- |
| Bruce-ladder | `--bruce` | All *Brucella* species, plus the S19 and RB51 vaccine strains |
| Suis-ladder | `--suis` | *B. suis* biovars 1 to 5, and *suis* from *canis* |
| AMOS PCR | `--amos` | *abortus* bv 1/2/4, *melitensis*, *ovis*, *suis* bv 1, S19 |
| RB51 junction | `--rb51` | RB51 from its parent 2308, which no other assay here does |

## Quick start

```bash
git clone https://github.com/duceppemo/brucetyper
cd brucetyper
pip install biopython

# Every assay, on a folder of assemblies
python bruce_ladder.py -i genomes/ -o results/

# One assay, on a single genome
python bruce_ladder.py -i GCF_000054005.1.fna.gz -o results/ --amos
```

Input is a FASTA file or a folder of them, gzipped or not, one sample per file. Assemblies, not
reads. The sample name is the file name with its extension stripped. With no assay flag, all four
are run.

Three reports are written to `--outdir`: `amplicons.tsv` (every product found), `calls.tsv` (one
row per sample and assay) and `report_<assay>.tsv` (sample and call, two columns).

```
$ cat results/calls.tsv
sample  assay  call          rough  bands                     sizes                 extra_bands
RB51    bruce  RB51          yes    152;450-451;590;774;2526  2526;774;590;450;152
RB51    amos   RB51_or_2308         178;387;494               494;387;178
RB51    rb51   RB51                 919                       1300;919              1300
2308    bruce  abortus       no     152;450-451;590;774;1682  1682;774;590;450;152
2308    rb51   not_RB51             456                       456
```

## Documentation

Everything else is in the [wiki](https://github.com/duceppemo/brucetyper/wiki), whose sources are
maintained in [`docs/wiki`](docs/wiki):

- [Usage](https://github.com/duceppemo/brucetyper/wiki/Usage) — requirements, inputs, every option,
  and the three reports
- [Methods](https://github.com/duceppemo/brucetyper/wiki/Methods) — how sizes become canonical
  bands, how a profile is matched, how the rough flag is read
- [The assays](https://github.com/duceppemo/brucetyper/wiki/Assays) — what each resolves, what it
  cannot, the RB51 junction assay in detail, and the publications
- [Validation](https://github.com/duceppemo/brucetyper/wiki/Validation) — every primer checked
  against its publication and every band size against 87 real genomes, including the primer sets
  that were checked and **rejected**
- [Large scale validation](https://github.com/duceppemo/brucetyper/wiki/Large-scale-validation) —
  all **2051** GenBank assemblies under taxid 234: 96.1% concordance with the submitted species
  label among genomes that produced a call, and **zero** calls on the 333 *Ochrobactrum*-derived
  genomes
- [Development](https://github.com/duceppemo/brucetyper/wiki/Development) — the test suite and the
  legacy 2023 scripts

All primers and expected band sizes live in [`assays.py`](assays.py) and nowhere else. Edit that
file to change how genotypes are called.

## Citation

If brucetyper is useful in your work, please cite it together with the assay publications listed
on the [wiki](https://github.com/duceppemo/brucetyper/wiki/Assays#assay-sources):

> Duceppe, M.-O. brucetyper: in silico Bruce-ladder, Suis-ladder, AMOS and RB51 junction PCR for
> *Brucella*. https://doi.org/10.5281/zenodo.23815138

The DOI above always resolves to the latest release; see
[`CITATION.cff`](CITATION.cff).

## License

[MIT](LICENSE). Copyright held by His Majesty the King in Right of Canada, as represented by the
Canadian Food Inspection Agency.
