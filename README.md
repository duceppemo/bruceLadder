<p align="center">
  <img src="docs/logo.svg" alt="" width="150">
</p>

<h1 align="center">bruceLadder</h1>

In silico **Bruce-ladder**, **Suis-ladder** and **AMOS PCR** for *Brucella*, run end to end from
assembled genomes.

For every FASTA file the primer pairs of each requested assay are located directly in the
nucleotide sequence, the resulting amplicon sizes are read the way a gel would be read, and a
genotype is called. No external PCR tool, no alignment step, no database.

## Requirements

- Python 3.7+ (profile order matters, so insertion-ordered dicts are required)
- [Biopython](https://biopython.org/)

```bash
pip install biopython
```

## Usage

```bash
# Every assay, on a folder of assemblies
python bruce_ladder.py -i genomes/ -o results/

# One assay, on a single genome
python bruce_ladder.py -i GCF_000054005.1.fna.gz -o results/ --amos
```

Input is a FASTA file or a folder of them, gzipped or not, one sample per file. Recognised
extensions are `.fasta`, `.fa`, `.fna`, `.fas`, `.seq` and `.ffn`, each optionally `.gz`. The
sample name is the file name with those extensions stripped. Assemblies, not reads.

`--bruce`, `--suis` and `--amos` select the assays; with none given, all three are run.

| Option | Default | What it does |
| --- | --- | --- |
| `-t`, `--threads` | all cores | FASTA files processed in parallel |
| `--mismatches` | 2 | Mismatches tolerated per primer |
| `--three-prime` | 3 | Bases at the 3' end of a primer that must match exactly |
| `--tolerance` | 5 | Size tolerance in bp when assigning an amplicon to an expected band |
| `--min-size` | 50 | Smallest product reported, primers included |
| `--max-size` | 3000 | Largest product reported, primers included |

## Output

Three kinds of report are written to `--outdir`:

- **`amplicons.tsv`** — every product found, with its locus, contig, strand, 1-based coordinates,
  size and the mismatch count of each primer.
- **`calls.tsv`** — one row per sample and assay: the call, the canonical bands that were present,
  the observed sizes, and any size belonging to no band of the assay.
- **`report_<assay>.tsv`** — sample and call, two columns, no header. This is the layout
  `bruce_rates.py` and `suis_rates.py` expect.

```
$ python bruce_ladder.py -i genomes/ -o results/ --amos
Found 6 FASTA file(s)
Running in silico PCR for assay(s): amos
  RB51_or_2308: 3 amplicon(s)
  ...
Calling genotypes
Reports written to /path/to/results

$ cat results/calls.tsv
sample          assay  call             bands         sizes         extra_bands
RB51_or_2308    amos   RB51_or_2308     178;387;494   494;387;178
S19             amos   S19              494           494
abortus_bv1_2_4 amos   abortus_bv1_2_4  178;494       494;178
melitensis      amos   melitensis       178;733       733;178
ovis            amos   ovis             178;976       976;178
suis_bv1        amos   suis_bv1         178;285       285;178
```

## How a call is made

The sizes an assay's profiles use are first merged into **canonical bands**: sizes within
`--tolerance` of each other become one band, because no gel separates them. Each sample is then
reduced to the presence or absence of every canonical band, and that pattern is compared to the
profiles in [`assays.py`](assays.py):

- **Bruce-ladder and Suis-ladder** match **exactly** — every canonical band must agree in presence
  and absence.
- **AMOS** matches on a **subset** — the bands of the profile must be present, extra bands are
  tolerated, and the first profile listed that fits wins.

A sample matching no profile is called `N/A`. When an assay genuinely cannot separate two
genotypes, both names are reported joined by `/`.

Repeated loci are collapsed before calling: a locus present several times in a genome yields
several products of the same size, which a gel shows as one band.

All primers and expected band sizes live in `assays.py` and nowhere else. Edit that file to change
how genotypes are called. Every primer has been checked against its publication and every band
size against the product it gives in 87 real genomes: see
[`docs/band-sizes.md`](docs/band-sizes.md).

## What the assays cannot resolve

These are properties of the assays, not of the implementation:

- **Bruce-ladder**: *B. microti* has the same profile as *B. suis* bv 2 and is not listed
  separately. *B. suis* bv 2 and bv 1/3/4 differ by a single bp (450 vs 451), which no gel and no
  size tolerance separates, so both are reported as `suis_bv2/suis_bv1_3_4`.
  *B. pinnipedialis*/*B. ceti* sometimes shows one extra band of roughly 1100 to 1700 bp.
- **AMOS**: the RB51 primer lies inside IS711 itself and pairs with the universal IS711 primer
  across the junction of two tandem elements, so the band reports that arrangement rather than the
  strain. Vaccine strain RB51 and its parent 2308 are indistinguishable here and are reported as
  `RB51_or_2308`. *B. ovis* carries a tandem IS711 pair of its own and gives the same product; the
  abortus band is what keeps the two apart. S19 is called on the absence of the 178 bp *eri* band,
  which its deletion removes.

## Tests

```bash
python -m unittest discover -s tests -v
```

The suite builds synthetic genomes whose products reproduce each expected band profile, runs the
CLI over them and checks the calls. It exercises the pipeline end to end; it is **not** a
validation of the assays against real *Brucella* sequence.

## Legacy scripts

`suis_ladder.py`, `bruce_rates.py` and `suis_rates.py` are the original 2023 scripts, kept as they
were. `suis_ladder.py` called genotypes from the consolidated report of a separate in silico PCR
tool, a step `bruce_ladder.py` now does itself. The two `*_rates.py` scripts score calls against
species or biovar names embedded in file names and have their input paths hard-coded at the top.

## Assay sources

**Bruce-ladder**

- García-Yoldi D, Marín CM, de Miguel MJ, Muñoz PM, Vizmanos JL, López-Goñi I. Multiplex PCR assay
  for the identification and differentiation of all *Brucella* species and the vaccine strains
  *Brucella abortus* S19 and RB51 and *Brucella melitensis* Rev1. *Clin Chem*.
  2006;52(4):779-81. [PMID 16595839](https://pubmed.ncbi.nlm.nih.gov/16595839/)
- López-Goñi I, García-Yoldi D, Marín CM, de Miguel MJ, Muñoz PM, Blasco JM, *et al.* Evaluation of
  a multiplex PCR assay (Bruce-ladder) for molecular typing of all *Brucella* species, including
  the vaccine strains. *J Clin Microbiol*.
  2008;46(10):3484-7. [PMID 18716225](https://pubmed.ncbi.nlm.nih.gov/18716225/)

**Suis-ladder**

- López-Goñi I, García-Yoldi D, Marín CM, de Miguel MJ, Barquero-Calvo E, Guzmán-Verri C, Albert D,
  Garin-Bastuji B. New Bruce-ladder multiplex PCR assay for the biovar typing of *Brucella suis*
  and the discrimination of *Brucella suis* and *Brucella canis*. *Vet Microbiol*.
  2011;154(1-2):152-5. [PMID 21782356](https://pubmed.ncbi.nlm.nih.gov/21782356/)

**AMOS**

- Bricker BJ, Halling SM. Differentiation of *Brucella abortus* bv. 1, 2, and 4, *Brucella
  melitensis*, *Brucella ovis*, and *Brucella suis* bv. 1 by PCR. *J Clin Microbiol*.
  1994;32(11):2660-6. [PMID 7852552](https://pubmed.ncbi.nlm.nih.gov/7852552/)
- Bricker BJ, Halling SM. Enhancement of the *Brucella* AMOS PCR assay for differentiation of
  *Brucella abortus* vaccine strains S19 and RB51. *J Clin Microbiol*.
  1995;33(6):1640-2. [PMID 7650203](https://pubmed.ncbi.nlm.nih.gov/7650203/)
