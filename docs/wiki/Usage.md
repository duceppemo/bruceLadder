# Usage

## Requirements

- Python 3.7+ (profile order matters, so insertion-ordered dicts are required). CI covers 3.9
  through 3.13; 3.7 and 3.8 are past end of life and are not exercised.
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

`--bruce`, `--suis`, `--amos` and `--rb51` select the assays; with none given, all four are run.

| Option | Default | What it does |
| --- | --- | --- |
| `-t`, `--threads` | all cores | FASTA files processed in parallel |
| `--mismatches` | 2 | Mismatches tolerated per primer |
| `--three-prime` | 3 | Bases at the 3' end of a primer that must match exactly |
| `--tolerance` | 5 | Size tolerance in bp when assigning an amplicon to an expected band |
| `--min-size` | 50 | Smallest product reported, primers included |
| `--max-size` | 3000 | Largest product reported, primers included |

`--tolerance` has a ceiling, and it is well below the gap you might expect. It is bounded by the
closest pair of band sizes in an assay: beyond half that gap the two windows overlap and a size in
the overlap is assigned to the lower band, and beyond the gap itself the two bands merge. The
closest pair is 120 bp for Bruce-ladder, 102 for AMOS and 463 for the RB51 junction assay, but
only **21 bp for Suis-ladder**, between its 278 and 299 bp bands — so `--tolerance 11` already
overlaps those two windows and 21 merges them. Raising it above 10 is not recommended; the
measured cost of widening even a single band is on the
[Large scale validation](https://github.com/duceppemo/brucetyper/wiki/Large-scale-validation#why-a-wider-tolerance-on-the-450-bp-band-was-rejected)
page.

## Output

Three kinds of report are written to `--outdir`:

- **`amplicons.tsv`** — every product found, with its locus, contig, strand, 1-based coordinates,
  size and the mismatch count of each primer.
- **`calls.tsv`** — one row per sample and assay: the call, the rough flag, the canonical bands
  that were present, the observed sizes, and any size belonging to no band of the assay.
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
sample          assay  call             rough  bands         sizes         extra_bands
RB51_or_2308    amos   RB51_or_2308            178;387;494   494;387;178
S19             amos   S19                     494           494
abortus_bv1_2_4 amos   abortus_bv1_2_4         178;494       494;178
melitensis      amos   melitensis              178;733       733;178
ovis            amos   ovis                    178;976       976;178
suis_bv1        amos   suis_bv1                178;285       285;178
```

The `rough` column is empty above because only Bruce-ladder carries the band it reads.
