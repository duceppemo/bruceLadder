# brucetyper

In silico **Bruce-ladder**, **Suis-ladder**, **AMOS PCR** and an **RB51 junction PCR** for
*Brucella*, run end to end from assembled genomes.

For every FASTA file the primer pairs of each requested assay are located directly in the
nucleotide sequence, the resulting amplicon sizes are read the way a gel would be read, and a
genotype is called. No external PCR tool, no alignment step, no database.

## Pages

- **[Installation](https://github.com/duceppemo/brucetyper/wiki/Installation)** — pip, pipx, from a clone, and what gets installed
- **[Usage](https://github.com/duceppemo/brucetyper/wiki/Usage)** — inputs, every command line option, and the three reports that
  come out
- **[Methods](https://github.com/duceppemo/brucetyper/wiki/Methods)** — how observed sizes become
  canonical bands, how a profile is matched, and how the rough flag is read
- **[The assays](https://github.com/duceppemo/brucetyper/wiki/Assays)** — what each assay resolves,
  what it cannot, the RB51 junction assay in detail, and the publications
- **[Validation](https://github.com/duceppemo/brucetyper/wiki/Validation)** — every primer checked
  against its publication and every band size against 87 real genomes, including the primer sets
  that were checked and rejected
- **[Large scale validation](https://github.com/duceppemo/brucetyper/wiki/Large-scale-validation)** — all 2051 GenBank assemblies under
  taxid 234 through all four assays: concordance with the submitted labels, what the
  disagreements actually are, and zero calls on the 333 *Ochrobactrum*-derived genomes
- **[Development](https://github.com/duceppemo/brucetyper/wiki/Development)** — the test suite and
  the legacy 2023 scripts

## At a glance

| Assay | Flag | Separates |
| --- | --- | --- |
| Bruce-ladder | `--bruce` | All *Brucella* species, plus the S19 and RB51 vaccine strains |
| Suis-ladder | `--suis` | *B. suis* biovars 1 to 5, and *suis* from *canis* |
| AMOS PCR | `--amos` | *abortus* bv 1/2/4, *melitensis*, *ovis*, *suis* bv 1, S19 |
| RB51 junction | `--rb51` | RB51 from its parent 2308, which no other assay here does |

These pages are generated from
[`docs/wiki`](https://github.com/duceppemo/brucetyper/tree/main/docs/wiki) in the repository. Edit
them there, not in the wiki, or the next sync will overwrite your changes.
