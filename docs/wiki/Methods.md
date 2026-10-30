# How a call is made

## How a call is made

The sizes an assay's profiles use are first merged into **canonical bands**: sizes within
`--tolerance` of each other become one band, because no gel separates them. Each sample is then
reduced to the presence or absence of every canonical band, and that pattern is compared to the
profiles in [`assays.py`](https://github.com/duceppemo/brucetyper/blob/main/assays.py):

- **Bruce-ladder and Suis-ladder** match **exactly** — every canonical band must agree in presence
  and absence.
- **AMOS** and the **RB51 junction PCR** match on a **subset** — the bands of the profile must be
  present, extra bands are tolerated, and the first profile listed that fits wins.

A sample matching no profile is called `N/A`. When an assay genuinely cannot separate two
genotypes, both names are reported joined by `/`.

Repeated loci are collapsed before calling: a locus present several times in a genome yields
several products of the same size, which a gel shows as one band.

All primers and expected band sizes live in `assays.py` and nowhere else. Edit that file to change
how genotypes are called. Every primer has been checked against its publication and every band
size against the product it gives in 87 real genomes: see
[Validation](https://github.com/duceppemo/brucetyper/wiki/Validation).

## The rough flag

Bruce-ladder rows carry a `rough` column, read from one band and one band only: the 1682 bp
product of `BMEI0998_BMEI0997`, which spans *wboA* and *wboB*, the two glycosyltransferases of the
O-antigen locus. `no` means the band is intact, `yes` means it is shifted or missing, and the
column is empty for Suis-ladder, AMOS and the RB51 junction PCR, which do not carry that band. The flag is reported
beside the call and never changes it.

It is a reading of one band, not a phenotype:

- **It misses *B. canis*.** *B. canis* is naturally rough but keeps an intact 1682 bp band, so all
  ten *canis* genomes tested are flagged `no`. Its roughness does not come from losing *wboA*.
- **An incomplete assembly can trip it.** A draft that simply does not contain the locus is
  indistinguishable from one that has lost it. In testing, a 49 contig *neotomae* draft 370 kb
  short of a full genome was flagged `yes` on that basis alone.

So read `yes` as "the O-antigen band is not where it should be, look closer", and do not read `no`
as "smooth". What it does catch is shown in [Validation](https://github.com/duceppemo/brucetyper/wiki/Validation).
