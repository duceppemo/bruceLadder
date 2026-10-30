# Primer and band size validation

Every primer in [`assays.py`](https://github.com/duceppemo/brucetyper/blob/main/assays.py) was checked against the sequence printed in the
publication that describes its assay, and every expected band size was checked against the product
the primer pair actually gives in real *Brucella* genomes.

## Method

87 RefSeq genomes were run through `bruce_ladder.py` at its defaults (`--mismatches 2`,
`--three-prime 3`, `--tolerance 5`): up to 10 per species for *B. melitensis*, *B. abortus*,
*B. suis*, *B. canis*, *B. ovis*, *B. ceti*, *B. pinnipedialis*, plus all available *B. neotomae*
and *B. microti*, plus the reference and vaccine strains 16M, 2308, S19, RB51, 1330 and the
*B. suis* biovar references. The full list is in
[`validation-genomes.tsv`](https://github.com/duceppemo/brucetyper/blob/main/docs/validation-genomes.tsv), which also records the call each genome
received from each of the four assays, so every count in this file can be checked against it.

"Observed" below is the product the pair gives in those genomes; `n` counts the genomes showing
the main size. A locus only amplifies in the species it was designed to mark, so `n` is well
below 87 for the discriminating loci.

## Primer sequences

**All 40 primer entries in `assays.py`, 33 distinct sequences, match the published sequences
exactly.** No substitutions, no typos. Entries outnumber sequences because several are reused: the
universal IS711 primer serves five AMOS loci, the `BMEI1426_BMEI1427` pair serves both Suis-ladder
and Bruce-ladder, and primer 3 serves both loci of the RB51 junction assay. The three sequences
added later for that assay were read from Fig. 1 of Vemulapalli *et al.* 1999 and are covered
under [RB51 junction PCR](#rb51-junction-pcr) below.
Two apparent discrepancies were checked and resolved in favour of `assays.py`:

| Primer | Printed elsewhere as | Resolution |
| --- | --- | --- |
| `BMEI0536` R | `CGCAGGGGAAAACAGCTATAA` | The `GGGG` form has **no** exact match in any genome; the form in `assays.py` matches and gives the published 450 bp. The other spelling is a typo in that secondary table. |
| `BMEI1688` F | `TCAATCGCGTGAACAATGCT` | That form has **no** exact match in 16M; the form in `assays.py` matches and gives the published 278 bp. |

One assay-level difference is **not** a typo and is described under Bruce-ladder below.

## Bruce-ladder

Published sizes from García-Yoldi *et al.* 2006 / López-Goñi *et al.* 2008.

| Locus | Published | In `assays.py` | Observed in silico | n | Note |
| --- | --- | --- | --- | --- | --- |
| `BMEI0998_BMEI0997` | 1682 | 1682 | 1682 | 74 | 2526 bp in RB51, where IS711 has inserted into *wboA* |
| `BMEI0535_BMEI0536` | 450 | 450 / 451 | 450 (50), 451 (14) | 50 | 449–452 across strains; 1049–1675 in *ceti*/*pinnipedialis* |
| `BMEII0843_BMEII0844` | 1071 | 1071 | 1071 | 71 | 1077 and 1094 in two *canis* |
| `BMEII0428_BMEII0428` | 587 | **590** | **590** | 85 | Published size is 3 bp low; 590 is unanimous across all 85 |
| `BR0953` | 272 | 272 | 272 | 37 | |
| `BMEI0752` | 218 | — | **no product** | 0 | See below |
| `BMEII0987` | 152 | 152 | 152 | 79 | 151 in two *ceti* |
| `BMEI1436_BMEI1435` | 794 | *not used* | 794 | — | See below |
| `BMEI1426_BMEI1427` | *not a Bruce-ladder locus* | 774 | 774 | 76 | See below |

### The fourth primer pair is the Suis-ladder locus, not the Bruce-ladder one

Published Bruce-ladder uses `BMEI1436F`/`BMEI1435R` (`ACGCAGACGACCTTCGGTAT` /
`TTTATCCATCGCCCTGTCAC`) for a 794 bp band. `assays.py` instead carries `BMEI1426`/`BMEI1427`
(774 bp), which is the Suis-ladder locus, and its `bruce` profiles encode 774 to match.

Both pairs were tested across the panel. The published pair gives exactly 794 bp, and the pair in
use gives exactly 774 bp; neither amplifies in *B. canis*. Because the profiles are written around
774 consistently, calls come out right — but this is a substituted locus, not an adjusted size,
and the assay as implemented is not the published Bruce-ladder at this position.

### BMEI0752 never amplifies

The `BMEI0752` pair yielded **no product in any of the 87 genomes**. The forward primer has one
exact match per genome; the reverse primer, `GATGTGGTAACGCACACCAA`, has none. The genomes carry
`GATGTGGTAACGCACACCCG` at that site — the published primer differs in its **last two bases**, at
the 3' end, which `--three-prime 3` requires to match exactly. Relaxing to `--three-prime 0`
recovers a 219 bp product, against a published 218.

This is not a transcription error in `assays.py`: the published tables print
`GATGTGGTAACGCACACCAA`, and the genome-matching spelling appears in no paper indexed by Europe PMC.
The published primer itself mismatches *Brucella* at its 3' end. The consequence is that
Bruce-ladder here runs on seven of its eight loci, which is consistent with 218 bp being absent
from every profile.

### Vaccine strains

`RB51` and `S19` are called from the two bands that separate them from abortus, and were added
after this validation:

| Genotype | Profile | Separated from abortus by | Confirmed in |
| --- | --- | --- | --- |
| `RB51` | 2526, 774, 590, 450, 152 | 1682 shifts to 2526, IS711 inside *wboA* | GCF_011801185.1 |
| `S19` | 1682, 774, 450, 152 | the 590 bp BMEII0428 band is absent | GCF_000018725.1, GCF_038420195.1 |

Adding them changed three calls in the 87 genome panel, all from `N/A` to the right vaccine
strain, and no other call. Strain 2308 and eight of the nine *abortus* field strains still call `abortus`, the ninth being
LBAB009, which calls `N/A` for the unrelated reason given below,
so the vaccine profiles do not shadow it. Note that 2308 is reported as `abortus` by this assay,
unlike AMOS, which groups it with RB51.

### The rough flag on the panel

The `rough` column of `calls.tsv` reads the 1682 bp *wboA* band alone. Across the 87 genomes it
flags 13:

| Flagged | Genomes | Why |
| --- | --- | --- |
| *B. ovis* | 10 of 10 | Naturally rough, no *wboA* product at all |
| RB51 | 1 | IS711 inside *wboA*, band shifted to 2526 bp |
| LBAB009 | 1 | ~15 kb deletion removing the locus, see below |
| *neotomae* NCTC10070 | 1 | **False positive**: a 49 contig draft, 2.94 Mb, ~370 kb short of a full genome. The locus is not in the assembly rather than not in the strain |

And the one it misses:

| Not flagged | Genomes | Why it matters |
| --- | --- | --- |
| *B. canis* | 10 of 10 | **False negative**: *B. canis* is naturally rough, but keeps an intact 1682 bp band. Its roughness does not come from losing *wboA*, so a band-level rule cannot see it |

`no` therefore does not mean smooth. It means this one band is where it should be. Of the two
naturally rough species, the flag catches *B. ovis* and misses *B. canis*.

### The one unexplained N/A: LBAB009

`abortus__LBAB009__GCF_030761395.1` is the only genome in the panel that calls `N/A` for a reason
other than the assay's own limits. It gives 774, 590, 450 and 152 but no 1682 bp *wboA* band.

It is not an assembly artefact: the genome is complete, two chromosomes, 3.29 Mb. Neither
`BMEI0998_BMEI0997` primer has any match in it, because roughly 15 kb of chromosome 1 is missing,
corresponding to 2308 `NC_007618.1` 953,840..968,920. That block carries `BAB_RS20745` and
`BAB_RS20750`, the two glycosyltransferases of the O-antigen locus (*wboA*, *wboB*), so the strain
is most likely a naturally rough isolate. Its sibling isolates LBAB008 and LBAB016 both carry the
band, and AMOS calls LBAB009 `abortus_bv1_2_4` normally.

`N/A` is the right answer here. The strain is *B. abortus* carrying a real deletion, and no profile
was added for it: calling it `abortus` would hide the finding, and a profile fitted to one isolate
would not generalise. It differs from RB51, where IS711 inserts into *wboA* and shifts the band to
2526 bp rather than removing it.

## Suis-ladder

Published sizes from López-Goñi *et al.* 2011.

| Locus | Published | In `assays.py` | Observed in silico | n | Note |
| --- | --- | --- | --- | --- | --- |
| `BMEI1426_BMEI1427` | 774 | 774 | 774 | 76 | |
| `BR1080` | 197 | 197 | 197 | 43 | |
| `BMEI1688_MEI1687` | 278 | 278 | 278 | 57 | |
| `BMEI0205` | 299 / 425 / 551 / 614 | 299 / 425 / 551 / 614 | 299 (19), 614 (13), 551 (12), 425 (11) | 55 | VNTR; also 173 (11, *ovis*), 236 (9, *abortus*/*melitensis*), 803 (5, *microti*), 549 (3) |

All five biovar reference strains were called correctly: `bv01` through `bv05`.

## AMOS

Published sizes from Bricker and Halling 1994 and 1995.

| Locus | Published | In `assays.py` | Observed in silico | n | Note |
| --- | --- | --- | --- | --- | --- |
| `abortus_IS711` | 498 | **494** | **494** | 10 | Published size is 4 bp high; 494 in all 10 *abortus*, S19 and RB51 included |
| `melitensis_IS711` | 731 | **733** | **732** (8), 733 (2) | 10 | Published size is 1–2 bp low; 733 is the 16M reference value, 732 the field-strain value. Both inside `--tolerance` |
| `ovis_IS711` | 976 | 976 | 976 | 9 | |
| `suis_IS711` | 285 | 285 | 285 (2), 276 (1) | 3 | 285 confirmed in 1330; amplifies in very few genomes |
| `RB51_2308_IS711` | 364 | **387** | **387** | 10 | Published size is 23 bp low. 387 in 2308 and RB51, and in 7 *ovis* and 1 *pinnipedialis*, which carry their own tandem IS711 pair. 2112 bp in three *ceti* |
| `eri` | 178 | 178 | 178 | 85 | Present in every genome but S19 |

## RB51 junction PCR

Primers from Vemulapalli *et al.* 1999, Fig. 1. Validated over the same 87 genomes with `--rb51`.

### The three primers

| | Sequence | Role |
| --- | --- | --- |
| Primer 3 | `GCCAACCAACCCAAATGCTCACAA` | Forward, in *wboA* upstream of the insertion point. Shared by both pairs |
| Primer 1 | `TTAAGCGCTGATGCCATTTCCTTCAC` | Reverse, in *wboA* downstream of the insertion point |
| Primer 2 | `TTTAGTTTGCCGTAATATAGGT` + `CTAGAACCTGTC` | Reverse, **chimeric**: 22 bases of *wboA* fused to 12 bases of IS711 |

Primer 2 is what makes the assay specific. Its *wboA* half has an exact site in each of the 75
genomes that carry the locus and in none of the other 12. Its IS711 half is the **terminal** 12 bp
of the element's right end, not an internal stretch, so its site count does not track IS711 copy
number: one to three exact sites per genome, 14 genomes having only one, and in 2308 only one of
its two sites is a real IS711 end, the other a coincidental 12-mer. The **fused** 34-mer has an
exact site in RB51 and in no other genome. It cannot anneal unless IS711 is inside *wboA*.

The whole margin lives in primer 2's 3' terminal 8 bases. In every *wboA* bearing non-RB51 genome
the target reads `TTCAAGCG|CTAG|ACCTATATTACGGCAAACTAAA` against the primer's
`GACAGGTT|CTAG|ACCTATATTACGGCAAACTAAA`: the 22 base *wboA* half and the first 4 bases of the
IS711 half already match, and all five mismatches fall in the last 8 bases, the 3' terminal base
included. The mismatch budget and the 3' anchor therefore block it independently.

### The insertion

Aligning the RB51 and 2308 *wboA* products places the insertion 342 bp into the Bruce-ladder
amplicon, at RB51 `NZ_CP046720.1:964,987`, spanning 964,987–965,830. The insertion adds **844 bp**,
which is the whole difference between the two sizes of the `wboA` locus below. IS711 is usually
quoted as 842 bp in the literature; 844 bp is what these assemblies carry, and a second full
length copy in 2308 (`NC_007618.1:541,972`) measures 844 bp as well.

The breakpoint is ambiguous over a 4 bp window from sequence alone, because the *wboA* target and
the element both end in `CTAG`: a naive common prefix and suffix walk returns 964,983 and an
840 bp insert. 964,987 is the placement at which the insert begins with the conserved IS711 left
end `GGCGTGTCTGCA`, which is what fixes it.

### Band sizes

| Locus | Primers | In `assays.py` | Observed in silico | n | Note |
| --- | --- | --- | --- | --- | --- |
| `wboA` | 3 + 1 | 456 | 456 (74), 1300 (1) | 74 | 1300 bp is RB51, 456 + the 844 bp insertion. It is deliberately not a profile size: it is tolerated as an extra band so that an RB51 draft missing it still calls RB51. No product in the 12 genomes without the locus |
| `wboA_IS711` | 3 + 2 | 919 | **919** | **1** | RB51 only, 0 mismatches on both primers. No product in the other 86 genomes, 2308 and S19 included |

The sizes Vemulapalli *et al.* print for these products are **not** recorded here, unlike every
other assay in this file: the figure gives the primer sequences but no sizes, and the sizes quoted
for this assay in secondary sources could not be reconciled with any genome, so only the observed
values are stated. See [Two sets of primers that are not this assay](#two-sets-of-primers-that-are-not-this-assay).

### Calls over the 87 genomes

| Genomes | Call | Why |
| --- | --- | --- |
| RB51 | `RB51` | 919 bp junction band, plus 1300 bp in `extra_bands` |
| 74 others | `not_RB51` | *wboA* intact at 456 bp, no junction band |
| 10 *B. ovis* | `N/A` | Naturally no *wboA* product at all |
| *abortus* LBAB009 | `N/A` | O-antigen locus deleted, see above |
| *neotomae* NCTC10070 | `N/A` | 49 contig draft that does not contain the locus |

No false positive, and every call the assay declined to make was a genome with no *wboA* locus to
read rather than a silent failure.

### What the mismatch budget is carrying

The `N/A` set is the set of genomes in which primer 3 has no site **within the tolerated mismatch
budget**, which is not the same as having no exact site. Primer 3 has no exact site in **19** of
the 87 genomes: the 12 above, plus seven *B. suis* biovar 2 genomes that carry
`GCCAACCA`**`C`**`CCCAAATGCTCACAA`, an A→C transversion at primer position 9. Primer 1 matches
those seven exactly, so at the default `--mismatches 2` they give a clean 456 bp product and call
`not_RB51`; at `--mismatches 0` they fall to `N/A` and the panel splits 19 / 67 / 1 instead of
12 / 74 / 1. The SNP sits far from the 3' end, so the `--three-prime` anchor does not block it.

In the other direction the assay has room to spare. The next best placement of primer 2 anywhere
in the panel, on either strand and with no 3' anchoring, is **5 mismatches**, so the calls are
unchanged from `--mismatches 2` through `--mismatches 6`. That 5 mismatch site is real and sits
75 bp from the primer 3 site in all 74 *wboA* bearing non-RB51 genomes — it is the intended target
minus the element. Forced through with `--mismatches 8 --three-prime 0`, it yields a 75 bp product
in all 74, and the calls still come out 1 / 74 / 12, because 75 bp is nowhere near the 919 bp
band.

### Two sets of primers that are not this assay

Both were checked against the genomes and rejected, which is why the sequences above were taken
from the figure rather than from the surrounding text:

| Sequences | Verdict |
| --- | --- |
| `TTATCGCAATCCATCGCGAT`, `AGTGCGACGGAGAGTCCCAA`, quoted for this assay by secondary sources | **Not *Brucella*.** The best placement of either anywhere in RB51 or 2308 is 4 mismatches, and 6 or more in the O-antigen region itself |
| `GGATGTCGACCAGCCCTCCACATCAATAGC`, `TTGCGGATCCTTTACTCGTCCGTCTCTTAC`, from the paper's methods section | **Genuine but not diagnostic.** They match exactly once their engineered SalI and BamHI tails are set aside, but they are the cloning primers used to characterise the insertion: they flank the whole gene and give 2114 bp in 2308 against 2958 bp in RB51, both strains amplifying |

## Summary of differences

Sizes in `assays.py` that deliberately depart from the published figure, each confirmed by the
genomes:

| Assay | Locus | Published | Used | Difference |
| --- | --- | --- | --- | --- |
| amos | `RB51_2308_IS711` | 364 | 387 | +23 |
| amos | `abortus_IS711` | 498 | 494 | −4 |
| bruce | `BMEII0428_BMEII0428` | 587 | 590 | +3 |
| amos | `melitensis_IS711` | 731 | 733 | +2 |

Open points, neither of which changes a call today:

- `BMEI1426_BMEI1427` (774 bp) stands where published Bruce-ladder uses `BMEI1436_BMEI1435`
  (794 bp).
- `BMEI0752` cannot amplify at the default 3' anchor, so the published 218 bp band is unreachable.
