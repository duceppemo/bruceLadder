# Primer and band size validation

Every primer in [`assays.py`](../assays.py) was checked against the sequence printed in the
publication that describes its assay, and every expected band size was checked against the product
the primer pair actually gives in real *Brucella* genomes.

## Method

87 RefSeq genomes were run through `bruce_ladder.py` at its defaults (`--mismatches 2`,
`--three-prime 3`, `--tolerance 5`): up to 10 per species for *B. melitensis*, *B. abortus*,
*B. suis*, *B. canis*, *B. ovis*, *B. ceti*, *B. pinnipedialis*, plus all available *B. neotomae*
and *B. microti*, plus the reference and vaccine strains 16M, 2308, S19, RB51, 1330 and the
*B. suis* biovar references. The full list is in
[`validation-genomes.tsv`](validation-genomes.tsv).

"Observed" below is the product the pair gives in those genomes; `n` counts the genomes showing
the main size. A locus only amplifies in the species it was designed to mark, so `n` is well
below 87 for the discriminating loci.

## Primer sequences

**All 36 primer sequences in `assays.py` match the published sequences exactly.** No
substitutions, no typos. Two apparent discrepancies were checked and resolved in favour of
`assays.py`:

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
