# Large scale validation

Every *Brucella* assembly in GenBank was run through all four assays, to ask three questions the
87 genome panel is too small to answer: does the assay agree with the submitted species label, how
often does it stay silent, and does it ever call something it should not.

## The panel

All **2051** current GenBank assemblies under NCBI taxid 234, downloaded 30 October 2026. GenBank
rather than RefSeq on purpose: every RefSeq assembly derives from a GenBank one, so 2051 is the
superset, and the assemblies RefSeq excludes are disproportionately the contaminated and
mislabelled ones this page is about.

Taxid 234 is no longer one organism. *Ochrobactrum* was folded into *Brucella* in 2020, so the
taxid now spans the classical pathogens and a set of environmental and opportunistic species that
these assays were never designed to see. They are counted separately throughout.

| Group | `group` in the table | Genomes | Bruce-ladder `N/A` |
| --- | --- | --- | --- |
| Classical named species | `classical named` | 1612 | 59 (4%) |
| *Ochrobactrum*-derived | `Ochrobactrum-derived` | 333 | **333 (100%)** |
| Unnamed classical clade: `Brucella sp.`, *amazoniensis* | `unnamed classical-clade` | 102 | 64 (63%) |
| *B. pecoris* | `pecoris` | 4 | 4 (100%) |

## Specificity: the assay stays silent where it should

**All 333 *Ochrobactrum*-derived genomes returned `N/A` on Bruce-ladder. Not one produced a call
of any kind.** The same holds for the four *B. pecoris*. That is the cleanest specificity result
on this page: across 337 genomes inside the same taxid but outside the assay's scope, there were
no false positives.

This is why the `N/A` count must not be read as a failure rate. Most of the silence in the full
2051 is the assay correctly declining to speak.

## Concordance with the submitted label

Among the 1612 genomes carrying a classical species name:

| | Genomes | Share |
| --- | --- | --- |
| Call agrees with the NCBI label | 1493 | 92.6% |
| Call disagrees | 60 | 3.7% |
| No call (`N/A`) | 59 | 3.7% |
| **Concordance among genomes that produced a call** | | **96.1%** |

## The disagreements are mostly one missing band

The 60 disagreements are not scattered. Three patterns account for 48 of them, and in each case
the two profiles differ by exactly one band:

| Label | Called | n | Mechanism |
| --- | --- | --- | --- |
| *melitensis* | `pinnipedialis_ceti` | 27 | The profiles differ **only** by the 450 bp band. **0 of 27** have it |
| *melitensis* | `abortus` | 11 | Differ **only** by the 1071 bp band. **0 of 11** have it |
| *abortus* | `melitensis` | 10 | These **gained** the 1071 bp band. **10 of 10** have it |

The first two are band dropout: lose one product and the pattern lands on a neighbouring profile.
Mostly drafts (24 of the 27 are contig or scaffold level), consistent with the locus falling in a
gap — though two of the eleven *melitensis* → `abortus` genomes are **complete**, where assembly
dropout cannot be the explanation and a primer site variant is more likely.

The third pattern is the opposite, and it was followed up: ten *abortus* labelled genomes carry a
1071 bp product that *abortus* should not have. All ten are contig level, and one is labelled
strain **2308**, the reference virulent *abortus* strain, which should not produce that band at
all.

### All ten are *B. melitensis*

Four independent lines of evidence agree, so these are mislabelled submissions and not an assay
failure:

1. **They carry the 1071 bp `BMEII0843_BMEII0844` band**, which the real 2308 reference
   (`GCA_000054005.1`) does not produce.
2. **They lack AMOS's 494 bp `abortus_IS711` band**, which the real 2308 does produce. Eight of
   the ten return `N/A` on AMOS entirely.
3. **Average nucleotide identity puts every one of them with *melitensis***. Against six complete
   *melitensis* and six complete *abortus* genomes:

| Genome | Deposited strain | mean ANI to *melitensis* | mean ANI to *abortus* |
| --- | --- | --- | --- |
| `GCA_018604785.1` | **2308** | **99.967** | 99.721 |
| `GCA_041861985.1` | BAIQ | **99.967** | 99.720 |
| `GCA_003966035.1` | C-612 | **99.964** | 99.718 |
| `GCA_016091965.1` | S-586 | **99.964** | 99.720 |
| `GCA_009733315.1` | 420 | **99.961** | 99.713 |
| `GCA_009733355.1` | C-577 | **99.961** | 99.712 |
| `GCA_009733515.1` | C-587 | **99.961** | 99.721 |
| `GCA_009733415.1` | I-2 | **99.960** | 99.718 |
| `GCA_016091955.1` | S-595 | **99.960** | 99.719 |
| `GCA_036881765.1` | Cam-1 MGL | **99.894** | 99.645 |

The separation is consistent at about 0.24 percentage points. The control runs the other way: the
genuine 2308 reference sits at 99.84–99.98% against those same *abortus* genomes and only
99.69–99.71% against the *melitensis* ones. *Brucella* species are so close that the window is
narrow, but the grouping is unambiguous and bimodal, with no overlap.

4. **A completely independent assay agrees.** The species-specific real-time PCR multiplex of
   Hoover *et al.* 2026 targets unique core sequences and shares no primer with Bruce-ladder or
   AMOS. Run in silico against `GCA_018604785.1`, it returns the *B. melitensis* target (`BM_10`,
   147 bp, probe binding inside the product) and **not** the *B. abortus* target (`BA`, 75 bp),
   which is the reverse of what the genuine 2308 reference gives.

   | Genome | `Bru_genus` | `BA` | `BM_10` | `BS` |
   | --- | --- | --- | --- | --- |
   | *melitensis* 16M | 93 bp | — | **147 bp** | — |
   | *suis* 1330 | 93 bp | — | — | **102 bp** |
   | *abortus* 2308 (genuine) | 93 bp | **75 bp** | — | — |
   | `GCA_018604785.1` (deposited as 2308) | 93 bp | — | **147 bp** | — |

   Every target fired only in its own species, and a *Brucella* (*Ochrobactrum*) genome was
   negative for all four, so the panel was behaving correctly on the controls.

The genome deposited as *B. abortus* **2308** is the clearest case: it is *B. melitensis*, and
anyone using it as the 2308 reference is using the wrong organism.

So **10 of the 60 disagreements are confirmed errors in the submitted labels, not errors in the
assay**. Scoring those as assay failures understates it: of the 1553 classically named genomes
that produced a call, the assay is demonstrably right in at least 1503.

### The 27 *melitensis* → `pinnipedialis_ceti` cases are the opposite: the assay is wrong

Checked the same way, every one of the 27 is genuinely *B. melitensis*: **99.80–99.97%** ANI to
complete *melitensis* genomes against **99.55–99.65%** to complete *ceti* and *pinnipedialis*
genomes, a consistent gap of about 0.3 points with no exceptions. Several are well known
*melitensis* strains, 16M, M5, M111 and 63/9 among them.

So the labels are right and the **call is wrong**. These are assay failures, not mislabelled
submissions, and the cause is always the missing 450 bp `BMEI0535_BMEI0536` band:

| Cause | n |
| --- | --- |
| Both primer sites present but on different contigs — assembly break | 14 |
| Product present but **shifted** off 450 bp | 12 |
| Only one primer site survives | 1 |

The shifted products measure 341, 356, 426, 437, 441, 476, 612, 612, 627, 627, 781 and 781 bp.
The default `--tolerance` is 5 bp, so none of them register as the 450 bp band, and losing that
one band moves the pattern onto the `pinnipedialis_ceti` profile, which differs from *melitensis*
by exactly that band.

This matters beyond these 27: the 87 genome panel saw this locus vary only between 449 and 452 bp
in non-marine strains, so **the locus is considerably more length-variable in *melitensis* than
that panel suggested**. Four of the twelve shifted products sit within 25 bp of 450, close enough
that a larger tolerance would recover them. That was implemented and measured against this panel,
and it does not work.

#### Why a wider tolerance on the 450 bp band was rejected

A per-band tolerance override was written and run over all 2051 genomes, widening only the 450/451
band and leaving every other band of every assay at the default 5 bp. The in silico PCR does not
depend on the tolerance, so the measurement re-reads the stored amplicon sizes: at the default it
reproduces the published calls in
[`docs/taxid234-calls.tsv`](https://github.com/duceppemo/brucetyper/blob/main/docs/taxid234-calls.tsv)
exactly, for all 2051 genomes, which makes every difference below attributable to the widened band
alone.

| Window | Genomes recovered | Correct calls destroyed | Net |
| --- | --- | --- | --- |
| ±5 bp (default) | 0 | 0 | 0 |
| ±10 bp | 1 | 0 | +1 |
| ±15 bp | 2 | 0 | +2 |
| ±20 bp | 3 | 0 | +3 |
| **±24 bp** | **4** | **0** | **+4** |
| **±25 bp** | **5** | **1** | **+4** |
| ±30 bp | 6 | 1 | +5 |
| ±40 bp | 6 | 2 | +4 |

At ±24 bp the window recovers three of the 27 — the 426, 437 and 441 bp products — plus one
unrelated genome, *B. abortus* bv. 2 `GCA_000157735.1`, whose 470 bp product had left it `N/A`.
Nothing is lost.

One more base pair is where it breaks. At ±25 bp the window reaches the fourth shifted product,
the 476 bp one in *melitensis* `GCA_003514625.1`, and in the same step **destroys a correct call
on `GCA_000590795.1`, *B. ceti* TE10759-12 — a complete genome, which gives a product of exactly
476 bp at the same primer pair**. The two products are the same size, so no choice of window
admits one and excludes the other.

Both genomes were checked against five complete *melitensis* and five complete *ceti* and
*pinnipedialis* references, by fastANI and independently by k-mer ANI. The two methods agree, and
the controls separate by the same margin and in the expected direction:

| Genome | Deposited as | mean ANI to *melitensis* | mean ANI to *ceti*/*pinnipedialis* | Gap |
| --- | --- | --- | --- | --- |
| `GCA_003514625.1` | *B. melitensis* | **99.886** | 99.670 | **+0.217** |
| `GCA_000590795.1` | *B. ceti* TE10759-12 | 99.664 | **99.860** | **−0.196** |
| `GCA_047324115.1` (control) | *B. ceti* | 99.625 | **99.885** | −0.260 |
| `GCA_002191335.1` (control) | *B. melitensis* | **99.958** | 99.679 | +0.279 |

So the recovery is genuine and so is the loss, and at ±25 bp they cancel: the same net as ±24 bp,
bought by trading a correct call for an incorrect one on a complete genome. Widening further to
±30 bp does show a better net, but it still destroys that *ceti* call, and its sixth recovery is a
*B. canis* genome, `GCA_001078335.1`, that is a separate matter from the 27.

#### The *ceti* allele at this locus is not confined to the large products

The collision is not bad luck, and the reason it was not anticipated is worth recording. The
`pinnipedialis_ceti` profile is built on the *absence* of a 450 bp band, on the understanding that
`BMEI0535_BMEI0536` gives a much larger product in the marine mammal species. Across the panel
that mostly holds, but not always:

| Product in *ceti* / *pinnipedialis* | n | Genomes |
| --- | --- | --- |
| 1049–1675 bp | 13 | the expected range |
| 762 bp | 1 | `GCA_050815205.1` |
| **476 and 487 bp** | **2** | `GCA_000590795.1` and `GCA_000590815.1`, both **complete** |

The two short ones are TE10759-12 and TE28753-12, a pair of complete *B. ceti* genomes, and they
agree with each other. They put the *ceti* minimum at this locus at **476 bp, not 1049**, close
enough to 450 that any window wide enough to rescue the shifted *melitensis* products reaches them
too.

The 450 bp band is therefore doing two jobs that cannot both be made more forgiving: separating
*melitensis* from the marine mammal species, and absorbing *melitensis* length variation. The
default 5 bp tolerance keeps the first one correct, which is the job the profile depends on, and
the missed calls among the 27 are the price. **No per-band tolerance was adopted.**

## Reading the disagreements as a whole

All 60 disagreements have been checked against ANI:

| | n |
| --- | --- |
| **The assay is wrong** — the label is right, a band was lost or shifted | 41 |
| **The label is wrong** — the submitted name is not supported by the sequence | 19 |

So of the 1553 classically named genomes that produced a call, the call matches the organism in
**1512, or 97.4%**: the 1493 that agree with the label outright, plus 19 where the assay
contradicted the label and was right to.

A disagreement can mean a mislabelled submission **or** an assay limit, and only sequence evidence
separates the two. Here it fell both ways, which is why neither should be assumed.

### Independent corroboration

Thirteen *Brucella* assemblies were independently reported as misnamed on the basis of MLST-9
(PubMLST *brucella* scheme), the species-specific real-time PCR of Hoover *et al.* 2026 scored in
silico, and Mash distances over a dereplicated reference set — none of which shares a method or a
primer with Bruce-ladder.

**Bruce-ladder had independently flagged 12 of the 13 from band patterns alone**, and the
thirteenth is not a contradiction:

| Accession | Deposited as | Bruce-ladder | Reported as |
| --- | --- | --- | --- |
| `GCA_041861985.1` | *B. abortus* | `melitensis` | *B. melitensis* |
| `GCA_016091955.1` | *B. abortus* | `melitensis` | *B. melitensis* |
| `GCA_016091965.1` | *B. abortus* | `melitensis` | *B. melitensis* |
| `GCA_018604785.1` | *B. abortus* (strain 2308) | `melitensis` | *B. melitensis* |
| `GCA_000250775.2` | *B. melitensis* S66 | `suis_bv2/suis_bv1_3_4` | *B. suis* |
| `GCA_000250835.2` | *B. melitensis* 16M13W | `abortus` | *B. abortus* |
| `GCA_000292125.1` | *B. suis* S2-30 | `melitensis` | *B. melitensis* |
| `GCA_055389085.1` | *B. melitensis* | `canis` | *B. canis* |
| `GCA_009823695.1` | *B. melitensis* | `suis_bv2/suis_bv1_3_4` | *B. suis* |
| `GCA_000292045.1` | *B. abortus* 104M | `suis_bv2/suis_bv1_3_4` | *B. suis* |
| `GCA_055388155.1` | *B. melitensis* | `suis_bv2/suis_bv1_3_4` | *B. suis* |
| `GCA_055389475.1` | *B. melitensis* | `abortus` | *B. abortus* |
| `GCA_000662015.2` | *B. inopinata* BO1 | `N/A` | classical, suis-canis group |

The last one is the assay declining rather than disagreeing: it produces no Bruce-ladder pattern
at all, which is the correct behaviour for a genome it cannot place.

### Where ANI is not enough

The species assignments above come from MLST and the real-time PCR assay, not from the ANI used
elsewhere on this page, because **ANI cannot separate *B. suis* from *B. canis*** at the distances
involved. Measured against three references per species:

| Assembly | ANI to *canis* | ANI to *suis* | Gap |
| --- | --- | --- | --- |
| `GCA_000292045.1` | 99.894 | 99.891 | **0.003** |
| `GCA_000250775.2` | 99.891 | 99.884 | **0.007** |
| `GCA_009823695.1` | 99.925 | 99.872 | 0.053 |

For comparison, the *melitensis* against *abortus* separation elsewhere on this page is about
0.25 points, and a known *B. canis* reference sits 0.118 above *suis*. A gap of 0.003 carries no
information. *B. canis* is a lineage within the *B. suis* group, so *suis* is the more general and
better supported reading for all three, and the nine-locus MLST profiles place them with deposited
*B. suis* outright.

The same limit applies to *suis* against *microti*: a known *B. suis* reference sits only 0.011
above *microti*. That is the same indistinguishability Bruce-ladder itself reports by calling
`suis_bv2/suis_bv1_3_4` rather than naming one.

## `N/A` is not an assembly quality gradient

| Assembly level | `N/A` rate among classical |
| --- | --- |
| Complete genome | 4.6% (10/216) |
| Chromosome | 3.3% (5/150) |
| Scaffold | 9.6% (18/187) |
| Contig | 2.5% (26/1059) |

Scaffold level assemblies fail most often, but complete genomes fail nearly twice as often as
contig level ones, so fragmentation alone does not explain silence. Ten complete genomes carrying
a classical name produced no Bruce-ladder call at all.

## The RB51 assay at scale

Across all 2051 genomes the `--rb51` junction assay called **RB51 exactly twice**:

| Accession | Organism | Strain | Bruce-ladder | AMOS |
| --- | --- | --- | --- | --- |
| `GCA_011801185.1` | *B. abortus* RB51-AHVLA | RB51 | `RB51` | `RB51_or_2308` |
| `GCA_052512005.2` | *B. abortus* | RB-51 C | `RB51` | `RB51_or_2308` |

Both are genuine RB51 vaccine strains, and both were independently confirmed by Bruce-ladder's
shifted *wboA* band. There was **no disagreement in either direction** across the whole panel: no
genome called `RB51` by one assay and not the other.

AMOS, as expected, could not separate either of them from strain 2308, returning `RB51_or_2308` —
which is the gap the junction assay was added to close.

The 430 `N/A` results from this assay are the genomes with no readable *wboA* locus, which is the
assay declining to answer rather than answering wrongly.

## Vaccine strains found

| Call | Assay | n |
| --- | --- | --- |
| `RB51` | Bruce-ladder and `--rb51` | 2 |
| `S19` | Bruce-ladder | 13 |
| `S19` | AMOS | 5 |

Bruce-ladder and AMOS disagree on S19 because they read different features: Bruce-ladder calls it
on the missing 590 bp band, AMOS on the absence of the 178 bp *eri* band that the S19 deletion
removes. The two criteria do not pick out the same genomes, and neither is being treated here as
ground truth.

## Reproducing this

```bash
# the accession list is the group column of the results table
brucetyper -i genomes/ -o results/
```

The real-time PCR cross-check above is not part of brucetyper, which reads amplicon sizes the way
a gel is read and cannot represent a hydrolysis probe. Probe assays are a first-class mode in
[insilicoPCR](https://github.com/duceppemo/insilicoPCR); the primers and probes are published in
Hoover E, Chmara J, Duceppe M-O, *et al.* Two novel real-time PCR assays for *Brucella* detection.
*Microbiol Spectr*. 2026;14:e0297425.
[PMID 41885434](https://pubmed.ncbi.nlm.nih.gov/41885434/).

Per genome calls for all 2051 assemblies are in
[`docs/taxid234-calls.tsv`](https://github.com/duceppemo/brucetyper/blob/main/docs/taxid234-calls.tsv):
accession, organism, strain, assembly level, group, and the call from each of the four assays.

The counts on this page are derived from that file, so every one of them can be checked.
