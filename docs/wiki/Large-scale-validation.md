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

| Group | Genomes | Bruce-ladder `N/A` |
| --- | --- | --- |
| Classical named species | 1612 | 59 (4%) |
| *Ochrobactrum*-derived | 333 | **333 (100%)** |
| Unnamed classical clade (`Brucella sp.`, *amazoniensis*) | 102 | 64 (63%) |
| *B. pecoris* | 4 | 4 (100%) |

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

This is the honest reading throughout: a disagreement can mean a mislabelled submission **or** an
assay limit, and the band level evidence is what separates them. Neither should be assumed — the
27 *melitensis* → `pinnipedialis_ceti` cases above have **not** been checked this way, and a
missing band is the likelier explanation there.

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
python bruce_ladder.py -i genomes/ -o results/
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
