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

The third pattern is the opposite and more interesting: ten *abortus* labelled genomes carry a
1071 bp product that *abortus* should not have. All ten are contig level. One is labelled strain
**2308**, the reference virulent *abortus* strain, which should not produce that band at all — a
genome worth treating as mixed, contaminated or mislabelled rather than as an assay failure.

This is the honest reading throughout: a disagreement can mean a mislabelled submission **or** an
assay limit, and the band level evidence is what separates them. Neither should be assumed.

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

Per genome calls for all 2051 assemblies are in
[`docs/taxid234-calls.tsv`](https://github.com/duceppemo/brucetyper/blob/main/docs/taxid234-calls.tsv):
accession, organism, strain, assembly level, group, and the call from each of the four assays.

The counts on this page are derived from that file, so every one of them can be checked.
