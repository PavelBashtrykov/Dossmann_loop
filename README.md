# Dosmann et. al

Two bash scripts that take raw single-end Illumina reads of a deep enzymology library
and turn them into a flanking-sequence / methylation report.

```
raw FASTQ ──► clean_data.sh ──► <...>_unique_seq.txt.gz ──► analyse_flanks.sh ──► report.txt
```

Both scripts are plain bash + awk. The only external programs they use are
`zcat`/`gzip` (and `pigz` when it is installed, for faster compression).

---

## Pipeline overview

### Stage 1 — `clean_data.sh`

Cleans the raw FASTQ and returns a list of unique read sequences.

| Step | What it does | Controlled by |
|------|--------------|---------------|
| 1 | **Hard trimming** — trim randomized nucleotides at the start and the end of a read, keep bases `-c`..`-k` of every read (a fixed window). Reads shorter than the window are discarded. | `-c`, `-k` |
| 2 | **Quality filtering / adapter screening** — discard reads containing the adapter (anywhere, allowing one substitution, insertion or deletion) and reads that fail the Phred `-q` cutoff. | `-q`, `-Q`, `-a` |
| 3 | **Read selection** — keep only reads that start with the forward barcode and end with the reverse barcode, each joined with the constant part of the substrate next to it (see [Anchors](#anchors-selecting-on-the-constant-parts-of-the-substrate)). Skipped when neither barcode is given. | `-f`, `-r` |
| 4 | **Barcode trimming** — cut `-b` bases off the 5′ end and `-e` bases off the 3′ end. Runs together with step 3. | `-b`, `-e` |
| 5 | **Duplicate removal** — removes reads with an identical full-length sequence, keeping the first record of each. Skipped with `-n`. | `-n` |

Steps 1 and 2 run as a single pass over the input. Quality trimming uses the
running-sum rule of cutadapt/BWA: a read is discarded when that rule would
shorten it below the full window (the window length is never reduced — reads
either pass intact or are dropped). With `-Q` the strict rule applies instead:
one base below the `-q` cutoff anywhere in the window discards the read.

**Outputs** (all in the output directory, one `.fq.gz` per step):

```
<prefix>_q<Q>_<win>bp.fq.gz              after steps 1+2
<prefix>_q<Q>_<win>bp_selected.fq.gz     after step 3
<prefix>_q<Q>_selected_<len>bp.fq.gz     after step 4
<prefix>_q<Q>_selected_<len>bp_unique.fq.gz       after step 5 (FASTQ)
<prefix>_q<Q>_selected_<len>bp_unique_seq.txt.gz  after step 5 (one sequence per line)
<prefix>.pipeline.log                    the shared pipeline log (see below)
```

The `*_unique_seq.txt.gz` file is the input of stage 2.

### Stage 2 — `analyse_flanks.sh`

Takes the plain sequence list and, for every read, reads out the methylation
state of a central motif together with its sequence context.

Each read carries **both strands** of the substrate:
the random part of the top strand starts at `-1`, the random part of the
bottom strand at `-3`; each is `-l` + motif + `-l` nucleotides long.

Per read:

1. **Motif lookup** on both strands, accepting all bisulfite states
   (for CG: `CG` or `TG` on top, `CG` or `TG` on bottom). Reads whose motif
   matches none of `-m` are discarded.
2. **Reconstitution** of the original, unconverted sequence of both random
   parts from the two strands, position by position. The central motif keeps
   its bisulfite state, so it still reports methylation.
3. **Mismatch filter** — reads with more than `-x` positions that no base pair
   matches are discarded.
4. **Dcm filter** — reads with a Dcm site (CCWGG, i.e. `CCAGG`, `CCTGG`,
   `CTAGG` or `CTTGG` on the reconstituted top strand) around the motif are
   discarded, unless `-D` is given.
5. **Counting** — reads are split into methylated and unmethylated by the
   state of the motif on the analysed strand `-s`, and the script counts
   - the frequency of every base at every position of the random part;
   - the frequency of every flanking context, for flanks of 1 to `-f`
     nucleotides on each side.


---

## Quick start

```bash
# Stage 1: clean one lane, select one sample by its barcodes
./clean_data.sh -i WT.fq.gz \
                -o WT_R1_pipeline \
                -Q \
                -f TCC -r GCT \
                -p 8

# Stage 2: flanking sequence analysis of the unique reads
./analyse_flanks.sh -i WT_R1_pipeline/..._unique_seq.txt.gz \
                    -o WT_R1_pipeline/WT_R1_flanks.txt \
                    -g WT_R1_pipeline/<prefix>.pipeline.log
```

`-g` appends the counts of the flank analysis to the log that `clean_data.sh`
wrote, so both stages end up in one file.

---

## `clean_data.sh` — all options

```
Usage: clean_data.sh -i <STR> [-o <STR>] [-p <INT>] [-q <INT>] [-Q] [-c <INT>] [-k <INT>]
                     [-f <STR>] [-r <STR>] [-b <INT>] [-e <INT>] [-a <STR>] [-n] [-h]
```

| Flag | Argument | Default | Meaning |
|------|----------|---------|---------|
| `-i` | STR | — (**required**) | `/path/to/reads.fq.gz` — single-end FASTQ, plain or gzipped. |
| `-o` | STR | `./<prefix>_YYMMDD_HHMM` | Output directory. The default carries a run timestamp. |
| `-p` | INT | `2` | Threads used for compression. Only effective when `pigz` is installed; with plain `gzip` the value is ignored. |
| `-q` | INT | `20` | Quality trimming cutoff (Phred) for step 2. |
| `-Q` | flag | off (running sum) | Strict quality filtering in step 2: a read is discarded as soon as **one** base of the kept window falls below the `-q` cutoff, instead of the running-sum rule of cutadapt/BWA. Reads are never shortened either way. The output files of steps 2–4 carry `q<Q>strict` instead of `q<Q>`, so a strict run never overwrites a running-sum run in the same directory. |
| `-c` | INT | `7` | First base to keep (1-based) in step 1. |
| `-k` | INT | `135` | Last base to keep (1-based) in step 1. The defaults keep positions 7..135 = **129 nt**. |
| `-b` | INT | `3` | Bases trimmed from the 5′ end in step 4. |
| `-e` | INT | `3` | Bases trimmed from the 3′ end in step 4. |
| `-f` | STR | unset | Forward barcode — the read must *start* with it, followed by the 5′ constant part of the substrate. |
| `-r` | STR | unset | Reverse barcode — the read must *end* with the 3′ constant part of the substrate followed by it. Both barcodes are matched exactly, as given, on the read as it is after step 2; reverse-complement the reverse barcode yourself if your library requires it. Only `A C G T N` are allowed. The constant parts are added by the script, see [Anchors](#anchors-selecting-on-the-constant-parts-of-the-substrate). **Step 3 (and with it step 4) is skipped when neither `-f` nor `-r` is given.** |
| `-a` | STR | `AGATCGGAAGAGC` | Adapter screened for in step 2 (Illumina default). A read is discarded when the adapter occurs anywhere in it, allowing one substituted, missing or inserted base. |
| `-n` | flag | off (dedup runs) | Skip step 5, do **not** remove PCR duplicates. |
| `-h` | flag | — | Print help and exit. |

All parameters, per-step counts, file sizes and runtimes go into
`<prefix>.pipeline.log` in the output directory; nothing but the final path is
printed on screen.

Step 5 keeps the first record of every read sequence seen (exact, full-length
match) and reports the duplicate fraction in the log. It writes
`<last step output>_unique.fq.gz` and `<last step output>_unique_seq.txt.gz`
(one unique read sequence per line, same order as the FASTQ).

Defaults in combination: window 7..135 = 129 nt → minus 3 + 3 barcode bases
= **123 nt**, which is exactly the read length `analyse_flanks.sh` expects.

### Anchors: selecting on the constant parts of the substrate

Step 3 does not match the barcodes alone. The substrate carries 24 constant
nucleotides on each side of the two random parts, directly inside the barcodes:

```
5′ <-f> GAGTGTGATTAGGTTTTTATTGTT ... random ... GGTAGTGAGAGTTTAGTTATATTT <-r> 3′
```

Those two constant parts are built into the script, as `ANCHOR_5` and
`ANCHOR_3` at the top of `clean_data.sh`. On the command line you still pass
the 3 nt barcodes only — the script joins each barcode with the constant part
next to it and selects on the result, so every read is anchored on the barcode
plus 24 more nucleotides on that side:

```bash
./clean_data.sh -i AJ053_..._L5_1.fq.gz -o WT_R1_pipeline -Q \
                -f TCC -r GCT \
                -p 8
# STEP 3 then matches  TCCGAGTGTGATTAGGTTTTTATTGTT  at the 5′ end
#                and   GGTAGTGAGAGTTTAGTTATATTTGCT  at the 3′ end
```

---

## `analyse_flanks.sh` — all options

```
Usage: analyse_flanks.sh -i <STR> [-o <STR>] [-g <STR>] [-m <STR>] [-1 <INT>] [-3 <INT>] [-l <INT>]
                         [-s <STR>] [-L <INT>] [-x <INT>] [-f <INT>] [-D] [-h]
```

| Flag | Argument | Default | Meaning |
|------|----------|---------|---------|
| `-i` | STR | — (**required**) | `/path/to/sequences.txt(.gz)` — one read per line, LF or CRLF line endings. This is the `*_unique_seq.txt.gz` of stage 1. |
| `-o` | STR | `./<prefix>_flanks_<strand>_YYMMDD_HHMM.txt` | Report file. The default carries a run timestamp; a value given here is used verbatim. |
| `-g` | STR | unset | Shared pipeline log written by `clean_data.sh`. The counts and parameters of this run are appended to it as *STEP 6*, together with a run metadata footer. The file is created if it does not exist. Counts still go to stderr as well. |
| `-m` | STR | `CG,CA,CT,CC` | Central motifs, comma separated, each starting with `C`. A read is assigned to the first motif it matches. |
| `-1` | INT | `25` | First base of the **top**-strand random part, 1-based. |
| `-3` | INT | `78` | First base of the **bottom**-strand random part, 1-based. |
| `-l` | INT | `10` | Flank length of the random part on each side of the motif. |
| `-s` | STR | `top` | Strand to analyse: `top` or `bottom`. Methylation must be readable on it — `CA`, `CT` and `CC` carry no cytosine on the bottom strand, and the script aborts if a requested motif/strand combination cannot report methylation. |
| `-L` | INT | `123` | Expected read length; reads of any other length are discarded (counted as *length is not N nt*). `0` disables the check — reads are then only dropped when the random parts do not fit into them (*random parts do not fit*). |
| `-x` | INT | `0` | Mismatch limit — the number of positions of the random part for which no base pair on the complementary strand matches. |
| `-f` | INT | `3` | Largest flank reported in the flank tables, 1 to 5 (and never larger than `-l`). The default produces tables for 1, 2 and 3 nucleotides. |
| `-D` | flag | off (Dcm reads dropped) | Keep reads with a Dcm site (CCWGG) around the motif. |
| `-h` | flag | — | Print help and exit. |

## The shared pipeline log

`clean_data.sh` writes `<prefix>.pipeline.log` into its output directory: the
run parameters, one section per step with its own parameters, read counts,
output file size and runtime, a summary table over all steps, and a run
metadata footer (script, git commit, full command line). Skipped steps are
listed as `skipped`, so the log always shows all five steps.

Passing that file to `analyse_flanks.sh` with `-g` appends the flank analysis
to it as *STEP 6*, with its own parameters, counts, runtime and metadata
footer — so one file documents the whole analysis:

```
==========================================================================
 Deep enzymology pipeline
==========================================================================
started:   2026-11-07 12:28:04
input:     AJ053_MKDL260006056-1A_23HJYLLT4_L5_1.fq.gz  (110M)
output:    WT_R1_pipeline

Parameters
  window (-c/-k)          7..135  (129 nt)
  quality cutoff (-q)     Phred 20
  quality rule            strict (-Q)
  adapter (-a)            AGATCGGAAGAGC
  forward barcode (-f)    TCC
  reverse barcode (-r)    GCT
  5 prime match           TCCGAGTGTGATTAGGTTTTTATTGTT
  3 prime match           GGTAGTGAGAGTTTAGTTATATTTGCT
  barcode trim (-b/-e)    3 / 3  ->  123 nt
  duplicate removal       on
  compressor              pigz -p 8

--------------------------------------------------------------------------
 STEP 1+2  hard trimming, adapter screening, quality filtering
--------------------------------------------------------------------------
  window kept:      7..135  (129 nt)
  quality cutoff:   Phred 20
  quality rule:     strict (-Q), any base below the cutoff discards the read
  adapter screened: AGATCGGAAGAGC

  Total reads processed:              2354253
  Reads shorter than the window:      0
  Reads with adapters (discarded):    2158
  Reads failing quality (discarded):  296038
  Reads written (passing filters):    2056057 (87.3%)

  output:    WT_R1_pipeline/AJ053_MKDL260006056-1A_23HJYLLT4_L5_1_q20strict_129bp.fq.gz  (64M)
  duration:  1m 18s

--------------------------------------------------------------------------
 STEP 3  read selection by barcodes
--------------------------------------------------------------------------
  forward barcode:  TCC
  reverse barcode:  GCT
  matched at the 5 prime end:  TCCGAGTGTGATTAGGTTTTTATTGTT
  matched at the 3 prime end:  GGTAGTGAGAGTTTAGTTATATTTGCT

  Total reads processed: 2056057
  Reads selected:        256642 (12.5%)

  output:    WT_R1_pipeline/AJ053_MKDL260006056-1A_23HJYLLT4_L5_1_q20strict_129bp_selected.fq.gz  (7.6M)
  duration:  7s

--------------------------------------------------------------------------
 STEP 4  barcode trimming
--------------------------------------------------------------------------
  trimmed from the 5 prime end:  3 nt
  trimmed from the 3 prime end:  3 nt
  resulting read length:         123 nt

  Total reads processed: 256642
  Reads written:         256642
  Reads too short:       0

  output:    WT_R1_pipeline/AJ053_MKDL260006056-1A_23HJYLLT4_L5_1_q20strict_selected_123bp.fq.gz  (7.5M)
  duration:  3s

--------------------------------------------------------------------------
 STEP 5  duplicate removal
--------------------------------------------------------------------------
  input:  WT_R1_pipeline/AJ053_MKDL260006056-1A_23HJYLLT4_L5_1_q20strict_selected_123bp.fq.gz

  256642 reads analyzed. 6.13 % duplicates found.
  Unique reads written:  240920

  output:    WT_R1_pipeline/AJ053_MKDL260006056-1A_23HJYLLT4_L5_1_q20strict_selected_123bp_unique.fq.gz  (7.2M)
  duration:  3s
  sequences: WT_R1_pipeline/AJ053_MKDL260006056-1A_23HJYLLT4_L5_1_q20strict_selected_123bp_unique_seq.txt.gz  (4.3M)

--------------------------------------------------------------------------
 SUMMARY  read cleaning
--------------------------------------------------------------------------
  step                                         reads in    reads out     kept
  STEP 1+2 trimming and quality filtering       2354253      2056057    87.3%
  STEP 3   read selection                       2056057       256642    12.5%
  STEP 4   barcode trimming                      256642       256642   100.0%
  STEP 5   duplicate removal                     256642       240920    93.9%
  ----------------------------------------------------------------------
  total                                         2354253       240920    10.2%

  final output:  WT_R1_pipeline/AJ053_MKDL260006056-1A_23HJYLLT4_L5_1_q20strict_selected_123bp_unique.fq.gz  (7.2M)
  total runtime: 1m 31s

--------------------------------------------------------------------------
script:    clean_data.sh
commit:    0.9.1
command:   clean_data.sh -i AJ053_MKDL260006056-1A_23HJYLLT4_L5_1.fq.gz -o WT_R1_pipeline -Q -f TCC -r GCT -p 8

--------------------------------------------------------------------------
 STEP 6  flanking sequence analysis
--------------------------------------------------------------------------
  input:            WT_R1_pipeline/AJ053_MKDL260006056-1A_23HJYLLT4_L5_1_q20strict_selected_123bp_unique_seq.txt.gz  (4.3M)
  motifs (-m)       CG,CA,CT,CC
  random parts      top at 25, bottom at 78, flank 10 nt
  analysed strand   top
  expected length   123 nt
  mismatch limit    0
  flank tables      1..3 nt
  Dcm filter        on (CCWGG reads discarded)

  Sequences read:                     240920
  Discarded, length is not 123 nt:    0
  Discarded, random parts do not fit: 0
  Discarded, motif does not match:    6283
  Discarded, mismatches above 0:      2374
  Discarded, Dcm site (CCWGG):        1287
  Discarded, unreadable motif state:  0
  Reads analysed:                     230976 (95.9%)
    motif CG                          60835
    motif CA                          54954
    motif CT                          76539
    motif CC                          38648

  report:    WT_R1_pipeline/WT_R1_flanks.txt  (276K)
  duration:  8s

--------------------------------------------------------------------------
script:    analyse_flanks.sh
commit:    1.0.0
command:   analyse_flanks.sh -i WT_R1_pipeline/AJ053_MKDL260006056-1A_23HJYLLT4_L5_1_q20strict_selected_123bp_unique_seq.txt.gz -o WT_R1_pipeline/WT_R1_flanks.txt -g WT_R1_pipeline/AJ053_MKDL260006056-1A_23HJYLLT4_L5_1.pipeline.log
```

---
