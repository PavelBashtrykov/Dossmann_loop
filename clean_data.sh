#!/usr/bin/bash

MAIN_DIR=$(pwd)
CORES=2
RUN_TS=$(date +%y%m%d_%H%M)

DEDUP=true
ADAPTER=AGATCGGAAGAGC

FIRST_KEEP=7
LAST_KEEP=135
QUALITY=20
STRICT=false
ANCHOR_5=GAGTGTGATTAGGTTTTTATTGTT
ANCHOR_3=GGTAGTGAGAGTTTAGTTATATTT
TRIM_5=3
TRIM_3=3

SCRIPT_NAME=$(basename -- "${BASH_SOURCE[0]}")
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)
ARGS=( "$@" )
RUN_START=$(date +%s)

STEP_LABEL=()
STEP_IN=()
STEP_OUT=()

hr_time(){
  local S=$1
  if (( S < 60 ))
  then
    printf '%ds' "${S}"
  elif (( S < 3600 ))
  then
    printf '%dm %02ds' $(( S / 60 )) $(( S % 60 ))
  else
    printf '%dh %02dm %02ds' $(( S / 3600 )) $(( (S % 3600) / 60 )) $(( S % 60 ))
  fi
}

fmt_size(){
  if [ -f "$1" ]
  then
    du -h -- "$1" | cut -f1
  else
    echo 'n/a'
  fi
}

log(){
  printf '%s\n' "$*" >> "${LOG}"
}

log_rule(){
  log '--------------------------------------------------------------------------'
}

log_step(){
  log ''
  log_rule
  log " $1  $2"
  log_rule
}

log_stats(){
  grep -v '^#KV' -- "$1" \
  | awk '{ if (length($0)) print "  " $0; else print "" }' >> "${LOG}"
}

kv_get(){
  awk -v k="$2" '$1 == "#KV" { for (i = 2; i <= NF; i++) { split($i, a, "="); if (a[1] == k) print a[2] } }' "$1"
}

log_result(){
  # log_result <output file> <step start epoch>
  log ''
  log "  output:    $1  ($(fmt_size "$1"))"
  log "  duration:  $(hr_time $(( $(date +%s) - $2 )))"
}

record_step(){
  STEP_LABEL+=( "$1" )
  STEP_IN+=( "$2" )
  STEP_OUT+=( "$3" )
}

record_skipped(){
  STEP_LABEL+=( "$1" )
  STEP_IN+=( '-' )
  STEP_OUT+=( '-' )
}

select_reads(){
  local T0
  T0=$(date +%s)
  log_step 'STEP 3' 'read selection by barcodes'
  log "  forward barcode:  ${FORWARD:-(none)}"
  log "  reverse barcode:  ${REVERSE:-(none)}"
  log "  matched at the 5 prime end:  ${MATCH_5:-(nothing)}"
  log "  matched at the 3 prime end:  ${MATCH_3:-(nothing)}"
  log ''

  if [ ! -f "${QUAL_FILTERED}" ]
  then
    echo "Missing quality-filtered file: ${QUAL_FILTERED}. Terminating..." >&2
    exit 1
  fi

  local STATS STATUS
  STATS="${OUT_DIR}/.${PREFIX}_select_stats"

  zcat "${QUAL_FILTERED}" \
  | awk -v b="${MATCH_5}" -v e="${MATCH_3}" '
      BEGIN { nb = length(b); ne = length(e) }
      NR % 4 == 1 { h = $0; next }
      NR % 4 == 2 { s = $0; next }
      NR % 4 == 3 { p = $0; next }
      {
        total++
        if (nb && substr(s, 1, nb) != b) next
        if (ne && substr(s, length(s) - ne + 1) != e) next
        kept++
        print h "\n" s "\n" p "\n" $0
      }
      END {
        printf "Total reads processed: %d\n", total > "/dev/stderr"
        printf "Reads selected:        %d (%.1f%%)\n", kept, (total ? 100 * kept / total : 0) > "/dev/stderr"
        printf "#KV in=%d out=%d\n", total, kept > "/dev/stderr"
      }' 2> "${STATS}" \
  | "${COMPRESS[@]}" > "${SELECTED}"

  STATUS=( "${PIPESTATUS[@]}" )
  log_stats "${STATS}"

  if [ "${STATUS[0]}" -ne 0 ] || [ "${STATUS[1]}" -ne 0 ] || [ "${STATUS[2]}" -ne 0 ]
  then
    echo "Read selection failed (zcat ${STATUS[0]}, awk ${STATUS[1]}, compressor ${STATUS[2]}). Terminating..." >&2
    rm -f "${SELECTED}" "${STATS}"
    exit 1
  fi

  record_step 'STEP 3   read selection' "$(kv_get "${STATS}" in)" "$(kv_get "${STATS}" out)"
  rm -f "${STATS}"
  log_result "${SELECTED}" "${T0}"
}

trim_barcodes(){
  local T0
  T0=$(date +%s)
  log_step 'STEP 4' 'barcode trimming'
  log "  trimmed from the 5 prime end:  ${TRIM_5} nt"
  log "  trimmed from the 3 prime end:  ${TRIM_3} nt"
  log "  resulting read length:         ${TRIMMED_LEN} nt"
  log ''

  if [ ! -f "${SELECTED}" ]
  then
    echo "Missing selected file: ${SELECTED}. Terminating..." >&2
    exit 1
  fi

  local STATS STATUS
  STATS="${OUT_DIR}/.${PREFIX}_trim_stats"

  zcat "${SELECTED}" \
  | awk -v t5="${TRIM_5}" -v t3="${TRIM_3}" '
      NR % 4 == 1 { h = $0; next }
      NR % 4 == 2 { s = $0; next }
      NR % 4 == 3 { p = $0; next }
      {
        total++
        keep = length(s) - t5 - t3
        if (keep < 1) { dropped++; next }
        kept++
        print h "\n" substr(s, t5 + 1, keep) "\n" p "\n" substr($0, t5 + 1, keep)
      }
      END {
        printf "Total reads processed: %d\n", total > "/dev/stderr"
        printf "Reads written:         %d\n", kept > "/dev/stderr"
        printf "Reads too short:       %d\n", dropped > "/dev/stderr"
        printf "#KV in=%d out=%d\n", total, kept > "/dev/stderr"
      }' 2> "${STATS}" \
  | "${COMPRESS[@]}" > "${BARCODE_TRIMMED}"

  STATUS=( "${PIPESTATUS[@]}" )
  log_stats "${STATS}"

  if [ "${STATUS[0]}" -ne 0 ] || [ "${STATUS[1]}" -ne 0 ] || [ "${STATUS[2]}" -ne 0 ]
  then
    echo "Barcode trimming failed (zcat ${STATUS[0]}, awk ${STATUS[1]}, compressor ${STATUS[2]}). Terminating..." >&2
    rm -f "${BARCODE_TRIMMED}" "${STATS}"
    exit 1
  fi

  record_step 'STEP 4   barcode trimming' "$(kv_get "${STATS}" in)" "$(kv_get "${STATS}" out)"
  rm -f "${STATS}"
  log_result "${BARCODE_TRIMMED}" "${T0}"
}

remove_duplicates(){
  local SRC STATS STATUS T0
  SRC="$1"
  T0=$(date +%s)
  log_step 'STEP 5' 'duplicate removal'
  log "  input:  ${SRC}"
  log ''

  if [ ! -f "${SRC}" ]
  then
    echo "Missing input for duplicate removal: ${SRC}. Terminating..." >&2
    exit 1
  fi

  DEDUPED="${SRC%.fq.gz}_unique.fq.gz"
  SEQ_LIST="${SRC%.fq.gz}_unique_seq.txt.gz"
  STATS="${OUT_DIR}/.${PREFIX}_dedup_stats"

  zcat "${SRC}" \
  | awk -v seqlist="${SEQ_LIST}" '
      BEGIN { gz = "gzip > \"" seqlist "\"" }
      NR % 4 == 1 { h = $0; next }
      NR % 4 == 2 { s = $0; next }
      NR % 4 == 3 { p = $0; next }
      {
        total++
        if (s in seen) { dups++; next }
        seen[s] = 1
        kept++
        print h "\n" s "\n" p "\n" $0
        print s | gz
      }
      END {
        if (close(gz) != 0)
        {
          printf "Writing %s failed.\n", seqlist > "/dev/stderr"
          exit 1
        }
        printf "%d reads analyzed. %.2f %% duplicates found.\n", \
               total, (total ? 100 * dups / total : 0) > "/dev/stderr"
        printf "Unique reads written:  %d\n", kept > "/dev/stderr"
        printf "#KV in=%d out=%d\n", total, kept > "/dev/stderr"
      }' 2> "${STATS}" \
  | "${COMPRESS[@]}" > "${DEDUPED}"

  STATUS=( "${PIPESTATUS[@]}" )
  log_stats "${STATS}"

  if [ "${STATUS[0]}" -ne 0 ] || [ "${STATUS[1]}" -ne 0 ] || [ "${STATUS[2]}" -ne 0 ]
  then
    echo "Duplicate removal failed (zcat ${STATUS[0]}, awk ${STATUS[1]}, compressor ${STATUS[2]}). Terminating..." >&2
    rm -f "${DEDUPED}" "${SEQ_LIST}" "${STATS}"
    exit 1
  fi

  record_step 'STEP 5   duplicate removal' "$(kv_get "${STATS}" in)" "$(kv_get "${STATS}" out)"
  rm -f "${STATS}"
  log_result "${DEDUPED}" "${T0}"
  log "  sequences: ${SEQ_LIST}  ($(fmt_size "${SEQ_LIST}"))"
}

trim_and_filter(){
  # Steps 1 and 2: window extraction, adapter screening and quality trimming
  # in a single pass. Writes ${QUAL_FILTERED} directly.
  local STATS STATUS T0
  local -a READER
  T0=$(date +%s)
  STATS="${OUT_DIR}/.${PREFIX}_trim_stats_s12"

  log_step 'STEP 1+2' 'hard trimming, adapter screening, quality filtering'
  log "  window kept:      ${FIRST_KEEP}..${LAST_KEEP}  (${READ_LEN} nt)"
  log "  quality cutoff:   Phred ${QUALITY}"
  log "  quality rule:     $([ "${STRICT}" = true ] \
                              && echo 'strict (-Q), any base below the cutoff discards the read' \
                              || echo 'running sum of cutadapt/BWA')"
  log "  adapter screened: ${ADAPTER}"
  log ''

  if [[ "${INPUT_FILE}" == *.gz ]]
  then
    READER=(zcat)
  else
    READER=(cat)
  fi

  "${READER[@]}" "${INPUT_FILE}" \
  | awk -v first="${FIRST_KEEP}" -v win="${READ_LEN}" \
        -v cutoff="${QUALITY}" -v adapter="${ADAPTER}" \
        -v strict="${STRICT}" '
      # adapter aligned at start, at most one substituted base
      function fits(seq, start, n,   k, d) {
        if (start < 1 || start + n - 1 > length(seq)) return 0
        d = 0
        for (k = 1; k <= n; k++)
          if (substr(seq, start+k-1, 1) != substr(ad, k, 1)) if (++d > 1) return 0
        return 1
      }
      # adapter aligned at start, one base missing from or inserted into the read
      function fits_indel(seq, start, n,   k, d) {
        if (start < 1 || start + n - 2 > length(seq)) return 0
        d = 0
        for (k = 1; k <= n-1; k++) if (substr(seq, start+k-1, 1) != substr(ad, k+1, 1)) { d = 1; break }
        if (!d) return 1
        d = 0
        for (k = 1; k <= n-1; k++) if (substr(seq, start+k, 1) != substr(ad, k, 1)) { d = 1; break }
        return !d
      }
      function check(seq, start) {
        if (start < 1 || start + alen - 1 > length(seq)) return 0
        return fits(seq, start, alen) || fits_indel(seq, start, alen)
      }
      function has_adapter(seq,   pos, j, t, o, a) {
        if (index(seq, probe)) return 1
        # one error cannot hit both halves of the adapter, so an exact match of
        # either half brackets every candidate position (pigeonhole) and the
        # full comparison runs only there; the offsets cover indel shifts
        for (t = 1; t <= 2; t++) {
          o = (t == 1) ? 0 : h1len
          a = (t == 1) ? ah1 : ah2
          pos = 0
          while ((j = index(substr(seq, pos+1), a)) > 0) {
            pos += j
            if (check(seq, pos-o) || check(seq, pos-o-1) || check(seq, pos-o+1)) return 1
          }
        }
        return 0
      }
      BEGIN {
        for (i = 33; i <= 126; i++) Q[sprintf("%c", i)] = i - 33
        ad = adapter; alen = length(ad)
        probe = substr(ad, 1, alen - 1)
        h1len = int(alen / 2); ah1 = substr(ad, 1, h1len); ah2 = substr(ad, h1len + 1)
      }
      NR % 4 == 1 { h = $0; next }
      NR % 4 == 2 { s = $0; next }
      NR % 4 == 3 { p = $0; next }
      {
        total++
        if (length(s) < first + win - 1) { tooshort++; next }
        seq  = substr(s,  first, win)
        qual = substr($0, first, win)
        if (alen > 1 && has_adapter(seq)) { withadapter++; next }
        if (strict == "true")
        {
          # strict rule: a single base below the cutoff discards the whole read
          bad = 0
          for (i = 1; i <= win; i++)
            if (Q[substr(qual, i, 1)] < cutoff) { bad = 1; break }
          if (bad) { lowqual++; next }
        }
        else
        {
          # quality trimming, same running-sum rule as cutadapt/BWA use
          sum = 0; maxsum = 0; cut = win
          for (i = win; i >= 1; i--)
          {
            sum += cutoff - Q[substr(qual, i, 1)]
            if (sum < 0) break
            if (sum > maxsum) { maxsum = sum; cut = i - 1 }
          }
          if (cut < win) { lowqual++; next }
        }
        kept++
        print h "\n" seq "\n" p "\n" qual
      }
      END {
        printf "Total reads processed:              %d\n", total > "/dev/stderr"
        printf "Reads shorter than the window:      %d\n", tooshort > "/dev/stderr"
        printf "Reads with adapters (discarded):    %d\n", withadapter > "/dev/stderr"
        printf "Reads failing quality (discarded):  %d\n", lowqual > "/dev/stderr"
        printf "Reads written (passing filters):    %d (%.1f%%)\n", \
               kept, (total ? 100 * kept / total : 0) > "/dev/stderr"
        printf "#KV in=%d out=%d\n", total, kept > "/dev/stderr"
      }' 2> "${STATS}" \
  | "${COMPRESS[@]}" > "${QUAL_FILTERED}"

  STATUS=( "${PIPESTATUS[@]}" )
  log_stats "${STATS}"

  if [ "${STATUS[0]}" -ne 0 ] || [ "${STATUS[1]}" -ne 0 ] || [ "${STATUS[2]}" -ne 0 ]
  then
    echo "Trimming and filtering failed (reader ${STATUS[0]}, awk ${STATUS[1]}, compressor ${STATUS[2]}). Terminating..." >&2
    rm -f "${QUAL_FILTERED}" "${STATS}"
    exit 1
  fi

  record_step 'STEP 1+2 trimming and quality filtering' "$(kv_get "${STATS}" in)" "$(kv_get "${STATS}" out)"
  rm -f "${STATS}"
  log_result "${QUAL_FILTERED}" "${T0}"
}

help() {
  printf '\nUsage: clean_data.sh -i <STR> [-o <STR>] [-p <INT>] [-q <INT>] [-Q] [-c <INT>] [-k <INT>] [-f <STR>] [-r <STR>] [-b <INT>] [-e <INT>] [-a <STR>] [-n] [-h]\n\n'
  echo '
  Single-end pipeline for deep enzymology libraries. Steps, in order:

    STEP 1  hard trimming      keep bases <-c>..<-k> of every read; reads
                               shorter than the window are discarded
    STEP 2  quality filtering  Phred <-q> cutoff and adapter screening
    STEP 3  read selection     keep only reads starting with <-f> plus the 5
                               prime constant part of the substrate and ending
                               with the 3 prime constant part plus <-r>;
                               skipped when neither barcode is given
    STEP 4  barcode trimming   cut <-b> bases off the 5 prime and <-e> off the
                               3 prime end; runs together with STEP 3
    STEP 5  duplicate removal  collapse reads sharing an identical full-length
                               sequence; skipped with -n

  Steps 1 and 2 run as a single pass over the input and need nothing but awk,
  zcat and gzip (pigz when present). Quality trimming uses the running-sum rule
  of cutadapt/BWA: a read is discarded when that rule would shorten it below the
  full window. With -Q the strict rule is used instead: a single base below the
  cutoff anywhere in the window discards the whole read.

  Every step writes its own .fq.gz into the output directory. All parameters,
  per step counts, file sizes, runtimes and a summary table go into the shared
  pipeline log <prefix>.pipeline.log. Pass that log to analyse_flanks.sh with
  its -g option and the flank analysis is appended to the same file, so one
  log covers the whole analysis.

  Options:

  -i <STR>                  /path/to/reads.fq.gz (single-end FASTQ, plain or gzipped)
  -o <STR>                  /path/to/output, default: ./<prefix>_YYMMDD_HHMM
  -p <INT>                  e.g. 4, threads used for compression (pigz only),
                            default: 2
  -q <INT>                  Quality trimming cutoff (Phred), default: 20
  -Q                        Strict quality filtering in step 2: discard the read
                            when any base of the kept window is below the -q
                            cutoff, instead of the running-sum rule. Reads are
                            never shortened either way. The output files of
                            steps 2 to 4 carry q<Q>strict instead of q<Q>, so a
                            strict run never overwrites a running-sum run in
                            the same directory.
  -c <INT>                  First base to keep, default: 7
  -k <INT>                  Last base to keep, default: 135
                            Defaults keep positions 7..135 = 129 nt
  -b <INT>                  Bases trimmed from the 5 prime end in step 4, default: 3
  -e <INT>                  Bases trimmed from the 3 prime end in step 4, default: 3
  -f <STR>                  Forward barcode, the read must start with it
  -r <STR>                  Reverse barcode, the read must end with it
                            Both are matched exactly, as given, on the read as it is
                            after step 2. Reverse-complement the reverse barcode
                            yourself if that is what your library requires.
                            STEP 3 is skipped when neither is given.
                            Each barcode is joined with the constant part of the
                            substrate next to it, so the read is anchored on the
                            barcode and 24 more nucleotides on that side:
                              5 prime  <-f> + GAGTGTGATTAGGTTTTTATTGTT
                              3 prime  GGTAGTGAGAGTTTAGTTATATTT + <-r>
                            The constant parts are held in ANCHOR_5 and ANCHOR_3
                            at the top of the script; edit them there when the
                            substrate changes.
  -a <STR>                  Adapter screened for in step 2, default:
                            AGATCGGAAGAGC (Illumina). A read is discarded when
                            the adapter occurs anywhere in it, allowing one
                            substituted, missing or inserted base.
  -n                        Skip STEP 5, do not remove PCR duplicates
                            STEP 5 keeps the first record of every read sequence
                            seen (exact, full-length match) and reports the
                            duplicate fraction in the log. Outputs:
                            <last step output>_unique.fq.gz       (FASTQ)
                            <last step output>_unique_seq.txt.gz  (one unique
                            read sequence per line, same order as the FASTQ)
  -h                        Prints help info
'
}

OPTIND=1
while getopts ":i:o:p:q:Qc:k:b:e:f:r:a:nh" option; do
  case "${option}" in
    i ) INPUT_FILE=${OPTARG} ;;
    n ) DEDUP=false ;;
    Q ) STRICT=true ;;
    a ) ADAPTER=${OPTARG} ;;
    f ) FORWARD=${OPTARG} ;;
    r ) REVERSE=${OPTARG} ;;
    b ) TRIM_5=${OPTARG} ;;
    e ) TRIM_3=${OPTARG} ;;
    o ) OUT_DIR=${OPTARG} ;;
    p ) CORES=${OPTARG} ;;
    q ) QUALITY=${OPTARG} ;;
    c ) FIRST_KEEP=${OPTARG} ;;
    k ) LAST_KEEP=${OPTARG} ;;
    h ) help; exit;;
    \? ) echo "Unknown option: -${OPTARG}" >&2; help; exit 1;;
    :  ) echo "Missing option argument for -${OPTARG}" >&2; help; exit 1;;
    *  ) echo "Unimplemented option: -$option" >&2; help; exit 1;;
  esac
done

if ((OPTIND == 1))
then
    echo "No options specified. Terminating..." >&2
    help; exit 1
fi

if [ -z "${INPUT_FILE}" ]
then
  echo 'Input FASTQ file is missing. Use: -i <STR>' >&2
  help
  exit 1
fi

if [ ! -f "${INPUT_FILE}" ]
then
  echo "Input file not found: ${INPUT_FILE}. Terminating..." >&2
  exit 1
fi

READ_LEN=$(( LAST_KEEP - FIRST_KEEP + 1 ))
if (( READ_LEN < 1 ))
then
  echo "Last base to keep (${LAST_KEEP}) must not be smaller than the first one (${FIRST_KEEP}). Terminating..." >&2
  exit 1
fi

PREFIX="$(basename -- "${INPUT_FILE}")"
PREFIX="${PREFIX%.gz}"
PREFIX="${PREFIX%.fq}"
PREFIX="${PREFIX%.fastq}"

if [ -n "${FORWARD}${REVERSE}" ] && [[ "${FORWARD}${REVERSE}" =~ [^ACGTNacgtn] ]]
then
  echo "Barcodes may only contain A, C, G, T or N. Terminating..." >&2
  exit 1
fi

# what STEP 3 matches: the barcode together with the constant part next to it
MATCH_5="${FORWARD:+${FORWARD}${ANCHOR_5}}"
MATCH_3="${REVERSE:+${ANCHOR_3}${REVERSE}}"

OUT_DIR=${OUT_DIR:-${MAIN_DIR}/${PREFIX}_${RUN_TS}}
LOG="${OUT_DIR}/${PREFIX}.pipeline.log"

if command -v pigz > /dev/null 2>&1
then
  COMPRESS=(pigz -p "${CORES}")
else
  COMPRESS=(gzip)
fi
QTAG="q${QUALITY}"
[ "${STRICT}" = true ] && QTAG="${QTAG}strict"
QUAL_FILTERED="${OUT_DIR}/${PREFIX}_${QTAG}_${READ_LEN}bp.fq.gz"
SELECTED="${OUT_DIR}/${PREFIX}_${QTAG}_${READ_LEN}bp_selected.fq.gz"

TRIMMED_LEN=$(( READ_LEN - TRIM_5 - TRIM_3 ))
if (( TRIMMED_LEN < 1 ))
then
  echo "Trimming ${TRIM_5} + ${TRIM_3} bases leaves nothing of a ${READ_LEN} nt read. Terminating..." >&2
  exit 1
fi
BARCODE_TRIMMED="${OUT_DIR}/${PREFIX}_${QTAG}_selected_${TRIMMED_LEN}bp.fq.gz"

if (( ${#MATCH_5} + ${#MATCH_3} > READ_LEN ))
then
  echo "Barcodes plus constant parts need ${#MATCH_5} + ${#MATCH_3} nt, more than the ${READ_LEN} nt window. Terminating..." >&2
  exit 1
fi

mkdir -p "${OUT_DIR}"
: > "${LOG}"

log '=========================================================================='
log ' Deep enzymology pipeline'
log '=========================================================================='
log "started:   $(date '+%Y-%m-%d %H:%M:%S')"
log "input:     ${INPUT_FILE}  ($(fmt_size "${INPUT_FILE}"))"
log "output:    ${OUT_DIR}"
log ''
log 'Parameters'
log "  window (-c/-k)          ${FIRST_KEEP}..${LAST_KEEP}  (${READ_LEN} nt)"
log "  quality cutoff (-q)     Phred ${QUALITY}"
log "  quality rule            $([ "${STRICT}" = true ] \
                                  && echo 'strict (-Q)' \
                                  || echo 'running sum (cutadapt/BWA)')"
log "  adapter (-a)            ${ADAPTER}"
log "  forward barcode (-f)    ${FORWARD:-(none)}"
log "  reverse barcode (-r)    ${REVERSE:-(none)}"
log "  5 prime match           ${MATCH_5:-(nothing)}"
log "  3 prime match           ${MATCH_3:-(nothing)}"
log "  barcode trim (-b/-e)    ${TRIM_5} / ${TRIM_3}  ->  ${TRIMMED_LEN} nt"
log "  duplicate removal       $([ "${DEDUP}" = true ] && echo 'on' || echo 'off (-n)')"
log "  compressor              ${COMPRESS[*]}"

trim_and_filter

LAST_OUT="${QUAL_FILTERED}"

if [ -n "${FORWARD}" ] || [ -n "${REVERSE}" ]
then
  select_reads
  trim_barcodes
  LAST_OUT="${BARCODE_TRIMMED}"
else
  record_skipped 'STEP 3   read selection'
  record_skipped 'STEP 4   barcode trimming'
fi

if [ "${DEDUP}" = true ]
then
  remove_duplicates "${LAST_OUT}"
  LAST_OUT="${DEDUPED}"
else
  record_skipped 'STEP 5   duplicate removal'
fi

log ''
log_rule
log ' SUMMARY  read cleaning'
log_rule
printf '  %-40s %12s %12s %8s\n' 'step' 'reads in' 'reads out' 'kept' >> "${LOG}"
for I in "${!STEP_LABEL[@]}"
do
  if [ "${STEP_IN[I]}" = '-' ]
  then
    printf '  %-40s %12s %12s %8s\n' "${STEP_LABEL[I]}" '-' '-' 'skipped' >> "${LOG}"
  else
    printf '  %-40s %12d %12d %7.1f%%\n' "${STEP_LABEL[I]}" "${STEP_IN[I]}" "${STEP_OUT[I]}" \
           "$(awk -v a="${STEP_OUT[I]}" -v b="${STEP_IN[I]}" 'BEGIN { print (b ? 100 * a / b : 0) }')" >> "${LOG}"
  fi
done

TOTAL_IN=''
TOTAL_OUT=''
for I in "${!STEP_LABEL[@]}"
do
  if [ "${STEP_IN[I]}" != '-' ]
  then
    TOTAL_IN=${TOTAL_IN:-${STEP_IN[I]}}
    TOTAL_OUT=${STEP_OUT[I]}
  fi
done

log '  ----------------------------------------------------------------------'
printf '  %-40s %12d %12d %7.1f%%\n' 'total' "${TOTAL_IN}" "${TOTAL_OUT}" \
       "$(awk -v a="${TOTAL_OUT}" -v b="${TOTAL_IN}" 'BEGIN { print (b ? 100 * a / b : 0) }')" >> "${LOG}"
log ''
log "  final output:  ${LAST_OUT}  ($(fmt_size "${LAST_OUT}"))"
log "  total runtime: $(hr_time $(( $(date +%s) - RUN_START )))"

COMMIT=$(git -C "${SCRIPT_DIR}" describe --always --dirty --tags 2> /dev/null)
COMMIT=${COMMIT:-unknown}
log ''
log_rule
log "script:    ${SCRIPT_NAME}"
log "commit:    ${COMMIT}"
log "command:   ${SCRIPT_NAME} ${ARGS[*]}"

echo "Log written to: ${LOG}" >&2
