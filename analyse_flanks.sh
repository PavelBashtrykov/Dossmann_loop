#!/usr/bin/bash

RUN_TS=$(date +%y%m%d_%H%M)
SCRIPT_NAME=$(basename -- "${BASH_SOURCE[0]}")
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)
ARGS=( "$@" )

MOTIFS=CG,CA,CT,CC
POS1=25
POS3=78
FLANK=10
STRAND=top
READ_LEN=123
MISMATCH_LIMIT=0
MAX_FLANK=3
DCM=true
RUN_START=$(date +%s)

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

help() {
  printf '\nUsage: analyse_flanks.sh -i <STR> [-o <STR>] [-g <STR>] [-m <STR>] [-1 <INT>] [-3 <INT>] [-l <INT>] [-s <STR>] [-L <INT>] [-x <INT>] [-f <INT>] [-D] [-h]\n\n'
  echo '
  Flanking sequence analysis of a deep enzymology library. The input is the
  plain list of read sequences produced by clean_data.sh (one sequence per line,
  <...>_unique_seq.txt.gz, plain or gzipped).

  All defaults describe the library in use, a CN substrate of 123 nucleotides
  read on the top strand (random parts at 25 and 78, ten nucleotide flanks).
  Every default can still be overridden.

  Every read carries both strands of the substrate: the random part of the top
  strand starts at <-1>, the random part of the bottom strand at <-3>, both
  are <-l> + motif + <-l> nucleotides long. For every read the script

    1  looks up the central motif on both strands, accepting all bisulfite
       states of it (for CG: CG or TG on top, CG or TG on bottom); reads whose
       motif does not match any of <-m> are discarded
    2  reconstitutes the original, unconverted sequence of both random parts
       from the two strands, position by position; the central motif keeps its
       bisulfite state, so it still reports methylation
    3  discards reads with more than <-x> positions that no base pair explains
    4  discards reads with a Dcm site (CCWGG, read as CCAGG, CCTGG, CTAGG or
       CTTGG on the reconstituted top strand) around the motif, unless -D
    5  splits the remaining reads into methylated and unmethylated by the
       state of the motif on the analysed strand <-s> and counts

         - the frequency of every base at every position of the random part
         - the frequency of every flanking context, for flanks of 1 to <-f>
           nucleotides on each side

  The report is a tab separated table per motif and per flank length, ready to
  be pasted into a spreadsheet. Read counts go to stderr, not into it; with -g
  they are also appended to the shared pipeline log of clean_data.sh. Both the
  report and that log end with a metadata footer naming the script, its git
  commit and the command line that produced them.

  Options:

  -i <STR>                  /path/to/sequences.txt(.gz), one read per line,
                            LF or CRLF line endings
  -o <STR>                  /path/to/report.txt, default:
                            ./<prefix>_flanks_<strand>_YYMMDD_HHMM.txt
                            A value given here is used verbatim, no timestamp
                            is appended.
  -g <STR>                  /path/to/pipeline.log, the shared log written by
                            clean_data.sh. The read counts of this run are
                            appended to it as STEP 6, together with the
                            parameters used, so a single file documents the
                            whole analysis. The file is created when it does
                            not exist. Counts still go to stderr as well.
  -m <STR>                  Central motifs, comma separated, each starting
                            with C, default: CG,CA,CT,CC
                            A read is assigned to the first motif it matches.
  -1 <INT>                  First base of the top-strand random part, 1 based,
                            default: 25
  -3 <INT>                  First base of the bottom-strand random part,
                            1 based, default: 78
  -l <INT>                  Flank length of the random part on each side of
                            the motif, default: 10
  -s <STR>                  Strand to analyse, top or bottom, default: top
                            Methylation must be readable on it: CA, CT and CC
                            carry no cytosine on the bottom strand, and the
                            script aborts when a requested motif cannot report
                            methylation on the requested strand.
  -L <INT>                  Expected read length, reads of any other length
                            are discarded; 0 disables the check, default: 123
  -x <INT>                  Mismatch limit, positions of the random part that
                            no base pair explains, default: 0
  -f <INT>                  Largest flank reported in the flank tables, 1 to 5,
                            default: 3 (tables for 1, 2 and 3 nucleotides)
  -D                        Keep reads with a Dcm site (CCWGG) around the
                            motif, default: off, such reads are discarded
  -h                        Prints help info
'
}

OPTIND=1
while getopts ":i:o:g:m:1:3:l:s:L:x:f:Dh" option; do
  case "${option}" in
    i ) INPUT_FILE=${OPTARG} ;;
    o ) OUT_FILE=${OPTARG} ;;
    g ) LOG=${OPTARG} ;;
    m ) MOTIFS=${OPTARG} ;;
    1 ) POS1=${OPTARG} ;;
    3 ) POS3=${OPTARG} ;;
    l ) FLANK=${OPTARG} ;;
    s ) STRAND=${OPTARG} ;;
    L ) READ_LEN=${OPTARG} ;;
    x ) MISMATCH_LIMIT=${OPTARG} ;;
    f ) MAX_FLANK=${OPTARG} ;;
    D ) DCM=false ;;
    h ) help; exit ;;
    \? ) echo "Unknown option: -${OPTARG}" >&2; help; exit 1 ;;
    :  ) echo "Missing option argument for -${OPTARG}" >&2; help; exit 1 ;;
    *  ) echo "Unimplemented option: -${option}" >&2; help; exit 1 ;;
  esac
done

if ((OPTIND == 1))
then
  echo "No options specified. Terminating..." >&2
  help; exit 1
fi

if [ -z "${INPUT_FILE}" ]
then
  echo 'Input file with read sequences is missing. Use: -i <STR>' >&2
  help; exit 1
fi

if [ ! -f "${INPUT_FILE}" ]
then
  echo "Input file not found: ${INPUT_FILE}. Terminating..." >&2
  exit 1
fi

for PAIR in "POS1:-1" "POS3:-3" "FLANK:-l" "READ_LEN:-L" "MISMATCH_LIMIT:-x" "MAX_FLANK:-f"
do
  NAME=${PAIR%%:*}
  if [[ ! "${!NAME}" =~ ^[0-9]+$ ]]
  then
    echo "Option ${PAIR#*:} needs a non-negative integer, got: ${!NAME}. Terminating..." >&2
    exit 1
  fi
done

if (( POS1 < 1 )) || (( POS3 < 1 ))
then
  echo "Random part start positions are 1 based and must be at least 1. Terminating..." >&2
  exit 1
fi

if (( FLANK < 1 ))
then
  echo "Flank length (-l) must be at least 1. Terminating..." >&2
  exit 1
fi

if (( MAX_FLANK < 1 )) || (( MAX_FLANK > 5 )) || (( MAX_FLANK > FLANK ))
then
  echo "Largest reported flank (-f) must be between 1 and 5 and must not exceed -l (${FLANK}). Terminating..." >&2
  exit 1
fi

if [ "${STRAND}" != top ] && [ "${STRAND}" != bottom ]
then
  echo "Strand (-s) must be top or bottom, got: ${STRAND}. Terminating..." >&2
  exit 1
fi

MOTIFS=${MOTIFS^^}
if [[ ! "${MOTIFS}" =~ ^C[ACGT]*(,C[ACGT]*)*$ ]]
then
  echo "Motifs (-m) must be a comma separated list, each motif starting with C and containing A, C, G or T only, got: ${MOTIFS}. Terminating..." >&2
  exit 1
fi

PREFIX="$(basename -- "${INPUT_FILE}")"
PREFIX="${PREFIX%.gz}"
PREFIX="${PREFIX%.txt}"
PREFIX="${PREFIX%_unique_seq}"

OUT_FILE=${OUT_FILE:-$(pwd)/${PREFIX}_flanks_${STRAND}_${RUN_TS}.txt}

if [[ "${INPUT_FILE}" == *.gz ]]
then
  READER=(zcat)
else
  READER=(cat)
fi

STATS="${OUT_FILE}.counts.tmp"

LIB_TYPE="${MOTIFS}_pos1-${POS1}_pos3-${POS3}_flank-${FLANK}"

"${READER[@]}" "${INPUT_FILE}" \
| awk -v infile="${INPUT_FILE}" -v libtype="${LIB_TYPE}" -v motifs="${MOTIFS}" \
      -v pos1="${POS1}" -v pos3="${POS3}" -v fl="${FLANK}" -v strand="${STRAND}" \
      -v readlen="${READ_LEN}" -v sml="${MISMATCH_LIMIT}" -v maxf="${MAX_FLANK}" \
      -v dcm="${DCM}" '
  function complement(base) {
    if (base == "A") return "T"
    if (base == "T") return "A"
    if (base == "G") return "C"
    if (base == "C") return "G"
    return "N"
  }
  function revcomp(seq,   i, out) {
    out = ""
    for (i = length(seq); i >= 1; i--) out = out complement(substr(seq, i, 1))
    return out
  }
  # the bisulfite state of the motif as it reads on the strand that carries the
  # methylated cytosine: only the first C keeps its identity, every other C of
  # the motif is unmethylated and reads as T
  function top_state(motif, first,   i, c, out) {
    out = first
    for (i = 2; i <= length(motif); i++) {
      c = substr(motif, i, 1)
      out = out ((c == "C") ? "T" : c)
    }
    return out
  }
  BEGIN {
    nl = split("A,T,G,C", LETTER, ",")
    nm = split(motifs, MOTIF, ",")

    for (m = 1; m <= nm; m++) {
      mo = MOTIF[m]
      MLEN[m] = length(mo)
      RP[m] = 2 * fl + MLEN[m]

      TMETH[m] = top_state(mo, "C")
      TUNMETH[m] = top_state(mo, "T")

      # the opposite strand reads as the reverse complement of the motif, in
      # which every cytosine is either methylated (C) or converted (T)
      rc = revcomp(mo)
      nvar = 1; VAR[1] = ""
      for (i = 1; i <= length(rc); i++) {
        c = substr(rc, i, 1)
        if (c == "C") {
          for (v = 1; v <= nvar; v++) { VAR[nvar + v] = VAR[v] "T"; VAR[v] = VAR[v] "C" }
          nvar = nvar * 2
        } else {
          for (v = 1; v <= nvar; v++) VAR[v] = VAR[v] c
        }
      }
      for (v = 1; v <= nvar; v++) BOTTOM_OK[m, VAR[v]] = 1
      BMETH[m] = rc
      BUNMETH[m] = rc; gsub(/C/, "T", BUNMETH[m])

      if (strand == "top") { METH[m] = TMETH[m]; UNMETH[m] = TUNMETH[m] }
      else                 { METH[m] = BMETH[m]; UNMETH[m] = BUNMETH[m] }

      if (METH[m] == UNMETH[m]) {
        printf "Motif %s carries no cytosine on the %s strand, its methylation cannot be read there. Terminating...\n", mo, strand > "/dev/stderr"
        fatal = 1; exit 1
      }
      if (fl + MLEN[m] + maxf > RP[m]) {
        printf "Flank %d does not fit into the random part of motif %s. Terminating...\n", maxf, mo > "/dev/stderr"
        fatal = 1; exit 1
      }
    }
  }
  {
    sub(/\r$/, "")          # input may come from a DOS formatted file
    total++
    seqlen = length($0)
    if (readlen && seqlen != readlen) { wronglen++; next }

    matched = 0; fitted = 0
    for (m = 1; m <= nm; m++) {
      rp = RP[m]
      if (pos1 + rp - 1 > seqlen || pos3 + rp - 1 > seqlen) continue
      fitted++
      top = substr($0, pos1, rp)
      bot = substr($0, pos3, rp)
      tw = substr(top, fl + 1, MLEN[m])
      bw = substr(bot, fl + 1, MLEN[m])
      if ((tw == TMETH[m] || tw == TUNMETH[m]) && ((m, bw) in BOTTOM_OK)) { matched = m; break }
    }
    if (!fitted) { tooshort++; next }
    if (!matched) { nomotif++; next }

    # reconstitute both strands of the random part: the bottom strand read as
    # it is, the top strand read backwards, so position n of the bottom pairs
    # with position rp-n+1 of the top
    m = matched; rp = RP[m]; mism = 0
    for (n = 1; n <= rp; n++) {
      if (n >= fl + 1 && n <= fl + MLEN[m]) {
        TSEQ[n] = substr(top, n, 1)
        BSEQ[n] = substr(bot, n, 1)
        continue
      }
      b = substr(bot, n, 1)
      t = substr(top, rp - n + 1, 1)
      pair = b t
      if      (pair == "AT")                { tc = "T"; bc = "A" }
      else if (pair == "GT" || pair == "GC") { tc = "C"; bc = "G" }
      else if (pair == "TA")                { tc = "A"; bc = "T" }
      else if (pair == "TG" || pair == "CG") { tc = "G"; bc = "C" }
      else                                  { tc = "!"; bc = "!"; mism++ }
      TSEQ[rp - n + 1] = tc
      BSEQ[n] = bc
    }
    if (mism > sml) { mismatched++; next }

    tseq = ""; bseq = ""
    for (n = 1; n <= rp; n++) { tseq = tseq TSEQ[n]; bseq = bseq BSEQ[n] }

    # Dcm methylates the second C of CCWGG; that C is the motif cytosine, so
    # the site spans one base in front of the motif and three behind it
    if (dcm == "true" && fl >= 1 && fl + 4 <= rp) {
      if (substr(tseq, fl, 5) ~ /^C[CT][AT]GG$/) { dcmsite++; next }
    }

    seq = (strand == "top") ? tseq : bseq
    state = substr(seq, fl + 1, MLEN[m])
    if      (state == METH[m])   st = 0
    else if (state == UNMETH[m]) st = 1
    else { unreadable++; next }

    eligible++; ELIG[m]++

    canon = substr(seq, 1, fl) MOTIF[m] substr(seq, fl + MLEN[m] + 1)
    for (p = 1; p <= rp; p++) {
      c = substr(canon, p, 1)
      if (c == "A" || c == "T" || c == "G" || c == "C") NFREQ[m, st, c, p]++
      else nfskip++
    }

    for (f = 1; f <= maxf; f++) {
      left = substr(seq, fl - f + 1, f)
      right = substr(seq, fl + MLEN[m] + 1, f)
      if ((left right) ~ /^[ATGC]+$/) FLANKS[f, m, st, left, right]++
      else flskip++
    }
  }
  END {
    if (fatal) exit 1
    if (!total) {
      print "Input file holds no sequences. Terminating..." > "/dev/stderr"
      exit 1
    }

    printf "Filename: %s\n", infile
    printf "Library type: %s\n", libtype
    printf "Analysed strand: %s\n", strand
    printf "Number of Sequences found: %d\n", total

    sub_desc = ""
    for (m = 1; m <= nm; m++) sub_desc = sub_desc (m > 1 ? ", " : "") MOTIF[m]
    printf "\nCurrent substrate: [%s], pos1=%d, pos3=%d, flankLength=%d\n", sub_desc, pos1, pos3, fl

    # base frequencies per position, methylated and unmethylated separately
    for (m = 1; m <= nm; m++) {
      rp = RP[m]
      for (st = 0; st <= 1; st++) {
        printf "%s Reads Motif: %s\n", (st ? "Unmethylated" : "Methylated"), MOTIF[m]
        printf "Pos\t"
        for (p = 1; p <= rp; p++) printf "%d\t", p
        printf "\n"
        for (i = 1; i <= nl; i++) {
          c = LETTER[i]
          printf "%s", c
          for (p = 1; p <= rp; p++) printf "\t%d", NFREQ[m, st, c, p] + 0
          printf "\n"
        }
      }
    }

    # flank tables, one block per flank length, motifs side by side
    for (f = 1; f <= maxf; f++) {
      printf "\n"
      for (m = 1; m <= nm; m++) printf "flank\t5mC\tC\t"
      printf "\n"
      rows = 1
      for (i = 1; i <= 2 * f; i++) rows = rows * nl
      for (r = 0; r < rows; r++) {
        context = ""
        rest = r
        for (i = 2 * f; i >= 1; i--) {
          d = rest % nl
          rest = int(rest / nl)
          context = LETTER[d + 1] context
        }
        left = substr(context, 1, f)
        right = substr(context, f + 1)
        for (m = 1; m <= nm; m++) {
          printf "%s %s %s\t%d\t%d\t", left, MOTIF[m], right, \
                 FLANKS[f, m, 0, left, right] + 0, FLANKS[f, m, 1, left, right] + 0
        }
        printf "\n"
      }
    }

    printf "\nSequences read:                     %d\n", total > "/dev/stderr"
    if (readlen) printf "%-36s%d\n", sprintf("Discarded, length is not %d nt:", readlen), wronglen + 0 > "/dev/stderr"
    printf "%-36s%d\n", "Discarded, random parts do not fit:", tooshort + 0 > "/dev/stderr"
    printf "Discarded, motif does not match:    %d\n", nomotif + 0 > "/dev/stderr"
    printf "Discarded, mismatches above %d:      %d\n", sml, mismatched + 0 > "/dev/stderr"
    if (dcm == "true") printf "Discarded, Dcm site (CCWGG):        %d\n", dcmsite + 0 > "/dev/stderr"
    printf "Discarded, unreadable motif state:  %d\n", unreadable + 0 > "/dev/stderr"
    printf "Reads analysed:                     %d (%.1f%%)\n", eligible + 0, 100 * (eligible + 0) / total > "/dev/stderr"
    for (m = 1; m <= nm; m++) printf "  motif %-4s                        %d\n", MOTIF[m], ELIG[m] + 0 > "/dev/stderr"
  }' 2> "${STATS}" > "${OUT_FILE}"

STATUS=( "${PIPESTATUS[@]}" )
cat "${STATS}" >&2
if [ "${STATUS[0]}" -ne 0 ] || [ "${STATUS[1]}" -ne 0 ]
then
  echo "Flank analysis failed (reader ${STATUS[0]}, awk ${STATUS[1]}). Terminating..." >&2
  if [ -f "${OUT_FILE}" ]
  then
    rm -f "${OUT_FILE}"
  fi
  rm -f "${STATS}"
  exit 1
fi

COMMIT=$(git -C "${SCRIPT_DIR}" describe --always --dirty --tags 2> /dev/null)
COMMIT=${COMMIT:-unknown}

{
  printf -- '--------------------------------------------------------------------------\n'
  printf 'script:    %s\n' "${SCRIPT_NAME}"
  printf 'commit:    %s\n' "${COMMIT}"
  printf 'command:   %s %s\n' "${SCRIPT_NAME}" "${ARGS[*]}"
} >> "${OUT_FILE}"

if [ -n "${LOG}" ]
then
  mkdir -p -- "$(dirname -- "${LOG}")"
  {
    printf '\n'
    printf -- '--------------------------------------------------------------------------\n'
    printf ' STEP 6  flanking sequence analysis\n'
    printf -- '--------------------------------------------------------------------------\n'
    printf '  input:            %s  (%s)\n' "${INPUT_FILE}" "$(fmt_size "${INPUT_FILE}")"
    printf '  motifs (-m)       %s\n' "${MOTIFS}"
    printf '  random parts      top at %s, bottom at %s, flank %s nt\n' "${POS1}" "${POS3}" "${FLANK}"
    printf '  analysed strand   %s\n' "${STRAND}"
    printf '  expected length   %s\n' "$( ((READ_LEN)) && echo "${READ_LEN} nt" || echo 'not checked' )"
    printf '  mismatch limit    %s\n' "${MISMATCH_LIMIT}"
    printf '  flank tables      1..%s nt\n' "${MAX_FLANK}"
    printf '  Dcm filter        %s\n' "$([ "${DCM}" = true ] && echo 'on (CCWGG reads discarded)' || echo 'off (-D)')"
    awk '{ if (length($0)) print "  " $0; else print "" }' "${STATS}"
    printf '\n'
    printf '  report:    %s  (%s)\n' "${OUT_FILE}" "$(fmt_size "${OUT_FILE}")"
    printf '  duration:  %s\n' "$(hr_time $(( $(date +%s) - RUN_START )))"
    printf '\n'
    printf -- '--------------------------------------------------------------------------\n'
    printf 'script:    %s\n' "${SCRIPT_NAME}"
    printf 'commit:    %s\n' "${COMMIT}"
    printf 'command:   %s %s\n' "${SCRIPT_NAME}" "${ARGS[*]}"
  } >> "${LOG}"
  echo "Counts appended to: ${LOG}" >&2
fi

rm -f "${STATS}"

echo "Report written to: ${OUT_FILE}" >&2
