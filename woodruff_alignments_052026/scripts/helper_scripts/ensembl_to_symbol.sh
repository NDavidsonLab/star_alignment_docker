#!/usr/bin/env bash
# =============================================================================
# ensembl_to_symbol.sh
#
# Takes a counts TSV (first col = ENSMUSG IDs with version, rest = samples)
# and prepends a "symbol" column. No jq required — pure curl + awk.
#
# Usage:
#   bash ensembl_to_symbol.sh counts.tsv > counts_with_symbols.tsv
# =============================================================================

set -euo pipefail

ENSEMBL_REST="https://rest.ensembl.org"
BATCH_SIZE=1000
RETRY=3
RETRY_DELAY=5
INPUT="${1:--}"

# ---------------------------------------------------------------------------
# 1. Extract unique Ensembl IDs (strip version suffix)
# ---------------------------------------------------------------------------
mapfile -t IDS < <(
  awk 'NR>1 { split($1, a, "."); print a[1] }' "$INPUT" | sort -u
)
TOTAL=${#IDS[@]}
echo "# Translating $TOTAL Ensembl ID(s) ..." >&2

# ---------------------------------------------------------------------------
# 2. Batch lookup via Ensembl REST, parse JSON with awk (POSIX-safe)
# ---------------------------------------------------------------------------
LOOKUP_FILE=$(mktemp)
trap 'rm -f "$LOOKUP_FILE"' EXIT

ensembl_lookup() {
  local ids_json="$1"
  for attempt in $(seq 1 $RETRY); do
    response=$(curl -sf --retry 2 \
      -H "Content-Type: application/json" \
      -H "Accept: application/json" \
      -d "{\"ids\": ${ids_json}}" \
      "${ENSEMBL_REST}/lookup/id") && echo "$response" && return
    echo "# Attempt $attempt/$RETRY failed, retrying in ${RETRY_DELAY}s ..." >&2
    sleep "$RETRY_DELAY"
  done
  echo "{}"
}

# POSIX awk JSON parser — no 3-arg match, no gensub
# Splits on display_name": " to find symbols, then pairs with nearest ENSMUSG id
parse_json_awk() {
  awk '
  {
    line = $0
    # Replace all { and , with newlines to process one field at a time
    gsub(/[{},]/, "\n", line)
    n = split(line, fields, "\n")

    current_id = ""
    for (i = 1; i <= n; i++) {
      f = fields[i]

      # Match an Ensembl gene ID key: "ENSMUSG00000051951"
      if (f ~ /"ENSMUSG[0-9]+"/) {
        # Extract just the ID between quotes
        id = f
        sub(/.*"ENSMUSG/, "ENSMUSG", id)
        sub(/".*/, "", id)
        current_id = id
      }

      # Match display_name value: "display_name":"Xkr4"
      if (f ~ /"display_name"/ && current_id != "") {
        sym = f
        sub(/.*"display_name":"/, "", sym)
        sub(/".*/, "", sym)
        print current_id "\t" sym
        current_id = ""
      }
    }
  }'
}

for (( i=0; i<TOTAL; i+=BATCH_SIZE )); do
  batch=("${IDS[@]:$i:$BATCH_SIZE}")
  echo "# Batch $((i/BATCH_SIZE + 1)): ${#batch[@]} IDs ..." >&2

  ids_json=$(printf '"%s",' "${batch[@]}" | sed 's/,$//' | sed 's/^/[/' | sed 's/$/]/')

  ensembl_lookup "$ids_json" | parse_json_awk >> "$LOOKUP_FILE"

  sleep 0.3
done

echo "# Lookup complete. Joining symbols ..." >&2

# ---------------------------------------------------------------------------
# 3. Join symbol back into original TSV with awk
# ---------------------------------------------------------------------------
awk -v lookup="$LOOKUP_FILE" '
BEGIN {
  FS = OFS = "\t"
  while ((getline line < lookup) > 0) {
    split(line, a, "\t")
    sym[a[1]] = a[2]
  }
}
NR == 1 {
  print "symbol", $0
  next
}
{
  split($1, a, ".")
  clean_id = a[1]
  symbol = (clean_id in sym) ? sym[clean_id] : "NA"
  print symbol, $0
}
' "$INPUT"