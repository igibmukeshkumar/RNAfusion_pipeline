#!/usr/bin/env bash
set -euo pipefail

studies=()
out=""
pattern=""

usage() {
    echo "Usage: $0 -i study_path [-i study_path ...] -o output.tsv -p 'pattern'" >&2
}

while getopts "i:o:p:" opt; do
    case "$opt" in
        i) studies+=("$OPTARG") ;;
        o) out="$OPTARG" ;;
        p) pattern="$OPTARG" ;;
        *) usage; exit 2 ;;
    esac
done

if [[ ${#studies[@]} -eq 0 || -z "$out" || -z "$pattern" ]]; then
    usage
    exit 2
fi

for study in "${studies[@]}"; do
    if [[ ! -d "$study" ]]; then
        echo "ERROR: Study directory not found: $study" >&2
        exit 1
    fi
done

out_dir=$(dirname -- "$out")
if [[ ! -d "$out_dir" ]]; then
    echo "ERROR: Output directory not found: $out_dir" >&2
    exit 1
fi

tmp=$(mktemp "${out}.tmp.XXXXXX")
trap 'rm -f -- "$tmp"' EXIT

first=1
count=0

for study in "${studies[@]}"; do
    study=$(realpath -- "$study")

    while IFS= read -r -d '' f; do
        count=$((count + 1))
        f=$(realpath -- "$f")
        rel="${f#"$study"/}"
        sample="${rel%%/*}"

        echo "Adding: $sample" >&2

        if [[ "$first" -eq 1 ]]; then
            awk -v p="$f" -v s="$sample" \
                'BEGIN{FS="\t"; OFS="\t"}
                 NR==1 {print $0,"path","sample_name"; next}
                 {print $0,p,s}' \
                "$f" > "$tmp"
            first=0
        else
            awk -v p="$f" -v s="$sample" \
                'BEGIN{FS="\t"; OFS="\t"}
                 NR>1 {print $0,p,s}' \
                "$f" >> "$tmp"
        fi
    done < <(find "$study" -type f -name "$pattern" -print0 | sort -z)
done

if [[ "$count" -eq 0 ]]; then
    echo "ERROR: No input files found; existing output was not changed." >&2
    exit 1
fi

mv -- "$tmp" "$out"
trap - EXIT

echo "Files combined: $count"
echo "Output: $(realpath -- "$out")"
echo "Rows: $(($(wc -l < "$out") - 1))"
