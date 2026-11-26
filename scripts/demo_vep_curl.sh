#!/usr/bin/env bash
# Only the public BRCA2 c.9217G>T educational allele. No accounts or input files.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "$0")" && pwd)"
result_dir="$(mktemp -d "${TMPDIR:-/tmp}/brca-vep-demo.XXXXXX")"
cd -- "$result_dir"
printf '公开 BRCA2 示例；结果目录：%s\n' "$result_dir"
endpoint='https://rest.ensembl.org/vep/human'
options='hgvs=1&merged=1&transcript_version=1&numbers=1&mane=1&variant_class=1'
fetch_json() {
  local file="$1" status
  shift
  if ! status="$(curl -sS --max-time 30 -D "$file.headers" -o "$file.json" -w '%{http_code}' -H 'Accept: application/json' "$@")"; then
    printf '%s：网络未完成；请查看结果目录。\n' "$file" >&2
    return 1
  fi
  printf '%s：HTTP %s\n' "$file" "$status"
  if [[ "$status" != 200 ]]; then
    printf '已保留响应头和正文；HTTP %s 不能当作未检出。\n' "$status" >&2
    return 1
  fi
}
failed=0
fetch_json vep-vcf -H 'Content-Type: application/json' --data-binary '{"variants":["13 32380106 . G T . . ."]}' "$endpoint/region?$options" || failed=1
fetch_json vep-c "$endpoint/hgvs/NM_000059.4%3Ac.9217G%3ET?$options" || failed=1
if [[ "$failed" != 0 ]]; then exit 1; fi
python3 "$script_dir/check_public_vep.py"
