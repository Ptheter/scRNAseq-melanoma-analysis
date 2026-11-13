#!/usr/bin/env bash
# Reproduce the GSE120575 analysis: download the two GEO files (if missing), then run the numbered scripts in order.
# usage: bash run_all.sh   (runs 01-03 and 08-17; the exploratory scripts 04-07 live in scripts/exploratory/ and are not run)
# Logs go to logs/. Run inside the conda environment from environment.yml.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p data/raw data/processed results docs/figures logs

BASE=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120575/suppl
for f in GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz GSE120575_patient_ID_single_cells.txt.gz; do
  if [ ! -s "data/raw/$f" ]; then
    echo "downloading $f"
    curl -fL -o "data/raw/$f" "$BASE/$f"
  fi
  gzip -t "data/raw/$f"
done

steps="01 02 03 08 09 10 11 12 13 14 15 16 17"

for n in $steps; do
  script=$(ls scripts/${n}_*.py)
  echo "=== $script"
  python "$script" 2>&1 | tee "logs/${n}.log"
done
echo "done. Results in results/, figures in docs/figures/, logs in logs/"
