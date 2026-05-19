#!/usr/bin/env python3


# To run: python /storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/TWAS/collect_TWAS_hits.py

import csv
from pathlib import Path

INPUT_DIR = "/scratch/cfc5873/mycobiome_TWAS_outputs/all_genes/maize"          # change if needed
OUTPUT_FILE = "/storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/TWAS/maize_688_TWAS_all_sig_hits.csv"
FDR_THRESHOLD = 0.05


all_hits = []
fieldnames = None

for csv_file in Path(INPUT_DIR).glob("TWAS.CMLM*.csv"):
    phenotype = csv_file.stem.replace("TWAS.CMLM_", "")

    with open(csv_file, newline="") as fin:
        reader = csv.DictReader(fin)

        if fieldnames is None:
            fieldnames = reader.fieldnames + ["Phenotype"]

        for row in reader:
            try:
                if float(row["FDR"]) < FDR_THRESHOLD:
                    row["Phenotype"] = phenotype
                    all_hits.append(row)
            except (KeyError, ValueError):
                continue

# Sort alphabetically by Phenotype
all_hits.sort(key=lambda x: x["Phenotype"])

# Write output
with open(OUTPUT_FILE, "w", newline="") as fout:
    writer = csv.DictWriter(fout, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(all_hits)

print(f"Finished writing {len(all_hits)} significant TWAS hits to {OUTPUT_FILE}")
