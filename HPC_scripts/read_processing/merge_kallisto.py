#!/usr/bin/env python3
import os
import sys

#First run:
#   chmod +x /storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/read_processing/merge_kallisto.py
#Then:
#  python /storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/read_processing/merge_kallisto.py /scratch/cfc5873/kallisto_out/ne_2020_maize

def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: python merge_kallisto.py /path/to/kallisto_out")

    mydir = sys.argv[1]
    if not os.path.isdir(mydir):
        sys.exit(f"ERROR: {mydir} is not a valid directory.")

    gene_exp = {}     # transcript -> {sample: tpm}
    samples = []      # maintain original order

    # Scan each subdirectory
    for entry in os.scandir(mydir):
        if not entry.is_dir():
            continue

        sample = entry.name
        abundance_path = os.path.join(entry.path, "abundance.tsv")

        if not os.path.exists(abundance_path):
            continue

        samples.append(sample)

        with open(abundance_path) as fh:
            next(fh)  # skip header
            for line in fh:
                fields = line.rstrip().split("\t")
                transcript = fields[0]
                tpm = float(fields[-1])

                if transcript not in gene_exp:
                    gene_exp[transcript] = {}

                gene_exp[transcript][sample] = tpm

    # Write output
    output_file = "/storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/read_processing/ne_2020_maize_merged_gene_tpms.csv"
    samples_sorted = sorted(samples)

    with open(output_file, "w") as out:
        out.write("TranscriptID," + ",".join(samples_sorted) + "\n")

        for transcript in sorted(gene_exp):
            values = [str(gene_exp[transcript].get(s, 0.0)) for s in samples_sorted]
            out.write(transcript + "," + ",".join(values) + "\n")

    print(f"✔ Merged file written to: {output_file}")

if __name__ == "__main__":
    main()


#Now getting raw counts instead of TPM

def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: python merge_kallisto.py /path/to/kallisto_out")

    mydir = sys.argv[1]
    if not os.path.isdir(mydir):
        sys.exit(f"ERROR: {mydir} is not a valid directory.")

    gene_exp = {}     # transcript -> {sample: tpm}
    samples = []      # maintain original order

    # Scan each subdirectory
    for entry in os.scandir(mydir):
        if not entry.is_dir():
            continue

        sample = entry.name
        abundance_path = os.path.join(entry.path, "abundance.tsv")

        if not os.path.exists(abundance_path):
            continue

        samples.append(sample)

        with open(abundance_path) as fh:
            next(fh)  # skip header
            for line in fh:
                fields = line.rstrip().split("\t")
                transcript = fields[0]
                tpm = float(fields[-2])

                if transcript not in gene_exp:
                    gene_exp[transcript] = {}

                gene_exp[transcript][sample] = tpm

    # Write output
    output_file = "/storage/work/cfc5873/Mycobiome/scripts/ne_2020_maize/read_processing/ne_2020_maize_merged_gene_raw_counts.csv"
    samples_sorted = sorted(samples)

    with open(output_file, "w") as out:
        out.write("TranscriptID," + ",".join(samples_sorted) + "\n")

        for transcript in sorted(gene_exp):
            values = [str(gene_exp[transcript].get(s, 0.0)) for s in samples_sorted]
            out.write(transcript + "," + ",".join(values) + "\n")

    print(f"✔ Merged file written to: {output_file}")

if __name__ == "__main__":
    main()