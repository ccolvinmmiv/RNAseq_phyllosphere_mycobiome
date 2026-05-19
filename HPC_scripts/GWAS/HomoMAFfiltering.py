from collections import Counter
import sys

# Inputs
vcf_path = sys.argv[1]
output_path = sys.argv[2]
maf_threshold = float(sys.argv[3])
het_threshold = float(sys.argv[4])


def convert_gt(g, ref, alt):
    """
    Convert genotypes to numeric format (0/0, 0/1, 1/1).
    If already numeric, return unchanged.
    """

    # Numeric GT formats → return unchanged
    if g in {"0/0","0/1","1/0","1/1","./.","0|0","0|1","1|0","1|1",".|."}:
        return g

    # Missing or bad data
    if g == "./." or g == ".|.":
        return "./."

    # Split genotype
    if "/" in g:
        a1, a2 = g.split("/")
        sep = "/"
    elif "|" in g:
        a1, a2 = g.split("|")
        sep = "|"
    else:
        return "./."  # Fallback if format unexpected

    # Handle non-numeric → convert
    if a1 == ref and a2 == ref:
        return f"0{sep}0"
    if a1 == alt and a2 == alt:
        return f"1{sep}1"
    if {a1, a2} == {ref, alt}:
        return f"0{sep}1"

    # Anything else → missing
    return "./."


# Process VCF
with open(vcf_path) as vcf, open(output_path, 'w') as out:
    for line in vcf:
        line = line.rstrip()

        # Write headers directly
        if line.startswith("#"):
            out.write(line + "\n")
            continue

        fields = line.split("\t")
        genotypes = fields[9:]

        # Count genotype types
        g_counts = Counter(genotypes)
        hom_ref = g_counts.get("0/0", 0) + g_counts.get("0|0", 0)
        hom_alt = g_counts.get("1/1", 0) + g_counts.get("1|1", 0)
        het = (
            g_counts.get("0/1", 0) + g_counts.get("0|1", 0) +
            g_counts.get("1/0", 0) + g_counts.get("1|0", 0)
        )

        total = len(genotypes)
        hom_maf = min(hom_ref, hom_alt) / total if total > 0 else 0
        het_freq = het / total if total > 0 else 0

        # Apply filters
        if hom_maf > maf_threshold and het_freq < het_threshold:
            out.write(line + "\n")