from collections import Counter
import sys

vcf_path = sys.argv[1]
output_path = sys.argv[2]
maf_threshold = float(sys.argv[3])
het_threshold = float(sys.argv[4])

DEBUG_LIMIT = 20
debug_count = 0

def convert_gt(g, ref, alt):

    numeric = {"0/0","0/1","1/0","1/1",".|.","./.","0|0","0|1","1|0","1|1"}

    if g in numeric:
        return g

    if g == "./." or g == ".|.":
        return "./."

    # Determine separator
    if "/" in g:
        a1, a2 = g.split("/")
        sep = "/"
    elif "|" in g:
        a1, a2 = g.split("|")
        sep = "|"
    else:
        return "./."

    # Convert allele-coded → numeric
    if a1 == ref and a2 == ref:
        return f"0{sep}0"
    if a1 == alt and a2 == alt:
        return f"1{sep}1"
    if {a1, a2} == {ref, alt}:
        return f"0{sep}1"

    # Anything unexpected
    return "./."


with open(vcf_path) as vcf, open(output_path, "w") as out:

    for line in vcf:
        if line.startswith("#"):
            out.write(line)
            continue

        fields = line.rstrip().split("\t")
        chrom, pos, vid, ref, alt = fields[:5]
        raw_gts = fields[9:]

        # Convert genotypes
        gts = [convert_gt(g, ref, alt) for g in raw_gts]

        # Count categories
        gc = Counter(gts)
        hom_ref = gc["0/0"] + gc["0|0"]
        hom_alt = gc["1/1"] + gc["1|1"]
        het = gc["0/1"] + gc["1/0"] + gc["0|1"] + gc["1|0"]

        total = len(gts)
        hom_maf = min(hom_ref, hom_alt) / total if total else 0
        het_freq = het / total if total else 0

        # ---------------------
        # DEBUGGING OUTPUT
        # ---------------------
        if debug_count < DEBUG_LIMIT:
            print(f"\n--- DEBUG variant {debug_count+1} ---")
            print(f"ID: {vid}   REF={ref}   ALT={alt}")
            print(f"Raw GTs: {raw_gts[:10]} ...")
            print(f"Converted GTs: {gts[:10]} ...")
            print(f"Counts: hom_ref={hom_ref}, hom_alt={hom_alt}, het={het}, total={total}")
            print(f"hom_maf = {hom_maf:.4f}, het_freq = {het_freq:.4f}")
            print(f"PASS? {hom_maf > maf_threshold and het_freq < het_threshold}")
            print("--------------------------------------")
            debug_count += 1

        # Filtering condition
        if hom_maf > maf_threshold and het_freq < het_threshold:
            out.write("\t".join(fields[:9] + gts) + "\n")
