#!/usr/bin/env python3
"""
Build a snpEff-style annotations table for unitig coordinates by interval-overlap
against the ATCC_700669 (FM211187) gene model — the SAME GenBank snpEff was built
from — so unitig gene labels use snpEff's gene naming (/gene, else /locus_tag), matching
the SNP analysis and the genes-of-interest list.

The gwas_workflow pipeline joins genes to plot points by EXACT position
(`match(variant_positions, ann$POS)` in pipeline.R) and reads the column `ANN[*].GENE`.
So we emit one row per distinct unitig coordinate with POS = that coordinate and
ANN[*].GENE = the gene whose interval contains it (intergenic / unmapped -> MODIFIER).
Header matches genotype/fields_filtered_maf05_multiallelic.txt so read.delim mangles it
to the same column names the pipeline expects.
"""
import argparse
import bisect
import re
import sys

HEADER = "POS\tREF\tALT\tANN[*].EFFECT\tANN[*].IMPACT\tANN[*].GENE\tANN[*].HGVS_C\tANN[*].HGVS_P\tLOF[*].GENE"


def parse_gbff_genes(path, feature="gene"):
    """Return sorted list of (start, end, gene) from GenBank `feature` records.
    gene name = /gene if present else /locus_tag."""
    intervals = []
    in_features = False
    cur_loc = cur_gene = cur_locus = None
    collecting_loc = False

    def flush():
        nonlocal cur_loc, cur_gene, cur_locus, collecting_loc
        if cur_loc is not None:
            nums = [int(x) for x in re.findall(r"\d+", cur_loc)]
            g = cur_gene or cur_locus
            if nums and g:
                intervals.append((min(nums), max(nums), g))
        cur_loc = cur_gene = cur_locus = None
        collecting_loc = False

    with open(path) as fh:
        for line in fh:
            if line.startswith("FEATURES"):
                in_features = True
                continue
            if not in_features:
                continue
            if line.startswith("ORIGIN") or (line and line[0] not in " \t"):
                flush()
                if line.startswith("ORIGIN"):
                    break
                continue
            is_feature_line = line[:5] == "     " and len(line) > 5 and line[5] != " "
            if is_feature_line:
                flush()
                key = line[5:21].strip()
                if key == feature:
                    cur_loc = line[21:].strip()
                    collecting_loc = True
            else:
                stripped = line.strip()
                if cur_loc is not None and collecting_loc and not stripped.startswith("/"):
                    cur_loc += stripped
                if stripped.startswith("/"):
                    collecting_loc = False
                    m = re.match(r'/gene="([^"]+)"', stripped)
                    if m and cur_gene is None:
                        cur_gene = m.group(1)
                    m = re.match(r'/locus_tag="([^"]+)"', stripped)
                    if m and cur_locus is None:
                        cur_locus = m.group(1)
    flush()
    intervals.sort()
    return intervals


def make_lookup(intervals):
    starts = [iv[0] for iv in intervals]

    def lookup(pos):
        i = bisect.bisect_right(starts, pos) - 1
        # genes are sorted by start and mostly non-overlapping; check a few back
        # to handle minor overlaps / nested features.
        j = i
        while j >= 0 and j > i - 6:
            s, e, g = intervals[j]
            if s <= pos <= e:
                return g
            j -= 1
        return None

    return lookup


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gbff", required=True)
    ap.add_argument("--varindex", required=True, help="spn_unitigs_mapped_variant_index.csv (variant_name,position)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--feature", default="gene")
    ap.add_argument("--validate", help="snpEff fields_*.txt (POS<tab>...GENE col6) to check naming against", default=None)
    args = ap.parse_args()

    intervals = parse_gbff_genes(args.gbff, args.feature)
    sys.stderr.write("Parsed %d '%s' features from gbff\n" % (len(intervals), args.feature))
    lookup = make_lookup(intervals)

    if args.validate:
        n = match = mism = nogene = 0
        examples = []
        with open(args.validate) as fh:
            header = fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 6 or not f[5]:
                    continue
                try:
                    pos = int(f[0])
                except ValueError:
                    continue
                truth = f[5]
                got = lookup(pos)
                n += 1
                if got is None:
                    nogene += 1
                    if len(examples) < 8:
                        examples.append((pos, truth, "None"))
                elif got == truth:
                    match += 1
                else:
                    mism += 1
                    if len(examples) < 8:
                        examples.append((pos, truth, got))
        sys.stderr.write("VALIDATION vs %s: n=%d match=%d (%.2f%%) mismatch=%d no-gene=%d\n"
                         % (args.validate, n, match, 100.0 * match / max(n, 1), mism, nogene))
        for pos, truth, got in examples:
            sys.stderr.write("  POS %d: snpEff=%s  lookup=%s\n" % (pos, truth, got))

    # distinct positions from the variant index
    positions = set()
    with open(args.varindex) as fh:
        header = fh.readline()
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            p = line.rsplit(",", 1)[-1]
            try:
                positions.add(int(p))
            except ValueError:
                pass

    n_gene = n_mod = 0
    with open(args.out, "w") as out:
        out.write(HEADER + "\n")
        for pos in sorted(positions):
            gene = None if pos < 0 else lookup(pos)
            if gene is None:
                gene = "MODIFIER"
                n_mod += 1
            else:
                n_gene += 1
            out.write("%d\t\t\t\t\t%s\t\t\t\n" % (pos, gene))
    sys.stderr.write("Wrote %d annotation rows (%d in-gene, %d MODIFIER/unmapped) to %s\n"
                     % (n_gene + n_mod, n_gene, n_mod, args.out))


if __name__ == "__main__":
    main()
