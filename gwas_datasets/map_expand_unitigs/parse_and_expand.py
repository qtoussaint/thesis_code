#!/usr/bin/env python3
"""
Parse bwa-mem SAM of unitigs aligned to FM211187 and build a coordinate-expanded
presence/absence matrix.

For every unitig we collect ALL reference coordinates it maps to (primary + secondary
alignments, near full-length, de-duplicated). Then we build a presence/absence matrix in
which a unitig mapping to k coordinates becomes k identical columns named
``<unitig_id>_1 .. <unitig_id>_k``, each tagged with one coordinate. Unitigs that do not
map anywhere are KEPT as a single column ``<unitig_id>_1`` with the sentinel coordinate
-10000 (so they stay in the GWAS and cluster at the far left of the Manhattan).

The expanded matrix is built from the MAF-filtered unitig set (the set that feeds the
GWAS); the per-unitig coordinate table is written for every unitig regardless of MAF.

Outputs:
  - coords.tsv          : unitig_id <tab> n_hits <tab> comma-separated sorted coordinates
  - mapped.rtab         : expanded presence/absence (expanded_name <tab> sample 0/1 ...)
  - variant_index.csv   : expanded_name,position   (position = coordinate, -10000 if unmapped)
"""
import argparse
import sys


def aligned_query_len(cigar):
    """Length of the query consumed by alignment ops M/I/=/X (excludes clips S/H/D/N/P)."""
    if cigar == "*":
        return 0
    n = 0
    num = 0
    for ch in cigar:
        if ch.isdigit():
            num = num * 10 + (ord(ch) - 48)
        else:
            if ch in "MI=X":
                n += num
            num = 0
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sam", required=True)
    ap.add_argument("--unitig_ids", required=True, help="spn_unitigs.unitig_ids.tsv (id<tab>seq)")
    ap.add_argument("--rtab", required=True, help="spn_unitigs.rtab (Unitig_sequence + samples)")
    ap.add_argument("--out_coords", required=True)
    ap.add_argument("--out_rtab", required=True)
    ap.add_argument("--out_varindex", required=True)
    ap.add_argument("--min_af", type=float, default=0.05)
    ap.add_argument("--max_af", type=float, default=0.95)
    ap.add_argument("--min_cov", type=float, default=0.8,
                    help="min fraction of the unitig that must align to count a hit")
    ap.add_argument("--unmapped_pos", type=int, default=-10000)
    args = ap.parse_args()

    # 1. unitig_id <-> sequence and unitig length
    id2seq = {}
    id2len = {}
    seq2id = {}
    with open(args.unitig_ids) as fh:
        for line in fh:
            uid, seq = line.rstrip("\n").split("\t")
            id2seq[uid] = seq
            id2len[uid] = len(seq)
            seq2id[seq] = uid
    sys.stderr.write("Loaded %d unitig ids\n" % len(id2seq))

    # 2. Parse SAM -> per-unitig set of leftmost reference positions
    id2pos = {}
    n_aln = 0
    with open(args.sam) as fh:
        for line in fh:
            if line[0] == "@":
                continue
            f = line.split("\t", 6)  # qname,flag,rname,pos,mapq,cigar,rest
            if len(f) < 6:
                continue
            flag = int(f[1])
            if flag & 0x4:  # unmapped
                continue
            qname = f[0]
            pos = int(f[3])
            cigar = f[5]
            ulen = id2len.get(qname)
            if ulen is None:
                continue
            if aligned_query_len(cigar) / ulen < args.min_cov:
                continue
            id2pos.setdefault(qname, set()).add(pos)
            n_aln += 1
    sys.stderr.write("Kept %d alignments; %d unitigs have >=1 hit\n" % (n_aln, len(id2pos)))

    # 3. coords.tsv for EVERY unitig (n_hits = 0 if unmapped)
    with open(args.out_coords, "w") as out:
        out.write("unitig_id\tn_hits\tcoordinates\n")
        for uid in id2seq:
            coords = sorted(id2pos.get(uid, ()))
            out.write("%s\t%d\t%s\n" % (uid, len(coords), ",".join(map(str, coords))))

    # 4. Stream rtab: MAF-filter, expand each kept unitig into one column per coordinate
    n_in = n_kept = n_cols = n_unmapped_kept = 0
    with open(args.rtab) as fh, \
         open(args.out_rtab, "w") as out_mat, \
         open(args.out_varindex, "w") as out_vi:
        header = fh.readline().rstrip("\n")
        samples = header.split("\t")[1:]
        nsamp = len(samples)
        out_mat.write("unitig\t" + "\t".join(samples) + "\n")
        out_vi.write("variant_name,position\n")
        for line in fh:
            n_in += 1
            tab = line.find("\t")
            seq = line[:tab]
            rest = line[tab + 1:].rstrip("\n")  # genotype fields, tab-separated 0/1
            ac = rest.count("1")
            af = ac / nsamp
            if af < args.min_af or af > args.max_af:
                continue
            n_kept += 1
            uid = seq2id.get(seq)
            if uid is None:
                # sequence not in id map (shouldn't happen); skip
                continue
            coords = sorted(id2pos.get(uid, ()))
            if not coords:
                coords = [args.unmapped_pos]
                n_unmapped_kept += 1
            for k, coord in enumerate(coords, 1):
                name = "%s_%d" % (uid, k)
                out_mat.write(name + "\t" + rest + "\n")
                out_vi.write("%s,%d\n" % (name, coord))
                n_cols += 1

    sys.stderr.write(
        "rtab rows read: %d; MAF-passing unitigs: %d; expanded columns written: %d; "
        "MAF-passing unmapped (pos=%d): %d\n"
        % (n_in, n_kept, n_cols, args.unmapped_pos, n_unmapped_kept))


if __name__ == "__main__":
    main()
