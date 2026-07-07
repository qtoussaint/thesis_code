#!/usr/bin/env python3
"""Extract CDS protein translations from a GenBank flat file, keyed by /locus_tag.

Usage: gbff_to_protein_fasta.py <in.gbff> <out.faa>
Writes one FASTA record per CDS that has both /locus_tag and /translation.
"""
import sys, re

def main(gbff, out):
    locus = None
    trans_lines = None        # list while inside a /translation, else None
    in_cds = False
    records = {}              # locus_tag -> sequence (first CDS wins)
    with open(gbff) as fh:
        for line in fh:
            # Feature key starts in column 6 (after 5 spaces), e.g. "     CDS  ..."
            m = re.match(r"^ {5}(\S+)\s", line)
            if m:
                # New feature: flush any open translation
                if trans_lines is not None and locus and locus not in records:
                    records[locus] = "".join(trans_lines)
                trans_lines = None
                in_cds = (m.group(1) == "CDS")
                if in_cds:
                    locus = None
                continue
            if not in_cds:
                continue
            lt = re.search(r'/locus_tag="([^"]+)"', line)
            if lt:
                locus = lt.group(1)
                continue
            if "/translation=" in line:
                seq = line.split("/translation=", 1)[1].strip().lstrip('"')
                trans_lines = []
                if seq.endswith('"'):
                    records.setdefault(locus, seq[:-1]) if locus else None
                    if locus and locus not in records:
                        records[locus] = seq[:-1]
                    elif locus:
                        pass
                    trans_lines = None
                else:
                    trans_lines = [seq]
                continue
            if trans_lines is not None:
                seg = line.strip()
                if seg.endswith('"'):
                    trans_lines.append(seg[:-1])
                    if locus and locus not in records:
                        records[locus] = "".join(trans_lines)
                    trans_lines = None
                else:
                    trans_lines.append(seg)
        if trans_lines is not None and locus and locus not in records:
            records[locus] = "".join(trans_lines)

    with open(out, "w") as o:
        for lt, seq in records.items():
            o.write(f">{lt}\n")
            for i in range(0, len(seq), 60):
                o.write(seq[i:i+60] + "\n")
    sys.stderr.write(f"wrote {len(records)} protein records to {out}\n")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
