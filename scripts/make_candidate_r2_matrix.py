#!/usr/bin/env python3
import argparse, os, shutil, subprocess, sys, tempfile
from collections import defaultdict

def read_positions(path):
    positions, seen = [], set()
    with open(path, encoding="utf-8-sig") as f:
        for lineno, line in enumerate(f, 1):
            line = line.strip().replace("\r", "")
            if not line:
                continue
            token = line.split()[0]
            if lineno == 1 and not token.isdigit():
                continue
            if not token.isdigit():
                raise ValueError(f"{path}:{lineno}: not an integer: {token!r}")
            pos = int(token)
            if pos in seen:
                raise ValueError(f"Duplicate position: {pos}")
            seen.add(pos)
            positions.append(pos)
    if not positions:
        raise ValueError("No candidate positions found.")
    return positions

def detect_chrom(cooc):
    r = subprocess.run(["tabix","-l",cooc], check=True, text=True, capture_output=True)
    contigs = [x.strip() for x in r.stdout.splitlines() if x.strip()]
    if len(contigs) == 1:
        return contigs[0]
    raise RuntimeError("Specify --chrom; index contains: " + ",".join(contigs))

def fmt(x):
    return "NA" if x is None else f"{x:.12g}"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cooccurrence")
    ap.add_argument("positions")
    ap.add_argument("-o","--output-prefix", default="candidate_r2")
    ap.add_argument("--chrom")
    args = ap.parse_args()

    if shutil.which("tabix") is None:
        sys.exit("ERROR: tabix not found in PATH")
    if not os.path.exists(args.cooccurrence):
        sys.exit(f"ERROR: file not found: {args.cooccurrence}")
    if not (os.path.exists(args.cooccurrence+".tbi") or os.path.exists(args.cooccurrence+".csi")):
        sys.exit("ERROR: .tbi/.csi index not found")

    try:
        positions = read_positions(args.positions)
        chrom = args.chrom or detect_chrom(args.cooccurrence)
    except Exception as e:
        sys.exit(f"ERROR: {e}")

    pos_set = set(positions)
    pair_rows = defaultdict(list)
    variant_repr = defaultdict(set)

    with tempfile.NamedTemporaryFile("w", delete=False, suffix=".regions.tsv") as tmp:
        region_file = tmp.name
        for pos in positions:
            tmp.write(f"{chrom}\t{pos}\n")

    try:
        p = subprocess.Popen(["tabix","-R",region_file,args.cooccurrence],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        n_queried = n_pair_rows = 0
        for line in p.stdout:
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            c = line.split("\t")
            if len(c) < 13:
                continue
            try:
                p1, p2, r2 = int(c[2]), int(c[8]), float(c[12])
            except ValueError:
                continue
            n_queried += 1
            if p1 in pos_set:
                variant_repr[p1].add((c[0],c[3],c[4]))
            if p2 in pos_set:
                variant_repr[p2].add((c[6],c[9],c[10]))
            if p1 not in pos_set or p2 not in pos_set or p1 == p2:
                continue
            n_pair_rows += 1
            key = tuple(sorted((p1,p2)))
            pair_rows[key].append((r2,c))
        stderr = p.stderr.read()
        rc = p.wait()
        if rc != 0:
            sys.exit("ERROR: tabix failed:\n"+stderr)
    finally:
        try: os.unlink(region_file)
        except OSError: pass

    resolved, ambiguous = {}, {}
    for key, rows in pair_rows.items():
        vals = sorted({x[0] for x in rows})
        if len(vals) == 1:
            resolved[key] = vals[0]
        else:
            resolved[key] = None
            ambiguous[key] = rows

    matrix_path = args.output_prefix + ".r2_matrix.tsv"
    pairs_path = args.output_prefix + ".r2_pairs.tsv"
    summary_path = args.output_prefix + ".r2_summary.tsv"
    qc_path = args.output_prefix + ".variant_qc.tsv"
    amb_path = args.output_prefix + ".ambiguous_pairs.tsv"

    with open(matrix_path,"w",newline="\n") as out:
        out.write("Position\t" + "\t".join(map(str,positions)) + "\n")
        for a in positions:
            cells = []
            for b in positions:
                if a == b:
                    cells.append("1")
                else:
                    cells.append(fmt(resolved.get(tuple(sorted((a,b))))))
            out.write(str(a) + "\t" + "\t".join(cells) + "\n")

    with open(pairs_path,"w",newline="\n") as out:
        out.write("position1\tposition2\tr2\tn_source_rows\tstatus\n")
        for key in sorted(pair_rows):
            val = resolved[key]
            out.write(f"{key[0]}\t{key[1]}\t{fmt(val)}\t{len(pair_rows[key])}\t"
                      f"{'OK' if val is not None else 'AMBIGUOUS_POSITION_PAIR'}\n")

    partners = defaultdict(list)
    for (a,b),v in resolved.items():
        if v is None: continue
        partners[a].append((b,v))
        partners[b].append((a,v))

    with open(summary_path,"w",newline="\n") as out:
        out.write("position\tn_observed_candidate_partners\tn_missing_candidate_pairs\tmax_r2\tmax_r2_partner\tmean_observed_r2\n")
        n = len(positions)
        for pos in positions:
            ps = partners.get(pos,[])
            if ps:
                mp,mr = max(ps,key=lambda x:x[1])
                mean = sum(v for _,v in ps)/len(ps)
            else:
                mp=mr=mean=None
            out.write(f"{pos}\t{len(ps)}\t{n-1-len(ps)}\t{fmt(mr)}\t{mp if mp is not None else 'NA'}\t{fmt(mean)}\n")

    with open(qc_path,"w",newline="\n") as out:
        out.write("position\tn_variant_representations\tstatus\trepresentations\n")
        for pos in positions:
            reps = sorted(variant_repr.get(pos,set()))
            status = "NOT_OBSERVED_IN_QUERY" if len(reps)==0 else ("OK" if len(reps)==1 else "MULTIPLE_ID_REF_ALT")
            txt = ";".join(f"{vid}:{ref}>{alt}" for vid,ref,alt in reps)
            out.write(f"{pos}\t{len(reps)}\t{status}\t{txt}\n")

    if ambiguous:
        with open(amb_path,"w",newline="\n") as out:
            out.write("matrix_position1\tmatrix_position2\tr2\tsource_row\n")
            for key in sorted(ambiguous):
                for r2,c in ambiguous[key]:
                    out.write(f"{key[0]}\t{key[1]}\t{r2:.12g}\t" + "\t".join(c) + "\n")

    print(f"Chromosome label: {chrom}")
    print(f"Candidate positions: {len(positions)}")
    print(f"Rows returned by tabix: {n_queried}")
    print(f"Candidate-candidate source rows: {n_pair_rows}")
    print(f"Resolved candidate pairs: {sum(v is not None for v in resolved.values())}")
    print(f"Ambiguous position pairs: {len(ambiguous)}")
    print(f"Wrote: {matrix_path}")
    print(f"Wrote: {pairs_path}")
    print(f"Wrote: {summary_path}")
    print(f"Wrote: {qc_path}")
    if ambiguous:
        print(f"Wrote: {amb_path}")

if __name__ == "__main__":
    main()
