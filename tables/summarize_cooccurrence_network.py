#!/usr/bin/env python3
import argparse
import re
from collections import defaultdict

MISSING = {"", ".", "NA", "N/A", "na", "NaN", "nan"}

def natural_key(x):
    parts = re.split(r"(\d+)", str(x))
    return tuple(int(p) if p.isdigit() else p for p in parts)

def looks_like_header(first_field):
    x = first_field.strip().lower()
    return "position" in x or "variant" in x or x in {"pos", "position_of_variant"}

def read_input(path):
    rows = []
    with open(path, "r", encoding="utf-8-sig") as f:
        for line_no, line in enumerate(f, 1):
            line = line.rstrip("\r\n")
            if not line.strip():
                continue
            cols = line.split("\t")
            target = cols[0].strip()
            if not rows and looks_like_header(target):
                continue
            if target in MISSING:
                continue
            partners = []
            seen = set()
            for value in cols[1:]:
                value = value.strip()
                if value in MISSING or value == target:
                    continue
                if value not in seen:
                    partners.append(value)
                    seen.add(value)
            rows.append((target, partners, line_no))
    if not rows:
        raise ValueError("No data rows were found.")

    order = []
    merged = defaultdict(list)
    merged_seen = defaultdict(set)
    for target, partners, _ in rows:
        if target not in merged:
            order.append(target)
        for p in partners:
            if p not in merged_seen[target]:
                merged[target].append(p)
                merged_seen[target].add(p)
    return order, dict(merged)

class UnionFind:
    def __init__(self, items):
        self.parent = {x: x for x in items}
        self.rank = {x: 0 for x in items}
    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x
    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1

def main():
    ap = argparse.ArgumentParser(description="Summarize a variable-width co-occurrence table and group target variants.")
    ap.add_argument("input", help="Input variable-width TSV")
    ap.add_argument("-o", "--output-prefix", default="cooccurrence_network", help="Output prefix")
    args = ap.parse_args()

    targets, outgoing_all = read_input(args.input)
    target_set = set(targets)
    input_index = {v: i for i, v in enumerate(targets)}

    outgoing_target = {}
    incoming_target = defaultdict(set)
    for target in targets:
        out = {p for p in outgoing_all.get(target, []) if p in target_set and p != target}
        outgoing_target[target] = out
        for p in out:
            incoming_target[p].add(target)

    direct = {}
    for target in targets:
        direct[target] = (outgoing_target.get(target, set()) | incoming_target.get(target, set())) - {target}

    uf = UnionFind(targets)
    for a in targets:
        for b in direct[a]:
            uf.union(a, b)

    components_by_root = defaultdict(list)
    for v in targets:
        components_by_root[uf.find(v)].append(v)
    component_lists = list(components_by_root.values())
    component_lists.sort(key=lambda members: min(input_index[x] for x in members))

    component_id = {}
    component_members = {}
    for i, members in enumerate(component_lists, 1):
        cid = f"C{i:03d}"
        members = sorted(members, key=lambda x: input_index[x])
        component_members[cid] = members
        for v in members:
            component_id[v] = cid

    edge_rows = []
    seen_edges = set()
    for a in targets:
        for b in direct[a]:
            key = tuple(sorted((a, b), key=natural_key))
            if key in seen_edges:
                continue
            seen_edges.add(key)
            a_to_b = b in outgoing_target.get(a, set())
            b_to_a = a in outgoing_target.get(b, set())
            if a_to_b and b_to_a:
                direction = "BOTH"
            elif a_to_b:
                direction = f"{a}->{b}"
            elif b_to_a:
                direction = f"{b}->{a}"
            else:
                direction = "UNKNOWN"
            edge_rows.append((key[0], key[1], direction))

    edge_count_by_component = defaultdict(int)
    for a, b, _ in edge_rows:
        cid = component_id[a]
        if component_id[b] == cid:
            edge_count_by_component[cid] += 1

    summary_path = args.output_prefix + ".summary.tsv"
    components_path = args.output_prefix + ".components.tsv"
    edges_path = args.output_prefix + ".edges.tsv"

    with open(summary_path, "w", encoding="utf-8", newline="\n") as out:
        out.write("variant\tn_outgoing_all\tn_outgoing_target\tn_incoming_target\tn_direct_target_partners\tdirection_status\tcomponent_id\tcomponent_size\toutgoing_to_targets\tincoming_from_targets\tdirect_target_partners\tall_outgoing_partners\n")
        for v in targets:
            out_t = outgoing_target.get(v, set())
            inc_t = incoming_target.get(v, set())
            partners = direct[v]
            if out_t and inc_t:
                status = "BOTH"
            elif out_t:
                status = "OUTGOING_ONLY"
            elif inc_t:
                status = "INCOMING_ONLY"
            else:
                status = "ISOLATED_WITHIN_TARGETS"
            cid = component_id[v]
            out.write(
                f"{v}\t{len(outgoing_all.get(v, []))}\t{len(out_t)}\t{len(inc_t)}\t{len(partners)}\t{status}\t{cid}\t{len(component_members[cid])}\t"
                f"{','.join(sorted(out_t, key=natural_key))}\t{','.join(sorted(inc_t, key=natural_key))}\t{','.join(sorted(partners, key=natural_key))}\t{','.join(outgoing_all.get(v, []))}\n"
            )

    with open(components_path, "w", encoding="utf-8", newline="\n") as out:
        out.write("component_id\tcomponent_size\tn_edges\tcomponent_status\tmembers\n")
        for cid, members in component_members.items():
            status = "SINGLETON" if len(members) == 1 else "MULTI"
            out.write(f"{cid}\t{len(members)}\t{edge_count_by_component[cid]}\t{status}\t{','.join(members)}\n")

    with open(edges_path, "w", encoding="utf-8", newline="\n") as out:
        out.write("variant1\tvariant2\tdirection\n")
        for a, b, direction in sorted(edge_rows, key=lambda x: (input_index.get(x[0], 10**9), input_index.get(x[1], 10**9))):
            out.write(f"{a}\t{b}\t{direction}\n")

    print(f"Wrote: {summary_path}")
    print(f"Wrote: {components_path}")
    print(f"Wrote: {edges_path}")

if __name__ == "__main__":
    main()
