#!/usr/bin/env python3
"""Order a FOF by taxonomic lineage so related genomes are adjacent.

usage: sort_fof_by_lineage.py LINEAGE_TSV IN.fof OUT.fof

LINEAGE_TSV is a tab-separated table whose first column is the sample id and
whose columns 5-10 are phylum, class, order, family, genus and species (the
661k-bacteria File2_taxid_lineage layout; a header line is skipped). The FOF
has one "name<TAB>path" line per dataset. Samples without a lineage go last.

The streaming planner merges each dataset into the best existing group, so an
input order that keeps similar genomes together gives smaller groups, fewer
candidate evaluations and a smaller index.
"""
import sys


def main():
    lineage_tsv, fof_in, fof_out = sys.argv[1:4]
    lineage = {}
    with open(lineage_tsv) as f:
        next(f)
        for line in f:
            cols = line.rstrip("\n").split("\t")
            lineage[cols[0]] = tuple(cols[4:10])
    rows = [line.rstrip("\n") for line in open(fof_in) if line.strip()]

    def key(row):
        name = row.split("\t")[0]
        return (name not in lineage, lineage.get(name, ()), name)

    rows.sort(key=key)
    with open(fof_out, "w") as f:
        f.write("\n".join(rows) + "\n")
    missing = sum(1 for r in rows if r.split("\t")[0] not in lineage)
    print(f"{fof_out}: {len(rows)} datasets, {missing} without lineage")


if __name__ == "__main__":
    main()
