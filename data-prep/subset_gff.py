# Keep the gene hierarchies (gene -> transcript -> exon/CDS/UTR) of genes ending before END.
import gzip, sys
src, dst, end = sys.argv[1], sys.argv[2], int(sys.argv[3])
rows = []
for l in gzip.open(src, "rt"):
    if l.startswith("#") or not l.strip():
        continue
    r = l.rstrip("\n").split("\t")
    attrs = dict(kv.split("=", 1) for kv in r[8].split(";") if "=" in kv)
    rows.append((r, attrs))
by_id = {a["ID"]: a for r, a in rows if "ID" in a}
def root(a):
    while "Parent" in a and a["Parent"] in by_id:
        a = by_id[a["Parent"]]
    return a.get("ID")
keep = {a["ID"] for r, a in rows if "ID" in a and "Parent" not in a and r[2] != "chromosome" and int(r[4]) <= end}
with open(dst, "w") as out:
    out.write("##gff-version 3\n")
    for r, a in rows:
        if root(a) in keep:
            out.write("\t".join(r) + "\n")
print("genes kept:", len(keep))
