"""
Pull Allen Mouse Brain Atlas ISH expression for genes underlying the
muscarinic M1 receptor pathway in OPC differentiation and remyelination —
the mechanism identified by Poon et al. (PNAS 2024) for the Green lab's
remyelinating therapy clemastine fumarate.

Gene groups:
  1. Muscarinic ACh receptors (M1–M4):
       Chrm1 (M1 — key therapeutic target, antagonism promotes OPC differentiation)
       Chrm2 (M2), Chrm3 (M3), Chrm4 (M4) — subtype context

  2. OPC / OL lineage markers:
       Olig2  (pan-OL lineage transcription factor)
       Pdgfra (OPC progenitor marker)
       Cspg4  (NG2 — OPC surface marker)
       Myrf   (master myelination TF — marks differentiation)
       Sox10  (OL lineage TF)

  3. Mature myelin:
       Mbp, Plp1, Cnp (mature OL / myelin)

  4. Neuroaxonal health:
       Nefl (neurofilament light — blood biomarker the lab uses to track neurodegeneration)

Structures: 8 brain regions with reliable ISH coverage, spanning
gray/white matter gradients and MS lesion-prone areas.
  Isocortex (315), Hippocampus (1089), Cerebellum (512), Striatum (672),
  Hypothalamus (1097), Midbrain (313), Thalamus (549), Olfactory bulb (507)
Note: Corpus callosum (776) excluded — ISH returns 0 for fiber tracts (no cell bodies).
"""
import requests, csv, time
from pathlib import Path

OUT  = Path(__file__).parent
BASE = "http://api.brain-map.org/api/v2/data/query.json"

GENES = {
    "muscarinic": ["Chrm1", "Chrm2", "Chrm3", "Chrm4"],
    "opc_lineage": ["Olig2", "Pdgfra", "Cspg4", "Myrf", "Sox10"],
    "myelin": ["Mbp", "Plp1", "Cnp"],
    "neuroaxonal": ["Nefl"],
}
ALL_GENES = [g for grp in GENES.values() for g in grp]

STRUCTURES = {
    "Isocortex":   315,
    "Hippocampus": 1089,
    "Cerebellum":  512,
    "Striatum":    672,
    "Hypothalamus": 1097,
    "Midbrain":    313,
    "Thalamus":    549,
    "Olfact_bulb": 507,
}

def get_dataset(gene):
    url = (BASE + "?criteria=model::SectionDataSet,"
           "rma::criteria,genes[acronym$eq" + gene + "],"
           "[failed$eqfalse],products[abbreviation$eqMouse]&num_rows=5")
    r = requests.get(url, timeout=20).json()
    msgs = r.get("msg", [])
    if not isinstance(msgs, list) or not msgs:
        return None
    for row in msgs:
        if row.get("plane_of_section_id") == 2:
            return row["id"]
    return msgs[0]["id"]

def get_expression(dataset_id, struct_ids):
    ids = ",".join(str(i) for i in struct_ids)
    url = (BASE + "?criteria=model::StructureUnionize,"
           "rma::criteria,[section_data_set_id$eq" + str(dataset_id) + "],"
           "[structure_id$in" + ids + "]&num_rows=50")
    r = requests.get(url, timeout=20).json()
    return {row["structure_id"]: row.get("expression_energy", 0.0) or 0.0
            for row in r.get("msg", [])}

struct_ids  = list(STRUCTURES.values())
id_to_name  = {v: k for k, v in STRUCTURES.items()}
gene_to_grp = {g: grp for grp, gs in GENES.items() for g in gs}

rows = []
for gene in ALL_GENES:
    print(f"  {gene}...", end=" ", flush=True)
    ds = get_dataset(gene)
    if ds is None:
        print("no dataset"); continue
    expr = get_expression(ds, struct_ids)
    row = {"gene": gene, "group": gene_to_grp[gene], "dataset_id": ds}
    for sid, name in id_to_name.items():
        row[name] = round(expr.get(sid, 0.0), 4)
    rows.append(row)
    print(f"ds={ds} ✓")
    time.sleep(0.3)

fieldnames = ["gene", "group", "dataset_id"] + list(STRUCTURES.keys())
with open(OUT / "expression.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    w.writerows(rows)

print(f"\nSaved expression.tsv ({len(rows)} genes)")
for r in rows:
    vals = {k: r[k] for k in STRUCTURES}
    top = max(vals, key=vals.get)
    print(f"  {r['gene']:10s}  peak={top} ({vals[top]:.3f})")
