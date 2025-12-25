import pandas as pd
import json

def union_to_network_json(edge_union,
                          node_attrs=None,
                          default_style=None,
                          default_source_occupancy=1,
                          default_other_occupancy=0):
    """
    edge_union: iterable of (source, target) pairs (tuples or lists)
    node_attrs: optional dict mapping node_id -> dict of extra attributes (overrides defaults)
    default_style: optional list for edge style (default ["solid","triangle"])
    Returns: dict {"nodes": [...], "edges": [...]}
    """
    if default_style is None:
        default_style = ["solid", "triangle"]
    node_attrs = node_attrs or {}

    # normalize edges to tuple form
    edges = [tuple(e) for e in edge_union]

    # collect nodes and which are sources
    nodes_set = set()
    sources = set()
    for s, t in edges:
        nodes_set.add(s)
        nodes_set.add(t)
        sources.add(s)

    # build nodes list (sorted for reproducibility)
    nodes = []
    for node in sorted(nodes_set):
        attrs = dict(node_attrs.get(node, {}))  # copy overrides if provided
        if "sua7Occupancy" not in attrs:
            attrs["sua7Occupancy"] = (default_source_occupancy
                                      if node in sources
                                      else default_other_occupancy)
        node_obj = {"id": node, "label": node}
        node_obj.update(attrs)
        nodes.append(node_obj)

    # build edges list in original union order
    edges_out = []
    for s, t in edges:
        edge_obj = {
            "source": s,
            "target": t,
            "label": "",
            "style": list(default_style),
        }
        edges_out.append(edge_obj)

    return {"nodes": nodes, "edges": edges_out}


infile  = "./data/Supplementary_Data_2.xlsx"
outfile = "./data/extracted_CT_to_IS_rows6-5386.xlsx"
sheet   = "Pol2_Promoters_x_TF_Motif_Cof"

df = pd.read_excel(infile, sheet_name=sheet, usecols="CT:IS", header=None, engine="openpyxl")
df2 = pd.read_excel(infile, sheet_name=sheet, usecols="A:B", header=None, engine="openpyxl")

# Now Excel rows 6..5386 correspond to df rows 5..5385
subset = df.iloc[5:5386]
subset2 = df2.iloc[5:5386]

# check lengths match
if len(subset) != len(subset2):
    raise ValueError("Row counts differ: {} vs {}".format(len(subset2), len(subset)))

# concat side-by-side, preserving row index
combined = pd.concat([subset2, subset], axis=1)

df = combined.reset_index(drop=True).copy()

# Ensure we have at least two columns (A and B)
if df.shape[1] < 2:
    raise ValueError("combined must have at least 2 columns (A and B)")

# Ensure at least two rows exist to move A2->A1 and B2->B1
if len(df) >= 2:
    # Move A2->A1 and B2->B1 (index 1 -> index 0)
    df.iat[0, 0] = df.iat[1, 0]   # column A (first column)
    df.iat[0, 1] = df.iat[1, 1]   # column B (second column)

    # Write "TF" into A2 and B2 (index 1)
    df.iat[1, 0] = "TFs"
    df.iat[1, 1] = "TFs"
else:
    # if only one row exists, just set A1/B1 to "TF" (or raise)
    df.iat[0, 0] = "TFs"
    df.iat[0, 1] = "TFs"

# Remove row3 entirely if it exists (Excel row 3 -> index 2)
if len(df) > 2:
    df = df.drop(index=2).reset_index(drop=True)

# now df is the modified DataFrame; write to Excel
df.to_excel(outfile, index=False, header=False)
print("Wrote", outfile)

FNAME = "./data/extracted_CT_to_IS_rows6-5386.xlsx"
SHEET = 0

# read file with two header rows -> MultiIndex columns
if FNAME.lower().endswith(".csv"):
    df = pd.read_csv(FNAME, header=[0, 1], dtype=object)
else:
    df = pd.read_excel(FNAME, header=[0, 1], sheet_name=SHEET, dtype=object)

# normalize column names (strip whitespace)
df.columns = pd.MultiIndex.from_tuples([(str(a).strip(), str(b).strip()) for a, b in df.columns])

# 1) find the column that holds the gene/feature IDs
def find_feature_column(columns):
    # prefer a first-level name containing 'feature' (case-insensitive)
    for col in columns:
        a, b = col
        if 'feature' in str(a).lower() or 'feature' in str(b).lower() or 'feature id' in str(a).lower():
            return col
    # fallback: second column (index 1) if nothing matched
    return columns[1]

feature_col = find_feature_column(df.columns)
print("Using feature column:", feature_col)

# 2) find TF presence columns (skip motif columns)
tf_cols = []
for col in df.columns:
    a, b = col
    a_low = str(a).lower()
    # include columns whose first-level header mentions 'tf' but not 'motif'
    if 'tf' in a_low and 'motif' not in a_low:
        # ensure TF name exists (b) and not something like 'TFs' placeholder
        if str(b).strip() and not str(b).strip().lower() in ('tfs', 'tf', 'tf motif'):
            tf_cols.append(col)

print(f"Found {len(tf_cols)} TF columns (presence columns). Sample:", tf_cols[:6])

# 3) collect edges (TF -> Feature) where the cell is non-zero (treat any non-zero as presence)
edges = []
for col in tf_cols:
    tf_name = str(col[1]).upper().strip()        # TF name like 'ABF1'
    # convert column to numeric if possible, treat NaN as 0
    series = pd.to_numeric(df[col].fillna(0), errors='coerce').fillna(0)
    mask = series != 0
    # get gene ids from feature_col
    genes = df.loc[mask, feature_col]
    for g in genes:
        edges.append((tf_name, str(g).strip()))

# 4) deduplicate and sort for reproducibility
edges_unique = sorted(set(edges))

print(f"Collected {len(edges)} raw edges, {len(edges_unique)} unique edges after dedup.")

ssTFs_SGDID = {}
infile = open('./data/All_gene_ID.txt', 'r')
for line in infile:
    ssTFs_SGDID[line.split()[1]] = line.split()[0]
infile.close()

converted_edges = []
for each_edge in edges_unique:
    if each_edge[1] in ssTFs_SGDID:
        converted_edges.append((each_edge[0], ssTFs_SGDID[each_edge[1]]))
    else:
        converted_edges.append((each_edge[0], each_edge[1]))


net = union_to_network_json(converted_edges)

# pretty-print to console
#print(json.dumps(net, indent=2))

# save to file
with open("./data/Rossi_TF_DNA_binding.json", "w") as f:
    json.dump(net, f, indent=2)
