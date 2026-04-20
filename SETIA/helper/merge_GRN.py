#!/usr/bin/env python3
import json
import glob
import sys

input_pattern = './result/GRN_nonTF_*.json'
output_file = './result/GRN_nonTF_merged.json'

files = sorted(glob.glob(input_pattern))
print(f'Found {len(files)} files to merge.')

nodes = {}
edges = {}

for fpath in files:
    with open(fpath, 'r', encoding='utf-8') as f:
        g = json.load(f)
    for node in g.get('nodes', []):
        if 'id' in node:
            nodes[node['id']] = node
    for edge in g.get('edges', []):
        k = (edge.get('source'), edge.get('target'))
        if k in edges:
            print(f'⚠️  Duplicate edge: {k[0]} → {k[1]} (in {fpath}), skipping.')
            continue
        edges[k] = edge

merged = {
    'nodes': list(nodes.values()),
    'edges': list(edges.values())
}

with open(output_file, 'w', encoding='utf-8') as f:
    json.dump(merged, f, indent=2, ensure_ascii=False)

print(f'✅ Done! Nodes: {len(merged["nodes"])}, Edges: {len(merged["edges"])}')
