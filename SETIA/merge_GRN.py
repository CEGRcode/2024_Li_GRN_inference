#!/usr/bin/env python3
"""
merge_GRN.py

Merge two JSON graph files (with "nodes" and "edges" lists).

Rules:
- Nodes: deduplicated by 'id'. If duplicate, later file overwrites earlier.
- Edges: if the same (source, target) exists in both → raise an error.
- Otherwise, edges are merged and written out cleanly.

Usage:
  python merge_graphs_strict.py file1.json file2.json -o merged.json
"""

import json
import sys
import argparse

def merge_graphs(file1, file2):
    # load files
    with open(file1, "r", encoding="utf-8") as f:
        g1 = json.load(f)
    with open(file2, "r", encoding="utf-8") as f:
        g2 = json.load(f)

    # --- merge nodes ---
    nodes = {}
    for node in g1.get("nodes", []):
        if "id" in node:
            nodes[node["id"]] = node
    for node in g2.get("nodes", []):
        if "id" in node:
            nodes[node["id"]] = node  # later file overwrites earlier

    # --- check and merge edges ---
    edges = {}
    for edge in g1.get("edges", []):
        k = (edge.get("source"), edge.get("target"))
        edges[k] = edge

    for edge in g2.get("edges", []):
        k = (edge.get("source"), edge.get("target"))
        if k in edges:
            raise ValueError(f"Duplicate edge found between {k[0]} → {k[1]}")
        edges[k] = edge

    # --- result ---
    merged = {
        "nodes": list(nodes.values()),
        "edges": list(edges.values())
    }

    return merged

def main():
    parser = argparse.ArgumentParser(description="Merge two JSON graphs with strict edge uniqueness.")
    parser.add_argument("file1", help="First JSON file")
    parser.add_argument("file2", help="Second JSON file")
    parser.add_argument("-o", "--output", required=True, help="Output JSON file")

    args = parser.parse_args()

    try:
        merged = merge_graphs(args.file1, args.file2)
    except ValueError as e:
        print(f"❌ Error: {e}", file=sys.stderr)
        sys.exit(1)

    with open(args.output, "w", encoding="utf-8") as out:
        json.dump(merged, out, indent=2, ensure_ascii=False)

    print(f"✅ Merged successfully! Nodes: {len(merged['nodes'])}, Edges: {len(merged['edges'])}")

if __name__ == "__main__":
    main()
