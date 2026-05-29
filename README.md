# Overview

This repository contains the implementation of the paper: **"Graph Reordering via Traversal Orders for Accelerating Graph Traversal Algorithms"**.
Given a graph in CSR format, we first generate traversal order for each vertex, compress them with k-cmer and MinHash sketches, group vertices
with a hierarchical reordering framework, and output a global vertex ordering.

## Requirements
- Linux server with GCC
- CMake
- OpenMP
- Linux `perf_event_open` support for hardware counter measurements

## Build
```bash
cd OCRH
make
```

## Input Format
Both tools read a binary CSR graph.

## Generate the Order
```bash
cd OCRH
./OCRH <graph.csr> <sampling> <k_clusters> <traversal_type> [-large]
```

Arguments:

| Argument | Description |
| --- | --- |
| `graph.csr` | Input graph in the binary CSR format above |
| `topk` | Number of high-degree vertices selected as traversal sources |
| `k_clusters` | Number of clusters used for vertex grouping |
| `traversal_type` | Signature traversal type: `0` = BFS, `1` = DFS, `2` = SSSP |
| `-large` | Optional memory-saving mode for large graphs |

Example:

```bash
./OCRH ../data/com-amazon.ungraph.csr 256 256 0
```

## Evaluate the Order
```bash
cd traversal
./eval <graph.csr> <order_dir> <traversal_num> [repeat_count]
```

Arguments:

| Argument | Description |
| --- | --- |
| `graph.csr` | Input graph in binary CSR format |
| `order_dir` | Directory containing one or more order files |
| `traversal_num` | Workload: `1` = BFS, `2` = DFS, `3` = SSSP, `4` = WCC, `5` = PageRank, `6` = Diameter |
| `repeat_count` | Optional number of repeated measurements; default is `1` |

Example:

```bash
./eval ./data/amazon.csr ../orders 0 1 5
```

For each workload, the evaluator compares:

- the original vertex order;
- a degree-sorted order;
- every order file found in `order_dir`;
- a workload-specific self order used as an upper-reference baseline.

## Dataset Notes

SNAP graphs can be used for quick experiments. For example:

```bash
wget http://snap.stanford.edu/data/com-amazon.ungraph.txt.gz
gunzip com-amazon.ungraph.txt.gz
```

This repository does not include a dataset converter. Before running the tools,
convert the edge list into the required binary CSR layout.

## Contact

For any questions or inquiries, please contact us at msd673@hnu.edu.cn.
