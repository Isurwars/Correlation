# ADR 0003: Machine Learning Interatomic Potentials (MLIP) & GNN Descriptors

## Status
Accepted

## Context
Modern atomic simulation workflows increasingly rely on Machine Learning Interatomic Potentials (e.g., MACE, CHGNet, ORB, SevenNet) and Graph Neural Networks (GNNs). To evaluate structural descriptors and local density of states (LDOS/TDOS) directly from ML models, Correlation requires native graph construction under periodic boundary conditions (PBC) and seamless tensor exchange with deep learning frameworks.

## Decision
1. **Periodic Graph Construction**: Implement `PeriodicGraphBuilder` in C++ with minimum image convention (MIC) handling to generate node feature matrices, edge indices, periodic displacement vectors, and radial Bessel basis projections directly from atomic structures.
2. **PyTorch Geometric (PyG) Integration**: Expose `to_torch_geometric()` in Python (`python/correlation/`) for zero-copy transfer of atomic graphs to `torch_geometric.data.Data`.
3. **Optional C++ LibTorch Inference**: Provide optional LibTorch linkage (`CORRELATION_ENABLE_LIBTORCH=ON`) allowing headless C++ inference of TorchScript-exported GNN models without requiring Python at runtime.
4. **TDOS Accumulation**: Implement `TDOSCalculator` to aggregate per-atom local densities of states predicted by MLIPs into global total densities of states.

## Consequences
### Positive
- Direct interoperability with state-of-the-art MLIP architectures (MACE, SevenNet, CHGNet, MatterSim).
- High-speed periodic graph extraction in C++ outperforming Python-only neighbor list builders.

### Negative / Trade-offs
- LibTorch dependencies are large (~1–2 GB) when enabled; therefore, LibTorch remains strictly optional (`CORRELATION_ENABLE_LIBTORCH=OFF` by default). Python PyG adapters provide zero-overhead integration for typical user workflows.
