# Station Velocity Methods in MCMC

This document describes the three methods for computing station velocities from triangular dislocation element (TDE) slip rates during MCMC sampling, configured via `mcmc_station_velocity_method` in `config.py` (implemented in `celeri/solve_mcmc.py`).

## The Problem

In crustal deformation modeling, we need to compute how slip on fault surfaces produces velocities at GPS stations. This relationship is linear:

$$\mathbf{v} = -\mathbf{G} \, \mathbf{s}$$

where:
- $\mathbf{v} \in \mathbb{R}^{3n_{\text{station}}}$ — velocity at all stations in the internal model vector ordering (x, y, z)
- $\mathbf{s} \in \mathbb{R}^{n_{\text{TDE}}}$ — slip rate on each triangular element for a **single slip component** (strike-slip or dip-slip)
- $\mathbf{G} \in \mathbb{R}^{3n_{\text{station}} \times n_{\text{TDE}}}$ — Green's function matrix for that slip component

**The challenge**: $\mathbf{G}$ can be very large (thousands of TDEs × hundreds of stations), and MCMC sampling evaluates this computation many times (>100k for short single-chain runs).

## Method Comparison

| Aspect | `direct` | `low_rank` | `project_to_eigen` |
|--------|----------|------------|-------------------|
| **Computation** | $\mathbf{v} = -\mathbf{G}\mathbf{s}$ | $\mathbf{v} = -\mathbf{U}_k \boldsymbol{\Sigma}_k (\mathbf{V}_k^\top \mathbf{s})$ | $\mathbf{v} = \mathbf{H} (\boldsymbol{\Phi}^\top \mathbf{s})$ |
| **Basis** | None | SVD of $\mathbf{G}$ (most observable patterns) | Mesh Laplacian eigenvectors (smoothest patterns) |
| **Rank selection** | Full rank | Threshold-based ($\sigma > 10^{-5}$) | User-specified (`n_modes_*`) |
| **Memory during build** | $O(n_{\text{sta}} \cdot n_{\text{TDE,total}})$ | $O(n_{\text{sta}} \cdot n_{\text{TDE,total}})$ | $O(\max_i n_{\text{sta}} \cdot n_{\text{TDE},i})$† |
| **Memory during MCMC** | $O(n_{\text{sta}} \cdot n_{\text{TDE,total}})$ | $O(n_{\text{sta}} \cdot n_{\text{TDE,total}})$ | $O(n_{\text{sta}} \cdot m_{\text{total}})$‡ |
| **Streaming support** | No | No | Yes |

† With streaming enabled (`discard_tde_to_velocities=True`). Here $i$ indexes meshes, so this is the size of the largest single mesh's $\mathbf{G}_i$ matrix.
‡ Here $m_{\text{total}} = \sum_i (n_{\text{modes\_ss},i} + n_{\text{modes\_ds},i})$ is the total number of eigenmodes across all meshes. The dominant term is H matrix storage; Φ matrix storage is typically smaller.

## What is "Streaming"?

**Streaming** refers to processing meshes one at a time during operator building, rather than loading all Green's function matrices simultaneously.

**Without streaming**:
1. Load $\mathbf{G}_1, \mathbf{G}_2, \ldots, \mathbf{G}_n$ for all meshes into memory
2. Build derived operators
3. Keep all $\mathbf{G}_i$ in memory for MCMC

**With streaming** (used by default in `solve_mcmc()`, and only compatible with `project_to_eigen`):
1. Load $\mathbf{G}_1$ → compute $\mathbf{H}_1 = -\mathbf{G}_1 \boldsymbol{\Phi}_1$ → delete $\mathbf{G}_1$
2. Load $\mathbf{G}_2$ → compute $\mathbf{H}_2 = -\mathbf{G}_2 \boldsymbol{\Phi}_2$ → delete $\mathbf{G}_2$
3. ... and so on
4. During MCMC, only the smaller $\mathbf{H}_i$ matrices are in memory

This reduces peak memory from $O(\sum_i |\mathbf{G}_i|)$ to $O(\max_i |\mathbf{G}_i|)$ during build, and from $O(\sum_i |\mathbf{G}_i|)$ to $O(\sum_i |\mathbf{H}_i|)$ during MCMC.

## Method Details

### `direct`

**What it does**: Direct matrix-vector multiplication with the full Green's function.

```
v = -G @ s
```

**Advantages**:
- Exact computation, no approximation
- Simple to understand and debug

**Disadvantages**:
- Requires full $\mathbf{G}$ matrix in memory during MCMC
- Slowest per-sample computation for large problems
- No streaming mode available

**When to use**: Small problems, debugging, or when exact computation is required.

---

### `low_rank`

**What it does**: Approximates $\mathbf{G}$ via truncated SVD, keeping only singular values above a threshold.

$$\mathbf{G} \approx \mathbf{U}_k \boldsymbol{\Sigma}_k \mathbf{V}_k^\top$$

The velocity computation becomes:
$$\mathbf{v} = -\mathbf{U}_k (\boldsymbol{\Sigma}_k (\mathbf{V}_k^\top \mathbf{s}))$$

**Advantages**:
- Data-driven rank selection based on numerical significance
- Captures the most "observable" slip patterns (those producing largest surface signals)
- Improved numerical stability (filters near-null-space components)

**Disadvantages**:
- **Current implementation** requires full $\mathbf{G}$ matrix in memory at build time (SVD is computed during PyMC model construction, not during operator building)
- SVD recomputed for each mesh/slip-type combination at model build time
- **Streaming is fundamentally incompatible**: Unlike `project_to_eigen`, where each mesh has its own independent Laplacian eigenvectors, an SVD of the full operator $\mathbf{G} = [\mathbf{G}_1 \ \mathbf{G}_2 \ \cdots]$ requires global orthogonality of the left singular vectors $\mathbf{U}_k$. If you instead computed SVDs of each $\mathbf{G}_i$ independently (compute-and-discard), the resulting $\mathbf{U}_{k,i}$ would *not* be mutually orthogonal across meshes, so concatenating them would not yield a valid SVD of the whole system
- Threshold is currently hardcoded (`1e-5`)

**When to use**: When you want data-driven dimensionality reduction and memory is not a constraint.

---

### `project_to_eigen` (Default)

**What it does**: Projects the slip field onto the mesh eigenmode basis, then uses a pre-computed eigen-to-velocity operator.

$$\boldsymbol{\alpha} = \boldsymbol{\Phi}^\top \mathbf{s}, \quad \mathbf{v} = \mathbf{H} \boldsymbol{\alpha}$$

where:
- $\boldsymbol{\Phi}$ — eigenvectors of mesh Laplacian (smooth spatial functions)
- $\mathbf{H} = -\mathbf{G}\boldsymbol{\Phi}$ — pre-computed eigen-to-velocity operator

**Advantages**:
- **Memory efficient**: Only stores $\boldsymbol{\Phi}$ and $\mathbf{H}$, not full $\mathbf{G}$
- **Streaming mode**: Can process one mesh at a time, never holding all $\mathbf{G}$ matrices simultaneously
- **Explicit control**: User specifies number of modes via `n_modes_strike_slip` and `n_modes_dip_slip`
- **Consistent with parameterization**: When slip is parameterized as eigenmodes without constraints, this is mathematically exact
- **Physics-based basis**: Eigenmodes represent smooth spatial functions, providing implicit regularization

**Disadvantages**:
- For coupling or constrained elastic slip, the projection is an approximation (see below)
- Eigenmodes may not capture all patterns that produce observable surface velocities

**When to use**: Default choice for most problems, especially with large meshes or memory constraints.

## Mathematical Relationship

Both `low_rank` and `project_to_eigen` are low-rank approximations of the same linear operator, but with different bases:

| Method | Approximation | Optimality criterion |
|--------|--------------|---------------------|
| `low_rank` | $\mathbf{G} \approx \mathbf{U}_k \boldsymbol{\Sigma}_k \mathbf{V}_k^\top$ | Minimizes $\|\mathbf{G} - \tilde{\mathbf{G}}\|_F$ (Eckart-Young theorem) |
| `project_to_eigen` | $\mathbf{G}\mathbf{s} \approx \mathbf{H}(\boldsymbol{\Phi}^\top \mathbf{s})$ | Exact for slip in eigenmode basis; approximates smoothest patterns |

## Interaction with Slip Parameterization

### Elastic Component Without Constraints (`_elastic_component`, no bounds)

Slip is parameterized directly as eigenmode coefficients:
$$\mathbf{s} = \boldsymbol{\Phi} \boldsymbol{\alpha}$$

For this case, `project_to_eigen` is **mathematically exact**:
$$\mathbf{v} = -\mathbf{G}\boldsymbol{\Phi}\boldsymbol{\alpha} = \mathbf{H}\boldsymbol{\alpha}$$

The code uses the eigen-to-velocity operator directly (i.e., $\mathbf{v}=\mathbf{H}\boldsymbol{\alpha}$) without going through `_station_vel_from_elastic_mesh`.

### Elastic Component With Constraints (`_elastic_component`, with bounds)

When elastic rate bounds are applied, the slip field undergoes a nonlinear transformation:
$$\mathbf{s} = f(\boldsymbol{\Phi} \boldsymbol{\alpha})$$

where $f$ is sigmoid (for two-sided bounds) or softplus (for one-sided bounds).

The transformed slip $\mathbf{s}$ is **no longer in the eigenmode subspace**. The `project_to_eigen` method projects it back:
$$\boldsymbol{\alpha}' = \boldsymbol{\Phi}^\top \mathbf{s}, \quad \mathbf{v} = \mathbf{H} \boldsymbol{\alpha}'$$

This **discards components outside the eigenmode basis** and is an approximation.

### Coupling Component (`_coupling_component`)

Slip is the element-wise product of kinematic slip and a coupling field:
$$\mathbf{s} = \mathbf{k} \odot \mathbf{c}$$

where:
- $\mathbf{k}$ — kinematic slip rate from block rotations (full resolution, $n_{\text{TDE}}$ components)
- $\mathbf{c}$ — coupling field, parameterized as eigenmodes (possibly with sigmoid/softplus constraints)
- $\odot$ — Hadamard (element-wise) product: $(\mathbf{k} \odot \mathbf{c})_i = k_i \cdot c_i$

**Even without coupling constraints**, the element-wise product $\mathbf{k} \odot \mathbf{c}$ is generally **not** in the eigenmode subspace, because multiplying a general vector by an eigenmode expansion does not yield another eigenmode expansion.

The `project_to_eigen` method projects this product back onto the eigenmode basis, which is an approximation. The approximation quality depends on how well the product can be represented by the available eigenmodes—typically good when both $\mathbf{k}$ and $\mathbf{c}$ are spatially smooth.

### Summary Table

| Case | Slip formula | `project_to_eigen` |
|------|-------------|-------------------|
| Elastic, no constraints | $\mathbf{s} = \boldsymbol{\Phi}\boldsymbol{\alpha}$ | **Exact** |
| Elastic, with constraints | $\mathbf{s} = f(\boldsymbol{\Phi}\boldsymbol{\alpha})$ | Approximate |
| Coupling, no constraints | $\mathbf{s} = \mathbf{k} \odot (\boldsymbol{\Phi}\boldsymbol{\beta})$ | Approximate |
| Coupling, with constraints | $\mathbf{s} = \mathbf{k} \odot f(\boldsymbol{\Phi}\boldsymbol{\beta})$ | Approximate |

## Memory Behavior

### Without Streaming (`discard_tde_to_velocities=False`)

All $\mathbf{G}$ matrices are stored in memory:
- Peak memory: $O(n_{\text{meshes}} \cdot n_{\text{sta}} \cdot n_{\text{TDE}})$

However, **only `direct` and `low_rank` actually use the $\mathbf{G}$ matrices during MCMC**. The `project_to_eigen` method only uses the pre-computed $\mathbf{H}$ matrices—the $\mathbf{G}$ matrices are stored but sit idle, wasting memory. This is why streaming is enabled by default in `solve_mcmc()` when using `project_to_eigen`.

### With Streaming (`discard_tde_to_velocities=True`, only for `project_to_eigen`)

Process one mesh at a time during operator building:
1. Load $\mathbf{G}_i$
2. Compute $\mathbf{H}_i = -\mathbf{G}_i \boldsymbol{\Phi}_i$
3. Discard $\mathbf{G}_i$

Memory usage:
- Peak during build: $O(\max_i |\mathbf{G}_i|)$
- During MCMC: $O(\sum_i (n_{\text{sta}} + n_{\text{TDE}}) \cdot m_i)$

## Real-World Example: Western North America (WNA) Large Model

Analysis of the `large_01` branch in `wna/` provides concrete numbers for a production-scale model:

### Model Size

| Component | Count |
|-----------|-------|
| Stations | 17,146 |
| Meshes | 198 |
| Total TDEs | 182,831 |
| Total eigenmodes (SS + DS) | 4,000 |

### Per-Mesh Statistics

| Metric | Value |
|--------|-------|
| Min TDEs per mesh | 40 |
| Max TDEs per mesh | 6,679 |
| Avg TDEs per mesh | 923 |
| Modes per mesh | 10 SS + 10 DS (typical) |

### Memory Comparison

| Method | Matrix Size | Memory (float64) |
|--------|-------------|------------------|
| `direct` / `low_rank` | $\sum_i \mathbf{G}_i$: 51,438 × 365,662 | **150.5 GB** |
| `project_to_eigen` | $\sum_i \mathbf{H}_i$: 34,292 × 4,000 | **1.1 GB** |

**Memory reduction: 137×**

This demonstrates why streaming with `project_to_eigen` is essential for large-scale models—the full Green's function matrices simply don't fit in memory on typical hardware.

### Streaming Benefit

With streaming enabled, peak memory during operator building is determined by the largest single mesh:
- Largest $\mathbf{G}_i$: 51,438 × 13,358 ≈ **5.5 GB**

Without streaming, all 198 meshes would need to be loaded simultaneously (150+ GB).

## Configuration

```json
{
  "mcmc_station_velocity_method": "project_to_eigen"
}
```

Options: `"direct"`, `"low_rank"`, `"project_to_eigen"` (default)

## Historical Development

1. **Oct 2025**: Original MCMC solver used direct multiplication
2. **Oct 2025**: `low_rank` added for speed (commit `75950f0`)
3. **Oct 2025**: `project_to_eigen` added and made default (commit `d09f223`)
4. **Jan 2026**: Streaming mode added for `project_to_eigen` (commit `2be0773`)

The memory optimization via streaming was added after the methods were created. The original motivation for both `low_rank` and `project_to_eigen` was computational speed, not memory efficiency.

## Potential Improvements

1. **Truncated SVD**: Use `scipy.sparse.linalg.svds` instead of full SVD for efficiency when only top $k$ singular vectors are needed.

2. **Configurable threshold**: Make the `low_rank` singular value threshold configurable instead of hardcoded at `1e-5`.
