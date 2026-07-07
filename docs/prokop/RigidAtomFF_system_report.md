# RigidAtomFF System Report

## Overview

Two implementations of rigid-body molecular dynamics with port-based bond constraints:

1. **RRsp3** — GPU (OpenCL) prototype in Python, located in `py/RRsp3/`
2. **RigidAtomFF** — CPU (pure C) port in Box3D, located in `src/rigid_atom_ff.{c,h}`

Both solve the same problem: simulate molecules as rigid atoms connected by spring-loaded "port" constraints, with collision handling and optional soft forces.

---

## 1. RRsp3 — GPU Prototype (`py/RRsp3/`)

### Architecture

OpenCL-accelerated XPBD/Jacobi solver with cluster-sorted topology and ghost atoms for workgroup-local memory access.

**Key files:**

- `@/home/prokop/git/box3d/py/RRsp3/RRsp3.py` — Python wrapper class `RRsp3`, manages OpenCL buffers, kernel dispatch, upload/download
- `@/home/prokop/git/box3d/py/RRsp3/RRsp3.cl` — 1748-line OpenCL kernel: bboxes, local topology, collisions, port constraints, corrections
- `@/home/prokop/git/box3d/py/RRsp3/RRsp3_utils.py` — molecule builders, packing, momentum analysis, debug log parsing
- `@/home/prokop/git/box3d/py/RRsp3/XPTB_utils.py` — molecule packing (`pack_molecules_contiguous`), geometry helpers, XYZ I/O

### Data Model

- Atoms stored as `float4` (xyz + invMass in `.w`)
- Grouped into clusters of `group_size=64` (one molecule per group, padded with dummy atoms)
- **Ghost atoms**: atoms from neighboring groups copied into local index space for workgroup-local collision/port sweeps
- **Node/cap split**: atoms with degree > 1 are "nodes" (have rotational inertia), degree ≤ 1 are "caps" (massless, follow node rotation)

### Solver Pipeline (per step)

```
1. update_bboxes_rigid(pos, radius) → per-group AABBs
2. build_local_topology_rigid(neighs, excl, bboxes) → local indices, ghost lists
3. compute_collision_cluster_rigid(pos, radius, excl_local) → dpos_coll
4. compute_ports_cluster_rigid(pos, quat, port_local) → dpos_node, drot_node, dpos_neigh
5. apply_corrections_rigid_ports(dpos_coll, dpos_node, dpos_neigh) → updated pos, quat
```

### Port Constraint Variants

| Kernel | Rotation model | Description |
|--------|---------------|-------------|
| `current` | Massfull (physical inertia) | Backward-Euler impulse on port tip displacement |
| `orig` | Massfull | Original formulation, slightly different damping |
| `substep_optimized` | Massless (geometric) | Newton iteration alignment, zero rotational inertia |
| `shapematch` | Massless | Kabsch/shape-matching optimal rotation |
| `eigen` | Massless | Two-pass: compute tips → eigen-decomposition optimal rotation |

### Momentum Acceleration

Heavy-ball momentum buffers (`dpos_mom`, `dquat_mom`) accelerate Jacobi convergence within a timestep. Reset between timesteps via `reset_momentum()`.

### Test Results (2026-07-07)

| Test | Status | Notes |
|------|--------|-------|
| `test_RRsp3_smoke.py` | PASS | Ghost discovery, local indexing, single step |
| `test_RRsp3_momentum.py --steps 5` | PASS | `\|dP\| < 3e-8`, `\|dL\| < 1e-8` |
| `test_RRsp3_debug.py` | PASS | Collision/exclusion logging, inter-molecular COLLIDE detected |
| `test_RRsp3_headless.py smoke` | PASS | No NaN, range checks OK |
| `test_RRsp3_headless.py momentum` | PASS | 50 steps within tolerance |
| `test_RRsp3_headless.py suite L0-L2` | PASS | All 3 levels pass |
| `test_RRsp3_headless.py convergence` | PASS | Converges at step 18 with `--stretch 1.05 --tol 1e-3` |
| `test_RRsp3_headless.py topology` | FAIL | Pre-existing: kernel debug print gating misses some pairs |
| `test_RRsp3_vispy.py` | PASS | GUI launches, renders, mouse interaction works (PyQt5 + vispy) |

### Dependencies

- `pyopencl` (2023.1.3), `numpy` (1.26.4), `vispy` (0.14.1), `PyQt5`
- GPU: NVIDIA GeForce RTX 3090

---

## 2. RigidAtomFF — CPU Port (`src/rigid_atom_ff.{c,h}`)

### Architecture

Pure C implementation using Sequential Impulses (SI) with local implicit (backward Euler) harmonic impulses. No OpenCL, no ghost atoms, no workgroup complexity — all atoms directly accessible.

**Key files:**

- `@/home/prokop/git/box3d/src/rigid_atom_ff.h` — data structures and API
- `@/home/prokop/git/box3d/src/rigid_atom_ff.c` — solver implementation (471 lines)
- `@/home/prokop/git/box3d/test/test_rigid_atom_ff.c` — unit tests (337 lines)
- `@/home/prokop/git/box3d/samples/sample_rigid_atoms.cpp` — interactive visualization (418 lines)

### Data Structures

```c
typedef struct {
    b3Vec3 pos, vel, omega;
    b3Quat quat;
    float mass, invMass, radius;
    b3Matrix3 invInertia;
    int nPorts;
} RAffAtom;

typedef struct {
    int neighIdx;
    b3Vec3 localDir;    // lever arm in body frame
    float kPort, dPort; // per-port stiffness/damping override (0 = world default)
} RAffPort;

typedef struct {
    RAffAtom* atoms;    // nAtoms
    RAffPort* ports;    // nAtoms * RAFF_MAX_PORTS (flat)
    int* incomingCount, *incomingOffset;
    int* incomingAtoms, *incomingPorts;  // reverse adjacency
    RAffCollisionPair* collisionPairs;   // dynamic broad-phase
    b3Vec3 *fsoft, *tsoft;               // cached soft forces
    b3Vec3 *Jown, *Lown, *Jrecoil, *Jcoll; // impulse accumulators
    float H, kBond, dBond, kColl, dColl, relaxation, softCutoff;
    int nsub;
} RAffWorld;
```

### Time Stepping (RESPA / Strang Splitting)

```
RAffWorld_Step(w):
    applySoftKick(w, 0.5*H, fsoft, tsoft)    // first half soft kick
    h = H / nsub
    for s in 0..nsub:
        hardMicrostep(w, h)                   // port impulses + drift
    evaluateSoftForces(w)                     // O(N²) Lennard-Jones
    applySoftKick(w, 0.5*H, fsoft, tsoft)    // second half soft kick
```

### Port Impulse (Local Implicit / Backward Euler)

For each directed port k of atom i toward neighbor j:

```
r = rotate(quat_i, localDir_k)           // world-frame lever arm
e = pos_j - (pos_i + r)                  // constraint error
C = |e|;  n = e/C                        // normal
u = n · (v_j - (v_i + ω_i × r))         // constraint velocity
s = r × n                                // torque arm
w = invM_i + invM_j + sᵀ invI_i s       // effective inverse mass

K = kBond/2  (per directed port; reciprocal ports double energy)
D = dBond/2

p = (hK·C + h·(D+hK)·u) / (1 + (hD+h²K)·w)
P = p·n

vel_i  += invM_i · P
omega_i += invI_i · (r × P)
vel_j  -= invM_j · P
```

### Collision Impulse

```
delta = R_ij - |pos_i - pos_j|    // penetration depth
if delta <= 0: skip
n = (pos_i - pos_j) / |...|
u_c = n · (v_i - v_j)
ω = sqrt(K · w_eff)               // natural frequency
ζ = D·w_eff / (2ω)               // damping ratio
a2 = h·ω·(2ζ + h·ω)
p = -(1/w_eff) · (a2/(1+a2)·u_c + ω/(2ζ+h·ω)·delta)
if p < 0: p = 0                   // unilateral clamp
vel_i += invM_i · p·n
vel_j -= invM_j · p·n
```

### Molecule Loading

`RAffWorld_LoadMol()` loads from `b3MolFile` (MDL MOL format):
- Atom positions, covalent radii (`b3GetCovalentRadius`), atomic masses (`b3GetAtomicMass`)
- Bonds → reciprocal ports with body-frame lever arms
- Spherical inertia: `I = 2/5 · m · r²`

### Visualization Sample

`@/home/prokop/git/box3d/samples/sample_rigid_atoms.cpp`:
- Loads `data/mol/nNonan.mol` (nonane C9H20, 29 atoms)
- Gravity + ground collision (position correction + velocity bounce)
- Ctrl+Click mouse grabbing via ray-sphere intersection
- ImGui controls for H, nsub, kBond, dBond, gravity, restitution, friction
- ImPlot KE history graph (512-frame ring buffer)
- Element coloring: H = white, C/other = blue, grabbed = orange, selected = green
- NaN detection with terminal reporting

### Test Results

| Test | Status | Notes |
|------|--------|-------|
| `RigidAtomFF_BondTest` | PASS | 2-atom bond converges in 50 steps, momentum conserved |
| `RigidAtomFF_CollisionTest` | DISABLED | Collisions commented out in `hardMicrostep` |
| `RigidAtomFF_MomentumTest` | PASS | 4-atom chain, 100 steps, `\|ΔP\| < 0.5` |
| `RigidAtomFF_ChainTest` | PASS | 6-atom chain, 200 steps, no explosion, bond lengths OK |
| Sample (nNonan) | RUNS | 1823 frames, 0 sokol errors, GPU: RTX 3090 |

### Current Limitations

- **Collisions disabled** in `hardMicrostep` (`@/home/prokop/git/box3d/src/rigid_atom_ff.c:345`) — port solver stability first
- **O(N²) broad phase** — no spatial hash yet
- **Spherical inertia only** — no anisotropic inertia tensors
- **No angle constraints** — only bond (port) constraints implemented
- **Single-threaded** — no parallel port sweeps or graph coloring

---

## 3. Comparison: RRsp3 (GPU) vs RigidAtomFF (CPU)

| Aspect | RRsp3 (GPU) | RigidAtomFF (CPU) |
|--------|-------------|-------------------|
| Language | Python + OpenCL C | Pure C |
| Parallelism | GPU workgroups (64 atoms/group) | Single-threaded |
| Ghost atoms | Yes (workgroup-local memory) | Not needed (direct access) |
| Topology | Cluster-sorted + local indexing | Flat global indexing |
| Port solver | Jacobi (5 variants) | Sequential Impulses (1 variant) |
| Rotation model | Massfull or massless (selectable) | Massfull only |
| Momentum accel | Heavy-ball | None |
| Collisions | Jacobi positional correction | SI with damping ratio (disabled) |
| Soft forces | Not in kernel (external) | Lennard-Jones O(N²) |
| Molecule loading | XYZ files (via `load_xyz`) | MOL files (via `b3LoadMolFile`) |
| Visualization | vispy + PyQt5 (3D scene) | Box3D samples (spheres + ImGui) |
| Test coverage | 7 test scripts, multiple systems | 4 unit tests |

---

## 4. Design Documents

- `@/home/prokop/git/box3d/docs/prokop/RigidAtomFF_SI_short.plan.md` — concise implementation plan (175 lines)
- `@/home/prokop/git/box3d/docs/prokop/RigidAtomFF_SI_long.plan.md` — detailed plan with code sketches (449 lines)
- `@/home/prokop/git/box3d/docs/prokop/RigidAtomFF_SI.chat.md` — full design discussion (4523 lines)
- `@/home/prokop/git/box3d/py/RRsp3/docs/` — 13 design/progress/chat docs for the GPU prototype

---

## 5. Key Physics

### Port Constraints

Each bond is modeled as two reciprocal "ports" — body-frame lever arms that rotate with the atom. The constraint enforces that the rotated port tip coincides with the neighbor atom's center. This gives:

- **Bond length** enforcement (radial component of port error)
- **Bond angle** enforcement (tangential component → torque on atom)
- **Momentum conservation**: equal-and-opposite impulses along the constraint normal; torque arm `r × P` is parallel to `e`, so angular momentum is conserved

### Parameters

| Parameter | RRsp3 | RigidAtomFF | Meaning |
|-----------|-------|-------------|---------|
| K (stiffness) | `K=200` (dimensionless) | `kBond=1e6` (physical) | Spring constant for port correction |
| D (damping) | implicit in relaxation | `dBond=1e4` | Radial velocity damping |
| H (macrostep) | `dt=0.1` | `H=0.02` | Soft force timestep |
| nsub | N/A (Jacobi iterations) | `nsub=16` | Hard microsteps per macrostep |
| relaxation | `0.5` (Jacobi) | `0.5` (unused in SI) | Damping for multi-constraint atoms |

### Momentum Conservation Proof

**Ports**: `P_i = +p·n` (at tip), `P_j = -p·n` (at center). Linear: `P_i + P_j = 0`. Angular: `(tip - x_j) × P = -e × (p·n) = 0` since `P ∥ e`.

**Collisions**: `P_i = +p·n`, `P_j = -p·n`. Linear: `P_i + P_j = 0`. Angular: sphere-sphere has no torque (force through center).

---

## 6. File Map

```
py/RRsp3/
├── RRsp3.py              OpenCL wrapper class (975 lines)
├── RRsp3.cl              OpenCL kernels (1748 lines)
├── RRsp3_utils.py        Molecule builders, packing, analysis (923 lines)
├── XPTB_utils.py         Packing, geometry, XYZ I/O (created during porting)
├── test_RRsp3_smoke.py   Basic smoke test
├── test_RRsp3_debug.py   Collision/exclusion debug logging
├── test_RRsp3_momentum.py  Momentum conservation
├── test_RRsp3_headless.py  Consolidated test runner (822 lines)
├── test_RRsp3_vispy.py     Interactive GUI debugger (861 lines)
├── test_RRsp3_convergence.py  Convergence testing
└── docs/                 13 design/progress documents

src/
├── rigid_atom_ff.h       API + data structures (114 lines)
├── rigid_atom_ff.c       Solver implementation (471 lines)
└── mol_loader.h          MOL file parser (included by rigid_atom_ff.h)

test/
└── test_rigid_atom_ff.c  Unit tests (337 lines)

samples/
└── sample_rigid_atoms.cpp  Interactive visualization (418 lines)

data/mol/
└── nNonan.mol            Nonane molecule (C9H20)
```
