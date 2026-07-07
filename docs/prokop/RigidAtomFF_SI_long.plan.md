# RigidAtomFF SI Implementation Plan
 
Concrete plan for reimplementing RRsp3 rigid-atom dynamics in pure C within Box3D,
using soft-step local implicit Sequential Impulses in an atom-centered gather-apply architecture.
 
---
 
## 1. Architecture Overview
 
Three-level time-stepping structure:
 
```
soft macrostep H
  -> evaluate expensive soft forces O(N^2) once
  -> s hard microsteps (h = H/s):
       frozen soft kick h
       gather all port + collision impulses on each atom
       apply one velocity/rotation update per atom
       drift by h
  -> new soft half-kick
```
 
Key principles:
- **Atom-centered**: only dynamic state is `pos[i]`, `quat[i]`, `vel[i]`, `omega[i]`. All forces/torques computed on the fly.
- **Local implicit harmonic impulse** for bonds, angles/ports, and collisions — no global linear solve, no exact constraints.
- **No warm starting**: no persistent Lagrange multipliers between timesteps.
- **Scalar radial reactions** along constraint normals to guarantee angular momentum conservation.
- **Momentum conservation**: equal-and-opposite impulses; ports produce torque on central atom only; collisions produce no torque (sphere-sphere).
 
---
 
## 2. Data Structures
 
```c
// Per-atom dynamic state (the only persistent arrays)
typedef struct {
    b3Vec3 pos;       // position of atom center
    b3Quat quat;      // orientation (normalized)
    b3Vec3 vel;       // linear velocity
    b3Vec3 omega;     // angular velocity
    float   mass;     // atomic mass
    float   invMass;  // 1/mass
    b3Matrix3 invInertia; // inverse inertia tensor (world frame)
    float   radius;   // van der Waals radius for collisions
    int     nPorts;   // number of outgoing ports (typically <=4)
} RigidAtom;
 
// Per-port topology (static, computed once from molecular graph)
typedef struct {
    int   neighIdx;   // index of neighbor atom this port points to
    b3Vec3 localDir;  // port direction in atom's local frame
    float  stiffness; // algorithmic stiffness parameter a = h^2 * K * w
    float  damping;   // algorithmic damping parameter d = h * C * w
} RigidPort;
 
// Per-atom port list (flat array for cache efficiency)
// ports[i * MAX_PORTS + k] gives port k of atom i
#define MAX_PORTS 4
 
// Collision pair (rebuilt each microstep or every few microsteps)
typedef struct {
    int   i, j;       // atom indices
    float Rij;        // sum of radii
} CollisionPair;
 
// Simulation context
typedef struct {
    RigidAtom*     atoms;
    int            nAtoms;
    RigidPort*     ports;       // flat: nAtoms * MAX_PORTS
    int*           portCount;   // per-atom port count
 
    // Reverse adjacency: for each atom, which (sourceAtom, sourcePort) pairs target it
    int*           incomingCount;
    int*           incomingAtoms;   // flat: total incoming refs
    int*           incomingPorts;   // flat: total incoming refs
 
    // Collision pairs (dynamic, from broad phase)
    CollisionPair* collisionPairs;
    int            nCollisionPairs;
 
    // Soft forces (cached, evaluated once per macrostep)
    b3Vec3*        fsoft;    // force per atom
    b3Vec3*        tsoft;    // torque per atom
 
    // Temporary accumulators (per microstep)
    b3Vec3*        Jown;     // own port linear impulse
    b3Vec3*        Lown;     // own port angular impulse
    b3Vec3*        Jrecoil;  // port recoil sent to neighbor (flat: nAtoms * MAX_PORTS)
    b3Vec3*        Jcoll;    // collision linear impulse
 
    // Parameters
    float H;             // soft macrostep
    int   nsub;          // number of hard microsteps
    float aBond;         // correction fraction for bonds (dimensionless)
    float dBond;         // damping for bonds (dimensionless)
    float aAngle;        // correction fraction for angles/ports
    float dAngle;        // damping for angles/ports
    float aCollision;    // correction fraction for collisions
    float dCollision;    // damping for collisions
} RigidAtomWorld;
```
 
### Box3D math types used
 
From `@/home/prokop/git/box3d/include/box3d/math_functions.h`:
- `b3Vec3`, `b3Quat`, `b3Matrix3` — core types
- `b3RotateVector(q, v)` — rotate vector by quaternion
- `b3MulQuat`, `b3NormalizeQuat`, `b3Conjugate` — quaternion algebra
- `b3MakeQuatFromAxisAngle(axis, angle)` — construct rotation
- `b3Cross`, `b3Dot`, `b3Length`, `b3Normalize` — vector ops
- `b3Solve3(M, b)` — solve 3x3 linear system
- `b3InvertMatrix(M)` — 3x3 inverse
 
From `@/home/prokop/git/box3d/src/math_internal.h`:
- `b3IntegrateRotation(q, deltaRotation)` — integrate quaternion from angular velocity * dt
- `b3QuatFromExponentialMap(v)` — exponential map quaternion
- `b3DeltaQuatToRotation(q, target)` — get angular displacement from quaternion difference
 
---
 
## 3. Core Algorithm
 
### 3.1 Macrostep (RESPA / Strang splitting)
 
```c
void RigidAtomWorld_Step(RigidAtomWorld* w) {
    // Half soft kick at beginning (using cached force from end of previous step)
    applySoftKick(w, 0.5f * w->H, w->fsoft, w->tsoft);
 
    float h = w->H / w->nsub;
 
    for (int isub = 0; isub < w->nsub; isub++) {
        // Distribute frozen soft force over microsteps
        for (int i = 0; i < w->nAtoms; i++) {
            w->atoms[i].vel   = b3Add(w->atoms[i].vel,   b3MulSV(h * w->atoms[i].invMass, w->fsoft[i]));
            w->atoms[i].omega = b3Add(w->atoms[i].omega, b3MulMV(w->atoms[i].invInertia, b3MulSV(h, w->tsoft[i])));
        }
 
        // Hard microstep: gather + apply
        hardMicrostep(w, h);
    }
 
    // Evaluate new soft forces at final configuration
    evaluateSoftForces(w, w->fsoft, w->tsoft);  // O(N^2) — Coulomb, vdW, torsion
 
    // Half soft kick at end with new force
    applySoftKick(w, 0.5f * w->H, w->fsoft, w->tsoft);
}
```
 
### 3.2 Hard Microstep
 
```c
void hardMicrostep(RigidAtomWorld* w, float h) {
    // 1. Clear accumulators
    for (int i = 0; i < w->nAtoms; i++) {
        w->Jown[i]  = b3Vec3_zero;
        w->Lown[i]  = b3Vec3_zero;
        w->Jcoll[i] = b3Vec3_zero;
    }
    for (int i = 0; i < w->nAtoms * MAX_PORTS; i++) {
        w->Jrecoil[i] = b3Vec3_zero;
    }
 
    // 2. Update collision pairs (broad phase — can be done every few microsteps)
    updateCollisionPairs(w);
 
    // 3. Gather port impulses (one pass over all atoms)
    gatherPortImpulses(w, h);
 
    // 4. Gather collision impulses (one pass over all atoms or all pairs)
    gatherCollisionImpulses(w, h);
 
    // 5. Apply: aggregate J_i = Jown + Jrecoil_incoming + Jcoll, then update
    applyAtomUpdate(w, h);
}
```
 
### 3.3 Port Impulse (Local Implicit Harmonic)
 
For port `k` of atom `i` pointing toward neighbor `j`:
 
```c
void gatherPortImpulses(RigidAtomWorld* w, float h) {
    for (int i = 0; i < w->nAtoms; i++) {
        RigidAtom* ai = &w->atoms[i];
        b3Vec3 Jacc = b3Vec3_zero;
        b3Vec3 Lacc = b3Vec3_zero;
 
        for (int k = 0; k < ai->nPorts; k++) {
            RigidPort* port = &w->ports[i * MAX_PORTS + k];
            int j = port->neighIdx;
            if (j < 0) continue;
 
            RigidAtom* aj = &w->atoms[j];
 
            // Port tip position in world
            b3Vec3 r = b3RotateVector(ai->quat, port->localDir);
            b3Vec3 tip = b3Add(ai->pos, r);
 
            // Error vector: from tip to neighbor center
            b3Vec3 e = b3Sub(aj->pos, tip);
            float C = b3Length(e);
            if (C < 1e-8f) continue;
 
            b3Vec3 n = b3MulSV(1.0f / C, e);
 
            // Relative velocity along constraint normal
            // u = n . (v_j - (v_i + omega_i x r))
            b3Vec3 tipVel = b3Add(ai->vel, b3Cross(ai->omega, r));
            float u = b3Dot(n, b3Sub(aj->vel, tipVel));
 
            // Effective inverse mass
            // w_port = invM_i + invM_j + (r x n)^T * invI_i * (r x n)
            b3Vec3 s = b3Cross(r, n);
            float wport = ai->invMass + aj->invMass + b3Dot(s, b3MulMV(ai->invInertia, s));
 
            // Local implicit impulse (backward-Euler spring)
            // a = h^2 * K * w,  d = h * D * w
            // p = -(a*C/h + (a+d)*u) / (w * (1 + a + d))
            // Equivalently using algorithmic parameters directly:
            float a = w->aBond;   // dimensionless correction fraction
            float d = w->dBond;   // dimensionless damping
            float p = -(a * C / h + (a + d) * u) / (wport * (1.0f + a + d));
 
            b3Vec3 P = b3MulSV(p, n);
 
            // Accumulate on central atom
            Jacc = b3Add(Jacc, P);
            Lacc = b3Add(Lacc, b3Cross(r, P));
 
            // Store recoil for neighbor
            w->Jrecoil[i * MAX_PORTS + k] = b3MulSV(-1.0f, P);
        }
 
        w->Jown[i] = Jacc;
        w->Lown[i] = Lacc;
    }
}
```
 
**Momentum conservation**:
- Linear: `J_i + J_j = P - P = 0` (exact, since recoil is `-P`)
- Angular: `(t_ik - x_j) x P = -e x (p*n) = 0` because `P || e` (scalar radial reaction)
- No torque on neighbor `j` (directed port acts on neighbor's center)
 
### 3.4 Collision Impulse (Soft Unilateral)
 
```c
void gatherCollisionImpulses(RigidAtomWorld* w, float h) {
    // Option A: iterate over collision pairs once, write to both atoms
    for (int p = 0; p < w->nCollisionPairs; p++) {
        int i = w->collisionPairs[p].i;
        int j = w->collisionPairs[p].j;
        float Rij = w->collisionPairs[p].Rij;
 
        RigidAtom* ai = &w->atoms[i];
        RigidAtom* aj = &w->atoms[j];
 
        b3Vec3 d = b3Sub(ai->pos, aj->pos);
        float r = b3Length(d);
        if (r < 1e-8f) continue;
 
        float delta = Rij - r;  // penetration (>0 means overlapping)
        if (delta <= 0.0f) continue;
 
        b3Vec3 n = b3MulSV(1.0f / r, d);  // from j to i
 
        // Separating velocity (positive = separating)
        float uc = b3Dot(n, b3Sub(ai->vel, aj->vel));
 
        // Effective inverse mass (sphere-sphere: no rotation contribution)
        float wc = ai->invMass + aj->invMass;
 
        // Local implicit impulse
        float a = w->aCollision;
        float d_param = w->dCollision;
        // For collision: C = -delta (negative because penetration is overlap)
        // We want repulsive impulse along +n (pushes i away from j)
        // p = -(a*(-delta)/h + (a+d)*(-uc)) / (wc * (1+a+d))
        //   = (a*delta/h - (a+d)*uc) / (wc * (1+a+d))
        float p = (a * delta / h - (a + d_param) * uc) / (wc * (1.0f + a + d_param));
 
        // Clamp to non-negative (no attractive contact)
        if (p < 0.0f) p = 0.0f;
 
        b3Vec3 P = b3MulSV(p, n);
 
        // Apply to both atoms
        w->Jcoll[i] = b3Add(w->Jcoll[i], P);
        w->Jcoll[j] = b3Add(w->Jcoll[j], b3MulSV(-1.0f, P));
    }
}
```
 
**Momentum conservation**:
- Linear: `P_c - P_c = 0` (exact)
- Angular: `(x_i - x_j) x P_c = d x (p*n) = 0` because `P || d` (no torque for sphere-sphere)
 
### 3.5 Atom Update (Gather-Apply)
 
```c
void applyAtomUpdate(RigidAtomWorld* w, float h) {
    for (int i = 0; i < w->nAtoms; i++) {
        RigidAtom* ai = &w->atoms[i];
 
        // Gather all linear impulses
        b3Vec3 J = w->Jown[i];  // own outgoing port reactions
 
        // Add incoming port recoils
        for (int s = 0; s < w->incomingCount[i]; s++) {
            int slot = w->incomingOffset[i] + s;
            int srcAtom = w->incomingAtoms[slot];
            int srcPort = w->incomingPorts[slot];
            J = b3Add(J, w->Jrecoil[srcAtom * MAX_PORTS + srcPort]);
        }
 
        // Add collision impulses
        J = b3Add(J, w->Jcoll[i]);
 
        // Angular impulse (only from own ports; collisions have no torque)
        b3Vec3 L = w->Lown[i];
 
        // Update velocities
        ai->vel   = b3Add(ai->vel,   b3MulSV(ai->invMass, J));
        ai->omega = b3Add(ai->omega, b3MulMV(ai->invInertia, L));
 
        // Integrate position and rotation
        ai->pos  = b3Add(ai->pos, b3MulSV(h, ai->vel));
        ai->quat = b3NormalizeQuat(b3IntegrateRotation(ai->quat, b3MulSV(h, ai->omega)));
    }
}
```

---

## 4. Soft Force Evaluation

The `evaluateSoftForces` function computes O(N^2) pairwise interactions:
- Coulomb electrostatics
- Lennard-Jones / Morse attractive part
- Soft torsion / dihedral terms

This is the expensive evaluation done **once per macrostep**. The force is cached and reused as a frozen kick across all `nsub` microsteps.

For minimization (not dynamics), replace the velocity update with a FIRE optimizer or heavy-ball momentum scheme. The hard SI projector remains the same.

---

## 5. Parameter Guide

All `a` and `d` parameters are **dimensionless algorithmic parameters**, not physical force constants.

| Parameter | Meaning | Starting Value | Tuning |
|-----------|---------|----------------|--------|
| `aBond` | Correction fraction per bond visit | 10 | Higher = stiffer bonds. `a=1` removes ~50% error, `a=9` removes ~90%, `a->inf` = hard constraint |
| `dBond` | Bond velocity damping | 1-3 | Prevents oscillation. `d=0` = undamped |
| `aAngle` | Correction fraction for ports/angles | 3-10 | Lower than bonds since ports are more coupled |
| `dAngle` | Port/angle damping | 1-3 | |
| `aCollision` | Correction fraction for contacts | 5-20 | Higher = harder spheres |
| `dCollision` | Contact damping | 2-5 | Prevents bouncing |
| `H` | Soft macrostep | Tune | Limited by soft force staleness and `max|H*v| < dx_max` |
| `nsub` | Hard microsteps per macrostep | 4-16 | Increase when velocities are high: `h*|v_ij| < 0.2*R_ij` |

### Dimensionless interpretation

For an isolated undamped constraint:
```
C' = (C + h*u) / (1 + a)
```
- `a = 1`: removes ~50% of predicted error
- `a = 3`: removes ~75%
- `a = 9`: removes ~90%
- `a -> inf`: hard projection

You can specify a desired correction fraction `rho = 1/(1+a)` and derive `a = (1-rho)/rho`.

---

## 6. Implementation Phases

### Phase 1: Core SI Solver (standalone C)

Create `src/rigid_atom_ff.h` and `src/rigid_atom_ff.c`:

1. **Data structures** — `RigidAtom`, `RigidPort`, `RigidAtomWorld` as above
2. **Port impulse gather** — `gatherPortImpulses()` with local implicit formula
3. **Collision impulse gather** — `gatherCollisionImpulses()` with unilateral clamp
4. **Atom update** — `applyAtomUpdate()` aggregating all impulses
5. **Hard microstep** — `hardMicrostep()` orchestrating gather + apply
6. **Macrostep** — `RigidAtomWorld_Step()` with RESPA splitting
7. **Soft force stub** — placeholder O(N^2) Lennard-Jones or Coulomb
8. **Broad phase** — simple O(N^2) pair search initially; upgrade to spatial hash later
9. **Reverse adjacency** — precompute `incomingCount`, `incomingAtoms`, `incomingPorts` from port topology

### Phase 2: Test Harness

Create `test/test_rigid_atom_ff.c`:

1. **Two-atom bond test** — verify bond length converges, momentum conserved
2. **Three-atom angle test** — verify angle converges, momentum conserved
3. **Collision test** — two atoms approaching, verify repulsion and momentum
4. **Chain relaxation** — small polymer chain, verify stability with large `H`
5. **Momentum conservation check** — assert `sum(p) = 0` and `sum(L) = 0` after N steps

### Phase 3: Box3D Visualization Integration

Create `samples/sample_rigid_atom_ff.cpp`:

1. **Inherit from `Sample`** — follow pattern from `sample_molecules.cpp`
2. **Load molecule** — reuse `b3LoadMolFile` from mol_loader
3. **Build `RigidAtomWorld`** — convert mol file atoms/bonds to `RigidAtom` + `RigidPort` arrays
4. **Custom `Step()`** — call `RigidAtomWorld_Step()` instead of `b3World_Step()`
5. **Draw atoms** — use `DrawSphere` from `gfx/draw.h` for each atom at `pos[i]` with `radius[i]`
6. **Draw bonds** — use `DrawLine` or `DrawCylinder` between bonded atom centers
7. **UI controls** — ImGui sliders for `H`, `nsub`, `aBond`, `aCollision`, damping
8. **Debug viz** — optional: draw port tips, collision normals, velocity arrows

Key integration points:
- `Sample` base class in `samples/sample.h`
- `DrawSphere`, `DrawLine` from `samples/gfx/draw.h`
- `b3LoadMolFile` from mol_loader (already used in `sample_molecules.cpp`)
- Camera setup via `m_camera->SetView()`

### Phase 4: Optimization

1. **Spatial hash** for broad-phase collisions — replace O(N^2) pair search
2. **Graph coloring** for parallel port sweeps — color atoms so no two bonded atoms share a color
3. **Collision pair update frequency** — rebuild every K microsteps instead of every one
4. **Optional: 6x6 VBD block solve** — upgrade path for coupled atoms, solving translation+rotation together as one block instead of separate scalar impulses
5. **Optional: FIRE optimizer** — replace simple velocity damping with FIRE (Fast Inertial Relaxation Engine) for minimization tasks

---

## 7. Key Files

- `include/box3d/math_functions.h` — `b3Vec3`, `b3Quat`, `b3Matrix3`, `b3RotateVector`, `b3Solve3`, `b3MulQuat`, `b3NormalizeQuat`, `b3MakeQuatFromAxisAngle`
- `src/math_internal.h` — `b3IntegrateRotation`, `b3QuatFromExponentialMap`, `b3DeltaQuatToRotation`
- `samples/sample.h` — `Sample` base class, `SampleContext`
- `samples/sample_molecules.cpp` — reference for mol loading + visualization pattern
- `samples/gfx/draw.h` — `DrawSphere`, `DrawLine`, `DrawCylinder` rendering API
- `docs/prokop/RigidAtomFF_SI.chat.md` — full discussion (4067 lines)
- `docs/prokop/SI_for_molecules.chat.md` — SI comparison and local implicit impulse derivation
- `docs/prokop/overview.md` — Box3D engine overview
- `docs/prokop/visualization.md` — Box3D visualization system


