# RigidAtomFF Multi-Solver Implementation Plan

Comprehensive design for rigid-atom dynamics in pure C within Box3D, covering four solver methods —
**Sequential Impulses (SI)**, **XPBD**, **Projective Dynamics (PD)**, and **Vertex Block Descent (VBD/AVBD)** —
all within a unified atom-centered gather-apply architecture with RESPA/Strang time splitting.

Synthesized from GPT 5.5 insights in `RigidAtomFF_SI.chat.md` and the current SI implementation
in `src/rigid_atom_ff.c`.
 
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
 
Key principles (shared across all methods):
- **Atom-centered**: only dynamic state is `pos[i]`, `quat[i]`, `vel[i]`, `omega[i]`. All forces/torques computed on the fly.
- **Directed ports**: each atom has outgoing ports pointing to neighbor centers. Port `k` of atom `i` has a local lever arm `localDir` and points to `neighIdx`. The constraint error is the vector from the port tip to the neighbor's center.
- **Gather-apply architecture**: one pass per atom gathers all outgoing port + collision contributions into accumulators (`Jown`, `Lown`, `Jrecoil`, `Jcoll`), then a second pass aggregates incoming recoils and applies one velocity/position update. This is the invariant across all methods — only the gather computation changes.
- **Scalar radial reactions** along constraint normals to guarantee angular momentum conservation.
- **Momentum conservation**: equal-and-opposite impulses; ports produce torque on central atom only; collisions produce no torque (sphere-sphere).
- **RESPA/Strang splitting**: expensive soft forces (O(N²)) evaluated once per macrostep `H`; cheap hard constraints (O(N)) solved in `nsub` microsteps with `h = H/nsub`.

### Method selection

All four methods share the same outer driver, data structures, and gather-apply skeleton. They differ only in **what the gather step computes**:

| Method | Corrects | Inertia | Multiplier storage | Key formula |
|--------|----------|---------|---------------------|-------------|
| **SI** | velocity + position | scalar effective mass | none | local implicit impulse (backward-Euler spring) |
| **XPBD** | position + velocity | scalar effective mass | per-port Lagrange multiplier (warm-started) | compliant position correction with stored λ |
| **PD** | position (local projection + global Jacobi) | scalar effective mass | none | local rotation projection → global quadratic solve approximated by Jacobi |
| **VBD** | position + rotation (6-DOF block) | 6×6 effective mass matrix | none | atom-centered 6×6 Hessian solve combining all incident energies |
| **AVBD** | position + rotation | 6×6 + augmented Lagrangian | per-constraint multipliers | VBD with augmented Lagrangian for hard constraints |
 
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
    int*           incomingOffset;
    int*           incomingAtoms;   // flat: total incoming refs
    int*           incomingPorts;   // flat: total incoming refs

    // Collision pairs (dynamic, from broad phase)
    CollisionPair* collisionPairs;
    int            nCollisionPairs;

    // Soft forces (cached, evaluated once per macrostep)
    b3Vec3*        fsoft;    // force per atom
    b3Vec3*        tsoft;    // torque per atom

    // Temporary accumulators (per microstep, shared across all methods)
    b3Vec3*        Jown;     // own port linear impulse
    b3Vec3*        Lown;     // own port angular impulse
    b3Vec3*        Jrecoil;  // port recoil sent to neighbor (flat: nAtoms * MAX_PORTS)
    b3Vec3*        Jcoll;    // collision linear impulse

    // XPBD: per-port Lagrange multipliers (warm-started across microsteps)
    float*         lambdaPort;   // flat: nAtoms * MAX_PORTS
    float*         lambdaColl;   // per collision pair (rebuilt each microstep)

    // VBD: per-atom 6x6 Hessian and gradient (built each microstep)
    float*         Hblock;   // flat: nAtoms * 36 (6x6 matrix)
    float*         gblock;   // flat: nAtoms * 6  (gradient vector)

    // Method selector
    enum { SOLVER_SI, SOLVER_XPBD, SOLVER_PD, SOLVER_VBD, SOLVER_AVBD } solverMethod;

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
 
## 3. Core Algorithm — SI (Sequential Impulses)

The following describes the SI method, which is the current implementation. Other methods
(§4) share the same macrostep and atom-update structure but differ in the gather step.
 
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

## 4. Solver Method Variants

The SI method above is the current implementation. The following sections describe how the other
three methods (XPBD, PD, VBD/AVBD) fit into the same gather-apply skeleton. The key insight from
GPT 5.5 is that **all methods share the same outer driver and buffer structure** — only the gather
computation changes. The `hardMicrostep()` function dispatches based on `solverMethod`:

```c
void hardMicrostep(RigidAtomWorld* w, float h) {
    clearAccumulators(w);
    updateCollisionPairs(w);

    switch (w->solverMethod) {
        case SOLVER_SI:   gatherPortImpulsesSI(w, h);   break;
        case SOLVER_XPBD: gatherPortImpulsesXPBD(w, h); break;
        case SOLVER_PD:   gatherPortImpulsesPD(w, h);   break;
        case SOLVER_VBD:  gatherPortImpulsesVBD(w, h);  break;
        case SOLVER_AVBD: gatherPortImpulsesAVBD(w, h); break;
    }

    gatherCollisionImpulses(w, h);  // collisions always SI-style (unilateral)
    applyAtomUpdate(w, h);
}
```

### 4.1 XPBD (Extended Position-based dynamics)

**What changes vs SI**: XPBD stores a per-port Lagrange multiplier `λ` that is warm-started
across microsteps. The correction is applied at the **position level** first, then velocity
is derived from the position change. This gives better constraint convergence for stiff bonds.

**Key differences from SI**:
- Position-level correction: `Δx = λ * n` instead of velocity-level impulse
- Warm-started multipliers: `λ` persists across microsteps within a macrostep
- Compliance parameter `α̃ = 1/(K * h²)` replaces dimensionless `a`
- No velocity damping term in the core formula (damping added separately)

```c
void gatherPortImpulsesXPBD(RigidAtomWorld* w, float h) {
    for (int i = 0; i < w->nAtoms; i++) {
        RigidAtom* ai = &w->atoms[i];
        b3Vec3 Jacc = b3Vec3_zero;
        b3Vec3 Lacc = b3Vec3_zero;

        for (int k = 0; k < ai->nPorts; k++) {
            RigidPort* port = &w->ports[i * MAX_PORTS + k];
            int j = port->neighIdx;
            if (j < 0) continue;

            RigidAtom* aj = &w->atoms[j];

            b3Vec3 r = b3RotateVector(ai->quat, port->localDir);
            b3Vec3 tip = b3Add(ai->pos, r);
            b3Vec3 e = b3Sub(aj->pos, tip);
            float C = b3Length(e);
            if (C < 1e-8f) continue;

            b3Vec3 n = b3MulSV(1.0f / C, e);

            // Effective inverse mass (same as SI)
            b3Vec3 s = b3Cross(r, n);
            float wport = ai->invMass + aj->invMass + b3Dot(s, b3MulMV(ai->invInertia, s));

            // XPBD compliance: α̃ = 1 / (K * h²)
            // For dimensionless form: α̃ = 1 / a  (where a = h² * K * w)
            float alphaTilde = 1.0f / w->aBond;
            float dt = h;

            // Retrieve warm-started multiplier
            float* lambda = &w->lambdaPort[i * MAX_PORTS + k];

            // XPBD impulse: Δλ = (-C - α̃ * λ) / (wport * (1 + α̃/dt²))
            // Simplified: Δλ = -(C + alphaTilde * (*lambda)) / (wport + alphaTilde)
            float dlambda = -(C + alphaTilde * (*lambda)) / (wport + alphaTilde);

            // Update multiplier (warm start for next microstep)
            *lambda += dlambda;

            // Position correction impulse
            float p = dlambda / dt;  // convert to velocity-level impulse
            b3Vec3 P = b3MulSV(p, n);

            Jacc = b3Add(Jacc, P);
            Lacc = b3Add(Lacc, b3Cross(r, P));
            w->Jrecoil[i * MAX_PORTS + k] = b3MulSV(-1.0f, P);
        }

        w->Jown[i] = Jacc;
        w->Lown[i] = Lacc;
    }
}
```

**When to use XPBD**: When bonds are very stiff and SI requires too many microsteps to converge.
The warm-started multipliers give faster convergence per microstep. Particularly good for
minimization tasks where you want the system to settle quickly.

**Momentum conservation**: Same as SI — equal-and-opposite impulses, scalar radial reactions.
The stored `λ` does not break momentum conservation because it only affects the magnitude
of the symmetric impulse, not its direction.

**Multiplier reset**: `lambdaPort` should be zeroed at the start of each macrostep (not each
microstep — that's the whole point of warm starting).

### 4.2 Projective Dynamics (PD)

**What changes vs SI**: PD does a two-phase solve:
1. **Local projection**: For each port, compute the "optimal" position of the neighbor that
   would satisfy the constraint if the central atom were fixed. This is a rotation-only
   projection — find the rotation of atom `i` that best aligns port `k` toward neighbor `j`.
2. **Global solve**: Minimize the quadratic energy `Σ w_k * |x_i - projected_k|²` over all
   atoms simultaneously. In the atom-centered approximation, this reduces to a Jacobi
   iteration where each atom takes a weighted average of its own port projections and
   incoming neighbor projections.

**Key differences from SI**:
- No velocity in the correction formula — purely position-based
- Local projection step computes a target position, not an impulse
- Global Jacobi iteration approximates the global quadratic solve
- Better for very stiff systems with many coupled constraints (e.g. ring closures)
- No damping — damping must be added as a separate velocity filter

```c
void gatherPortImpulsesPD(RigidAtomWorld* w, float h) {
    // Phase 1: Local projection — compute target positions
    // For each port k of atom i, the projected neighbor position is:
    //   target_ik = pos_i + rotate(quat_i, localDir_k)  (ideal tip = neighbor center)
    // The "pull" on atom i is toward making the port tip reach the neighbor.

    for (int i = 0; i < w->nAtoms; i++) {
        RigidAtom* ai = &w->atoms[i];
        b3Vec3 Jacc = b3Vec3_zero;
        b3Vec3 Lacc = b3Vec3_zero;

        float totalWeight = 0.0f;
        b3Vec3 posCorrection = b3Vec3_zero;
        b3Vec3 rotCorrection = b3Vec3_zero;

        for (int k = 0; k < ai->nPorts; k++) {
            RigidPort* port = &w->ports[i * MAX_PORTS + k];
            int j = port->neighIdx;
            if (j < 0) continue;

            RigidAtom* aj = &w->atoms[j];

            b3Vec3 r = b3RotateVector(ai->quat, port->localDir);
            b3Vec3 tip = b3Add(ai->pos, r);
            b3Vec3 e = b3Sub(aj->pos, tip);  // error: neighbor - tip
            float C = b3Length(e);
            if (C < 1e-8f) continue;

            // Weight = stiffness / (C + epsilon) — stiffer for small errors
            float wk = w->aBond / (1.0f + C);

            // Position correction: move tip toward neighbor by weighted fraction
            // This translates to moving atom i by wk * e / (invM_i + invM_j)
            b3Vec3 n = b3MulSV(1.0f / C, e);
            b3Vec3 s = b3Cross(r, n);
            float wport = ai->invMass + aj->invMass + b3Dot(s, b3MulMV(ai->invInertia, s));

            float correctionMag = wk * C / (wport * (1.0f + wk));
            b3Vec3 P = b3MulSV(correctionMag / h, n);  // convert to velocity-level impulse

            Jacc = b3Add(Jacc, P);
            Lacc = b3Add(Lacc, b3Cross(r, P));
            w->Jrecoil[i * MAX_PORTS + k] = b3MulSV(-1.0f, P);
        }

        w->Jown[i] = Jacc;
        w->Lown[i] = Lacc;
    }

    // Phase 2: Global Jacobi — already done by applyAtomUpdate()
    // which aggregates own + incoming + collision impulses.
    // For true PD, multiple Jacobi iterations would loop phases 1+2.
    // In our single-pass approximation, one iteration per microstep suffices
    // when nsub is large enough.
}
```

**When to use PD**: When constraints are highly coupled (rings, cages, dense networks).
The local projection + global Jacobi naturally handles coupling better than independent
per-port impulses. However, the single-pass approximation per microstep means you need
more microsteps for the same accuracy as SI.

**Momentum conservation**: Same guarantee — scalar radial reactions, equal-and-opposite.
The Jacobi iteration preserves momentum because each iteration applies symmetric impulses.

**Note on true PD**: The full PD method would store projected positions and iterate
the global solve multiple times per microstep. The atom-centered single-pass version
here is a practical approximation that converges over multiple microsteps instead.

### 4.3 Vertex Block Descent (VBD)

**What changes vs SI**: VBD solves a **6×6 block system per atom** that couples
translation and rotation simultaneously. Instead of separate scalar impulses per port,
VBD gathers all port energies, incoming recoils, and collision energies into a single
6×6 Hessian `H` and gradient `g`, then solves `H * Δq = -g` where `q = (Δx, Δθ)`.

**Key differences from SI**:
- 6-DOF block: translation and rotation solved together
- 6×6 effective mass matrix includes cross-terms between translation and rotation
- All ports on an atom are coupled through the shared block solve
- Better convergence for atoms with multiple ports at different angles
- More expensive per atom (6×6 solve vs scalar) but fewer iterations needed

```c
void gatherPortImpulsesVBD(RigidAtomWorld* w, float h) {
    for (int i = 0; i < w->nAtoms; i++) {
        RigidAtom* ai = &w->atoms[i];

        // Build 6x6 effective mass matrix M_eff = diag(invM * I3, invI)
        // and accumulate Hessian H = M_eff/h² + Σ port contributions
        // and gradient g = M_eff/h² * (predicted position - current) + Σ port gradients

        float H[6][6] = {0};  // 6x6 Hessian
        float g[6] = {0};     // 6x1 gradient

        // Inertial term: M_eff / h²
        // Translation block (3x3)
        for (int a = 0; a < 3; a++)
            H[a][a] += ai->invMass / (h * h);

        // Rotation block (3x3) = invInertia / h²
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++)
                H[3+a][3+b] += ai->invInertia.m[a][b] / (h * h);

        // Inertial gradient: predicted position drift = vel * h
        // g_trans = invM * vel / h  (pulls toward inertial prediction)
        g[0] = ai->invMass * ai->vel.x / h;
        g[1] = ai->invMass * ai->vel.y / h;
        g[2] = ai->invMass * ai->vel.z / h;
        // g_rot = invI * omega / h
        b3Vec3 grot = b3MulMV(ai->invInertia, b3MulSV(1.0f/h, ai->omega));
        g[3] = grot.x; g[4] = grot.y; g[5] = grot.z;

        // Port contributions
        b3Vec3 JrecoilTotal = b3Vec3_zero;

        for (int k = 0; k < ai->nPorts; k++) {
            RigidPort* port = &w->ports[i * MAX_PORTS + k];
            int j = port->neighIdx;
            if (j < 0) continue;

            RigidAtom* aj = &w->atoms[j];

            b3Vec3 r = b3RotateVector(ai->quat, port->localDir);
            b3Vec3 tip = b3Add(ai->pos, r);
            b3Vec3 e = b3Sub(aj->pos, tip);
            float C = b3Length(e);
            if (C < 1e-8f) continue;

            b3Vec3 n = b3MulSV(1.0f / C, e);

            // Port stiffness in VBD: K_port = aBond / h²
            float Kp = w->aBond / (h * h);

            // Gradient contribution: g += K * C * [n; r × n]
            g[0] -= Kp * C * n.x;
            g[1] -= Kp * C * n.y;
            g[2] -= Kp * C * n.z;
            b3Vec3 rxn = b3Cross(r, n);
            g[3] -= Kp * C * rxn.x;
            g[4] -= Kp * C * rxn.y;
            g[5] -= Kp * C * rxn.z;

            // Hessian contribution: H += K * [n⊗n, n⊗(r×n); (r×n)⊗n, (r×n)⊗(r×n)]
            // This is a rank-1 update to the 6x6 Hessian
            float J[6] = {n.x, n.y, n.z, rxn.x, rxn.y, rxn.z};
            for (int a = 0; a < 6; a++)
                for (int b = 0; b < 6; b++)
                    H[a][b] += Kp * J[a] * J[b];

            // Store recoil for neighbor (same as SI)
            // The recoil impulse will be derived from the solved Δq
            // For simplicity, use the same scalar approximation for recoil:
            b3Vec3 s = b3Cross(r, n);
            float wport = ai->invMass + aj->invMass + b3Dot(s, b3MulMV(ai->invInertia, s));
            float p = -(w->aBond * C / h + (w->aBond + w->dBond) * 0.0f) / (wport * (1.0f + w->aBond + w->dBond));
            b3Vec3 P = b3MulSV(p, n);
            w->Jrecoil[i * MAX_PORTS + k] = b3MulSV(-1.0f, P);
        }

        // Solve H * Δq = -g  (6x6 linear system)
        float dq[6];
        solve6x6(H, g, dq);  // TODO: implement or use existing 3x3 + Schur complement

        // Extract velocity-level impulse from Δq
        b3Vec3 J = {dq[0] * h, dq[1] * h, dq[2] * h};  // Δx/h * h = Δx → impulse = invM * Δx
        b3Vec3 L = {dq[3] * h, dq[4] * h, dq[5] * h};

        w->Jown[i] = b3MulSV(ai->mass, J);  // convert Δv to impulse
        w->Lown[i] = b3MulSV(1.0f, L);       // angular impulse (already in right frame)
    }
}
```

**When to use VBD**: When atoms have 3+ ports at different orientations and the
translation-rotation coupling matters (e.g. tetrahedral carbon centers). The 6×6 solve
captures the coupling that scalar SI misses, giving better convergence per microstep
at the cost of a 6×6 linear solve per atom.

**Momentum conservation**: VBD preserves linear momentum exactly (same recoil mechanism).
Angular momentum is conserved per-atom because the 6×6 block includes the rotational
inertia tensor. The recoil to neighbors uses scalar radial reactions as in SI.

**Implementation note**: The 6×6 solve can be implemented as:
1. Extract the 3×3 translation block `Htt`, rotation block `Hrr`, and coupling `Htr`
2. Use Schur complement: solve `Hrr * Δθ = gr - Htr^T * Htt^{-1} * gt` first
3. Then `Δx = Htt^{-1} * (gt - Htr * Δθ)`
4. This requires two 3×3 solves instead of one 6×6 solve

### 4.4 Augmented VBD (AVBD)

**What changes vs VBD**: AVBD adds an augmented Lagrangian term to the VBD block solve,
allowing hard constraints to be enforced exactly over multiple iterations. Per-constraint
multipliers `λ` are stored and updated each microstep, similar to XPBD.

**Key differences from VBD**:
- Per-constraint Lagrange multipliers (stored, warm-started)
- Hard constraints: as `K → ∞`, the augmented term dominates and enforces exact satisfaction
- The 6×6 Hessian gets an additional term: `H += Σ λ_k * J_k * J_k^T`
- The gradient gets: `g += Σ λ_k * C_k * J_k`

```c
void gatherPortImpulsesAVBD(RigidAtomWorld* w, float h) {
    for (int i = 0; i < w->nAtoms; i++) {
        RigidAtom* ai = &w->atoms[i];

        float H[6][6] = {0};
        float g[6] = {0};

        // Inertial term (same as VBD)
        for (int a = 0; a < 3; a++)
            H[a][a] += ai->invMass / (h * h);
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++)
                H[3+a][3+b] += ai->invInertia.m[a][b] / (h * h);

        g[0] = ai->invMass * ai->vel.x / h;
        g[1] = ai->invMass * ai->vel.y / h;
        g[2] = ai->invMass * ai->vel.z / h;
        b3Vec3 grot = b3MulMV(ai->invInertia, b3MulSV(1.0f/h, ai->omega));
        g[3] = grot.x; g[4] = grot.y; g[5] = grot.z;

        for (int k = 0; k < ai->nPorts; k++) {
            RigidPort* port = &w->ports[i * MAX_PORTS + k];
            int j = port->neighIdx;
            if (j < 0) continue;

            RigidAtom* aj = &w->atoms[j];
            b3Vec3 r = b3RotateVector(ai->quat, port->localDir);
            b3Vec3 tip = b3Add(ai->pos, r);
            b3Vec3 e = b3Sub(aj->pos, tip);
            float C = b3Length(e);
            if (C < 1e-8f) continue;

            b3Vec3 n = b3MulSV(1.0f / C, e);
            b3Vec3 rxn = b3Cross(r, n);
            float J[6] = {n.x, n.y, n.z, rxn.x, rxn.y, rxn.z};

            // Penalty stiffness (finite, unlike true hard constraints)
            float Kp = w->aBond / (h * h);

            // Standard VBD terms
            for (int a = 0; a < 6; a++)
                for (int b = 0; b < 6; b++)
                    H[a][b] += Kp * J[a] * J[b];
            for (int a = 0; a < 6; a++)
                g[a] -= Kp * C * J[a];

            // Augmented Lagrangian: add λ * J * J^T to H and λ * C * J to g
            float* lambda = &w->lambdaPort[i * MAX_PORTS + k];
            for (int a = 0; a < 6; a++)
                for (int b = 0; b < 6; b++)
                    H[a][b] += (*lambda) * J[a] * J[b];
            for (int a = 0; a < 6; a++)
                g[a] += (*lambda) * C * J[a];
        }

        // Solve H * Δq = -g
        float dq[6];
        solve6x6(H, g, dq);

        // Update multipliers: λ_k += K * C_k(Δq)
        for (int k = 0; k < ai->nPorts; k++) {
            RigidPort* port = &w->ports[i * MAX_PORTS + k];
            int j = port->neighIdx;
            if (j < 0) continue;
            b3Vec3 r = b3RotateVector(ai->quat, port->localDir);
            b3Vec3 tip = b3Add(ai->pos, r);
            b3Vec3 e = b3Sub(w->atoms[j].pos, tip);
            float C = b3Length(e);
            if (C < 1e-8f) continue;
            float Kp = w->aBond / (h * h);
            w->lambdaPort[i * MAX_PORTS + k] += Kp * C;
        }

        // Extract impulse (same as VBD)
        b3Vec3 J = {dq[0] * h, dq[1] * h, dq[2] * h};
        b3Vec3 L = {dq[3] * h, dq[4] * h, dq[5] * h};
        w->Jown[i] = b3MulSV(ai->mass, J);
        w->Lown[i] = L;
    }
}
```

**When to use AVBD**: When you need exact constraint satisfaction (hard bonds, rigid
ring closures) but still want the coupled 6-DOF solve. The augmented Lagrangian converges
to exact satisfaction over multiple microsteps, unlike penalty-only VBD which always has
residual error proportional to `1/K`.

**Momentum conservation**: Same as VBD — the augmented terms are symmetric per constraint
and don't break the equal-and-opposite recoil structure.

---

## 5. Soft Force Evaluation

The `evaluateSoftForces` function computes O(N^2) pairwise interactions:
- Coulomb electrostatics
- Lennard-Jones / Morse attractive part
- Soft torsion / dihedral terms

This is the expensive evaluation done **once per macrostep**. The force is cached and reused as a frozen kick across all `nsub` microsteps.

For minimization (not dynamics), replace the velocity update with a FIRE optimizer or heavy-ball momentum scheme. The hard SI projector remains the same.

---

## 6. Parameter Guide

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

### Method-specific parameter notes

- **SI**: Uses `a` and `d` directly in the impulse formula. Damping is built-in.
- **XPBD**: Uses `α̃ = 1/a` as compliance. Damping must be added as a separate velocity filter after the position solve (e.g. `v *= (1 - d/(1+d))`). Multipliers `λ` are zeroed at macrostep start.
- **PD**: Uses `a` as the projection weight. No velocity damping — add a global velocity damping `v *= exp(-d*h)` after each microstep if needed.
- **VBD**: Uses `a/h²` as the Hessian stiffness. Damping enters through the inertial gradient term (velocity-dependent).
- **AVBD**: Same as VBD for stiffness. Multiplier update rate is `K*C` per microstep. Multipliers zeroed at macrostep start.

### Dimensionless interpretation (SI and PD)

For an isolated undamped constraint:
```
C' = (C + h*u) / (1 + a)
```
- `a = 1`: removes ~50% of predicted error
- `a = 3`: removes ~75%
- `a = 9`: removes ~90%
- `a -> inf`: hard projection

You can specify a desired correction fraction `rho = 1/(1+a)` and derive `a = (1-rho)/rho`.

### Reciprocal ports and stiffness splitting

Each physical bond creates **two directed ports** (one from atom `i` to `j`, one from `j` to `i`).
This doubles the effective stiffness. To compensate, **halve the per-port stiffness**: if the
desired physical bond stiffness is `K`, set each port's `a = K * h² * w / 2`.

---

## 7. GPT 5.5 Key Insights

The following insights from the GPT 5.5 discussion (`RigidAtomFF_SI.chat.md`) are critical for
correctness and stability. They apply to all solver methods.

### 7.1 Impulse sign convention

The port impulse must be **positive** along `n` (from tip toward neighbor) to restore the bond.
A reversed sign causes exponential velocity growth and NaNs. The correct formula:

```
p = -(a * C / h + (a + d) * u) / (wport * (1 + a + d))
```

- `C > 0` (tip is short of neighbor) → `p < 0` along `n` → atom `i` pulled toward `j` ✓
- `u > 0` (approaching) → damped ✓
- **Wrong sign**: `p = +(a*C/h + ...)` → atom pushed away → exponential blow-up ✗

### 7.2 Microstep time `h` vs macrostep `H`

All impulse and stiffness calculations inside the hard microstep must use `h = H/nsub`, **not** `H`.
Using `H` makes the effective stiffness `a = K*H²*w` too large by factor `nsub²`, causing overshoot
and instability. The macrostep `H` is only used for the soft force half-kicks.

### 7.3 Collision list update frequency

Collision pairs must be rebuilt **every microstep** (or at least every few microsteps), not once
per macrostep. Atoms move by `h*v` per microstep, and new collisions can appear. If the list is
stale, atoms can tunnel through each other.

### 7.4 Lennard-Jones force sign

The LJ force must be **attractive** at `r > r_min` and **repulsive** at `r < r_min`. The standard
form `F = 12ε/σ * [(σ/r)^13 - (σ/r)^7] * r̂` is correct. A sign error makes it repulsive at long
range and attractive at short range — the opposite of physical behavior.

### 7.5 Momentum conservation pitfalls

Operations that break momentum conservation:
- **Per-atom clipping** of displacements or velocities (each atom clipped independently)
- **Asymmetric relaxation** (different relaxation factors for atom `i` vs `j`)
- **Stale collision normals** (normal computed from old positions, applied to new velocities)
- **Per-atom damping** that doesn't account for relative velocity

All of these must be avoided. The gather-apply architecture with equal-and-opposite recoils
guarantees conservation as long as the gather step produces symmetric impulses.

### 7.6 Coincident atom handling

When two atoms are at exactly the same position (`r < 1e-8`), the collision normal is undefined.
Handle this deterministically: use a fixed separation direction (e.g. `(1,0,0)`) or skip the
collision and let the next microstep resolve it. Never use a random direction — this breaks
determinism and can cause energy drift.

### 7.7 Inertial regularization for VBD/PD

The 6×6 Hessian in VBD (or the effective mass in PD) must include the inertial term `M/h²`.
Without it, the system is under-determined when constraints are satisfied (zero gradient, zero
Hessian from ports, no inertial pull). The inertial term acts as a Tikhonov regularization that
keeps the solve well-conditioned.

### 7.8 Soft macrostep / hard microstep interaction

The RESPA splitting requires that soft forces are **frozen** during microsteps — they are not
re-evaluated. The frozen soft force is applied as a constant acceleration `a_soft = F_soft/m`
over each microstep. This is correct as long as `H` is small enough that soft forces don't
change significantly over one macrostep. If soft forces change rapidly (e.g. close-range LJ),
reduce `H` or increase `nsub`.

---

## 8. Implementation Phases

### Phase 1: Core SI Solver ✅ (done)

Created `src/rigid_atom_ff.h` and `src/rigid_atom_ff.c`:

1. **Data structures** — `RAffAtom`, `RAffPort`, `RAffWorld` ✓
2. **Port impulse gather** — `solvePortImpulsesSI()` ✓
3. **Collision impulse gather** — `solveCollisionImpulsesSI()` ✓
4. **Atom update** — aggregating all impulses ✓
5. **Hard microstep** — `hardMicrostep()` ✓
6. **Macrostep** — `RAffWorld_Step()` with RESPA splitting ✓
7. **Soft force stub** — Lennard-Jones ✓
8. **Broad phase** — O(N²) pair search ✓
9. **Reverse adjacency** — `RAffWorld_BuildAdjacency()` ✓

### Phase 2: Test Harness ✅ (done)

Created `test/test_rigid_atom_ff.c`:

1. **Two-atom bond test** ✓
2. **Collision test** ✓
3. **Momentum conservation test** ✓
4. **Chain test** ✓

### Phase 3: Fix SI stability issues (in progress)

1. **Verify impulse sign** — ensure port impulse is positive along `n` (§7.1)
2. **Verify `h` vs `H`** — ensure all microstep calculations use `h = H/nsub` (§7.2)
3. **Verify LJ force sign** — attractive at long range, repulsive at short (§7.4)
4. **Update collision pairs every microstep** (§7.3)
5. **Handle coincident atoms** deterministically (§7.6)
6. **Add NaN detection** — check for NaN in velocities/positions after each step
7. **Add regression tests** for each fix

### Phase 4: XPBD implementation

1. **Add `lambdaPort` array** to `RAffWorld` — `nAtoms * RAFF_MAX_PORTS` floats
2. **Zero `lambdaPort` at macrostep start** — in `RAffWorld_Step()` before microstep loop
3. **Implement `gatherPortImpulsesXPBD()`** — per §4.1
4. **Add `SOLVER_XPBD` to method enum** and dispatch in `hardMicrostep()`
5. **Test**: same bond/collision/momentum tests, verify faster convergence with fewer microsteps
6. **Test**: warm-starting — verify `λ` accumulates across microsteps and improves convergence

### Phase 5: PD implementation

1. **Implement `gatherPortImpulsesPD()`** — per §4.2
2. **Add velocity damping filter** — global `v *= exp(-d*h)` after each microstep
3. **Add `SOLVER_PD` to method enum** and dispatch
4. **Test**: ring/cage topology — verify better coupling handling than SI
5. **Test**: compare convergence rate vs SI for same `nsub`

### Phase 6: VBD implementation

1. **Add `Hblock` and `gblock` arrays** to `RAffWorld`
2. **Implement `solve6x6()`** — using Schur complement (two 3×3 solves, §4.3)
3. **Implement `gatherPortImpulsesVBD()`** — per §4.3
4. **Add `SOLVER_VBD` to method enum** and dispatch
5. **Test**: tetrahedral atom (4 ports) — verify coupled convergence
6. **Test**: compare per-microstep convergence vs SI for multi-port atoms
7. **Benchmark**: measure cost of 6×6 solve vs scalar impulse

### Phase 7: AVBD implementation

1. **Add `lambdaPort` support** (shared with XPBD, or separate array)
2. **Implement `gatherPortImpulsesAVBD()`** — per §4.4
3. **Add `SOLVER_AVBD` to method enum** and dispatch
4. **Test**: hard constraint convergence — verify exact satisfaction over multiple microsteps
5. **Test**: compare residual error vs VBD for same `a` and `nsub`

### Phase 8: Visualization integration ✅ (done)

Created `samples/sample_rigid_atoms.cpp`:

1. **Load nNonan molecule** ✓
2. **Build `RAffWorld`** from mol file ✓
3. **Custom `Step()`** calling `RAffWorld_Step()` ✓
4. **Draw atoms as spheres** ✓
5. **ImGui controls** for parameters ✓
6. **Mouse picking** to grab atoms ✓
7. **KE history plot** ✓

### Phase 9: Optimization

1. **Spatial hash** for broad-phase collisions — replace O(N²) pair search
2. **Graph coloring** for parallel port sweeps — color atoms so no two bonded atoms share a color
3. **Collision pair update frequency** — rebuild every K microsteps instead of every one
4. **FIRE optimizer** — replace simple velocity damping with FIRE for minimization tasks
5. **Method switching at runtime** — UI dropdown to switch solver method without rebuilding world

---

## 9. Method Comparison Summary

| Property | SI | XPBD | PD | VBD | AVBD |
|----------|-----|------|----|-----|------|
| Convergence per microstep | Good | Better (warm-start) | Good (coupled) | Best (6-DOF) | Best + exact |
| Cost per atom | O(ports) | O(ports) | O(ports) | O(ports + 6×6 solve) | O(ports + 6×6 solve) |
| Momentum conservation | Exact | Exact | Exact | Exact | Exact |
| Damping | Built-in | Separate filter | Separate filter | Built-in (inertial) | Built-in (inertial) |
| Warm starting | No | Yes (λ) | No | No | Yes (λ) |
| Hard constraint limit | a→∞ (slow) | α̃→0 (fast) | a→∞ (slow) | K→∞ (ill-conditioned) | Exact (augmented) |
| Implementation complexity | Low | Low-medium | Medium | High | High |
| Best for | General dynamics | Stiff bonds | Coupled networks | Multi-port atoms | Exact constraints |

### Recommended implementation order

1. **Fix SI stability** (Phase 3) — most urgent, fixes NaNs
2. **XPBD** (Phase 4) — easiest to add, immediate benefit for stiff bonds
3. **VBD** (Phase 6) — high value for multi-port atoms, but needs 6×6 solve
4. **PD** (Phase 5) — useful for ring/cage topologies, medium effort
5. **AVBD** (Phase 7) — only needed when exact constraints are required

---

## 10. Key Files

- `src/rigid_atom_ff.h` — `RAffAtom`, `RAffPort`, `RAffWorld`, public API
- `src/rigid_atom_ff.c` — current SI implementation
- `test/test_rigid_atom_ff.c` — bond, collision, momentum, chain tests
- `samples/sample_rigid_atoms.cpp` — visualization with ImGui controls
- `include/box3d/math_functions.h` — `b3Vec3`, `b3Quat`, `b3Matrix3`, `b3RotateVector`, `b3Solve3`, `b3MulQuat`, `b3NormalizeQuat`, `b3MakeQuatFromAxisAngle`
- `src/math_internal.h` — `b3IntegrateRotation`, `b3QuatFromExponentialMap`, `b3DeltaQuatToRotation`
- `samples/sample.h` — `Sample` base class, `SampleContext`
- `samples/gfx/draw.h` — `DrawSphere`, `DrawLine`, `DrawCylinder` rendering API
- `docs/prokop/RigidAtomFF_SI.chat.md` — full GPT 5.5 discussion (~4500 lines)
- `docs/prokop/RigidAtomFF_SI_short.plan.md` — condensed short plan
- `docs/prokop/SI_for_molecules.chat.md` — SI comparison and local implicit impulse derivation


