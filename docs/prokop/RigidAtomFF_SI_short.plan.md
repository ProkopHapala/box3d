# RigidAtomFF SI Implementation Plan

## Architecture

Three-level time-stepping (RESPA/Strang splitting):

```
soft macrostep H:
  evaluate soft forces O(N^2) once
  for s hard microsteps (h=H/s):
    frozen soft kick h
    gather port + collision impulses per atom
    apply one velocity/rotation update per atom
    drift by h
  new soft half-kick
```

- **Atom-centered**: state is only `pos, quat, vel, omega` per atom
- **Local implicit harmonic impulse** for bonds/ports/collisions — no global solve
- **No warm starting**, no persistent multipliers
- **Scalar radial reactions** along normals → exact momentum conservation

## Data Structures (`src/rigid_atom_ff.h`)

```c
#define RAFF_MAX_PORTS 4

typedef struct {
    b3Vec3 pos; b3Quat quat; b3Vec3 vel; b3Vec3 omega;
    float mass, invMass; b3Matrix3 invInertia; float radius;
    int nPorts;
} RAtom;

typedef struct {
    int neighIdx; b3Vec3 localDir; float aBond; float dBond;
} RPort;

typedef struct {
    RAtom* atoms; int nAtoms;
    RPort* ports;               // flat nAtoms*RAFF_MAX_PORTS
    int* incomingCount; int* incomingOffset;
    int* incomingAtoms; int* incomingPorts;  // reverse adjacency
    // collision pairs
    int* collI; int* collJ; float* collR; int nColl;
    // soft forces (cached)
    b3Vec3* fsoft; b3Vec3* tsoft;
    // temp accumulators
    b3Vec3* Jown; b3Vec3* Lown; b3Vec3* Jrecoil; b3Vec3* Jcoll;
    // params
    float H; int nsub;
    float aBond, dBond, aAngle, dAngle, aColl, dColl;
} RAtomWorld;
```

Uses Box3D math: `b3RotateVector`, `b3MulQuat`, `b3NormalizeQuat`, `b3IntegrateRotation`, `b3Solve3`, `b3Cross`, `b3Dot` from `@/home/prokop/git/box3d/include/box3d/math_functions.h` and `@/home/prokop/git/box3d/src/math_internal.h`.

## Core Functions (`src/rigid_atom_ff.c`)

### Macrostep

```c
void RAtomWorld_Step(RAtomWorld* w) {
    applySoftKick(w, 0.5f*w->H, w->fsoft, w->tsoft);
    float h = w->H / w->nsub;
    for (int s=0; s<w->nsub; s++) {
        for (int i=0; i<w->nAtoms; i++) {
            w->atoms[i].vel = b3Add(w->atoms[i].vel,
                b3MulSV(h*w->atoms[i].invMass, w->fsoft[i]));
            w->atoms[i].omega = b3Add(w->atoms[i].omega,
                b3MulMV(w->atoms[i].invInertia, b3MulSV(h, w->tsoft[i])));
        }
        hardMicrostep(w, h);
    }
    evaluateSoftForces(w);
    applySoftKick(w, 0.5f*w->H, w->fsoft, w->tsoft);
}
```

### Port Impulse

```c
// For port k of atom i toward neighbor j:
// r = rotate(quat_i, localDir_k)
// e = pos_j - (pos_i + r);  C = |e|;  n = e/C
// u = n . (v_j - (v_i + omega_i x r))
// s = r x n;  w = invM_i + invM_j + s^T invI_i s
// a = aBond; d = dBond (dimensionless)
// p = -(a*C/h + (a+d)*u) / (w*(1+a+d))
// P = p*n
// Jown[i] += P;  Lown[i] += r x P;  Jrecoil[i*MAX+k] = -P
```

### Collision Impulse

```c
// For pair (i,j): delta = Rij - |pos_i - pos_j|;  if delta<=0 skip
// n = (pos_i - pos_j)/|...|;  uc = n.(v_i - v_j)
// wc = invM_i + invM_j
// p = max(0, (a*delta/h - (a+d)*uc) / (wc*(1+a+d)))
// Jcoll[i] += p*n;  Jcoll[j] -= p*n
```

### Atom Update

```c
// J = Jown[i] + sum(incoming Jrecoil) + Jcoll[i]
// L = Lown[i]
// vel += invM * J;  omega += invI * L
// pos += h*vel;  quat = normalize(integrateRotation(quat, h*omega))
```

## Momentum Conservation

- **Ports**: `J_i + J_j = P - P = 0`. Angular: `(tip - x_j) x P = -e x (p*n) = 0` since `P || e`.
- **Collisions**: `P_c - P_c = 0`. Angular: `d x (p*n) = 0` since `P || d`.
- No per-atom clipping. Same relaxation on both endpoints. Radial damping only.

## Parameters (dimensionless, not physical)

| Param | Meaning | Start |
|-------|---------|-------|
| aBond | bond correction fraction | 10 |
| dBond | bond damping | 1-3 |
| aAngle | port/angle correction | 3-10 |
| dAngle | port/angle damping | 1-3 |
| aColl | collision correction | 5-20 |
| dColl | collision damping | 2-5 |
| H | soft macrostep | tune |
| nsub | hard microsteps | 4-16 |

`a=1`→50% error removed, `a=9`→90%, `a→∞`→hard constraint. `rho=1/(1+a)`.

## Implementation Phases

### Phase 1: Core solver — `src/rigid_atom_ff.h` + `src/rigid_atom_ff.c`
1. Data structures + init/destroy
2. Reverse adjacency precompute from port topology
3. `gatherPortImpulses` — local implicit formula
4. `gatherCollisionImpulses` — unilateral with clamp
5. `applyAtomUpdate` — aggregate all impulses, integrate
6. `hardMicrostep` — orchestrate gather+apply
7. `RAtomWorld_Step` — RESPA macrostep
8. Soft force stub (O(N^2) Lennard-Jones)
9. Broad phase: simple O(N^2) pair search; upgrade to spatial hash later

### Phase 2: Tests — `test/test_rigid_atom_ff.c`
1. Two-atom bond: length converges, momentum conserved
2. Three-atom angle: angle converges, momentum conserved
3. Collision: two atoms approach, verify repulsion + momentum
4. Chain relaxation: small polymer, stable with large H
5. Momentum check: `sum(p)=0`, `sum(L)=0` after N steps

### Phase 3: Visualization — `samples/sample_rigid_atom_ff.cpp`
1. Inherit `Sample`, follow `sample_molecules.cpp` pattern
2. Load `.mol` file via `b3LoadMolFile`
3. Build `RAtomWorld` from atoms/bonds
4. Custom `Step()` calling `RAtomWorld_Step()` (not `b3World_Step`)
5. Draw: `DrawSphere` per atom, `DrawLine`/`DrawCylinder` per bond
6. ImGui sliders for H, nsub, aBond, aColl, damping
7. Optional debug: port tips, collision normals, velocity arrows

### Phase 4: Optimization
1. Spatial hash for broad-phase collisions
2. Graph coloring for parallel port sweeps
3. Collision pair update every K microsteps (not every one)
4. Optional: 6x6 VBD block solve for coupled atoms (upgrade path)

## Key Files

- `@/home/prokop/git/box3d/include/box3d/math_functions.h` — b3Vec3, b3Quat, b3Matrix3, b3RotateVector, b3Solve3
- `@/home/prokop/git/box3d/src/math_internal.h` — b3IntegrateRotation, b3QuatFromExponentialMap
- `@/home/prokop/git/box3d/samples/sample.h` — Sample base class
- `@/home/prokop/git/box3d/samples/sample_molecules.cpp` — reference for mol loading + visualization
- `@/home/prokop/git/box3d/docs/prokop/RigidAtomFF_SI.chat.md` — full discussion (4067 lines)
