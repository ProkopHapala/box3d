# USER

OK, how would you implement simple demo of molecular dynamics using this engine?

I mean molecules are atoms connected by sticks (distance contrasn but also free rotation), but there are angular constrains between the sticks, and most importantly the particles have collisions they are sphere collision objects.

Or do you think it would be better to ad angular constrain, something like "as rigid as possible" within the sequantial impulses framework ?

---

# GLM 5.2

Good question. Let me look at how joints are implemented to give concrete advice.

Great question. There's a clean mapping from MD to Box3D that works entirely within the existing API — no engine modifications needed. The key insight comes from computational chemistry itself.

## The SHAKE Approach: Angles as Distances

In MD, the **SHAKE** algorithm (Ryckaert et al., 1977) constrains bond angles by converting them to **1-3 distance constraints**. If atoms A-B-C have fixed bond lengths `L_AB` and `L_BC`, then the angle θ₀ at B uniquely determines the distance between A and C:

```
|AC|² = L_AB² + L_BC² - 2·L_AB·L_BC·cos(θ₀)
|AC|  = sqrt(L_AB² + L_BC² - 2·L_AB·L_BC·cos(θ₀))
```

So instead of a 3-body angular constraint, you add a **2-body distance joint between A and C** with the computed length. Box3D's existing distance joint handles this perfectly. This is exactly what SHAKE does — and it maps cleanly onto Box3D's 2-body joint framework.

## Concrete Implementation Plan

### 1. Atoms → Dynamic Bodies with Sphere Shapes

```c
b3BodyDef bodyDef = b3DefaultBodyDef();
bodyDef.type = b3_dynamicBody;
bodyDef.position = {x, y, z};

b3BodyId atomId = b3CreateBody(worldId, bodyDef);

b3ShapeDef shapeDef = b3DefaultShapeDef();
shapeDef.density = 1.0f;  // atomic mass / volume
shapeDef.material.friction = 0.1f;
shapeDef.material.restitution = 0.0f;

b3Sphere sphere = {.radius = atomRadius};
b3CreateShapeSphere(atomId, sphere, shapeDef);
```

Sphere-sphere collisions are automatic — non-bonded atoms will collide and push each other apart. This gives you **excluded volume** (the repulsive part of Lennard-Jones) for free.

### 2. Bonds (1-2) → Distance Joints

```c
b3DistanceJointDef bondDef = b3DefaultDistanceJointDef();
bondDef.bodyIdA = atomA;
bondDef.bodyIdB = atomB;
bondDef.localFrameA.p = {0, 0, 0};  // atom center
bondDef.localFrameB.p = {0, 0, 0};
bondDef.length = bondLength;          // e.g. 1.0 Å → scaled to engine units
bondDef.collideConnected = false;     // bonded atoms don't collide
bondDef.hertz = 0.0f;                 // rigid (no spring)
bondDef.dampingRatio = 0.0f;

b3JointId bondId = b3CreateDistanceJoint(worldId, bondDef);
```

The distance joint maintains the bond length exactly (within solver tolerance). Free rotation is the default — the joint only constrains distance, not orientation. This is a **ball-and-stick bond**.

### 3. Bond Angles (1-3) → Distance Joints Between A and C

```c
float compute_1_3_distance(float lab, float lbc, float angle_deg) {
    float theta = angle_deg * B3_PI / 180.0f;
    return sqrtf(lab*lab + lbc*lbc - 2.0f*lab*lbc*cosf(theta));
}

// For angle A-B-C = 109.5° (tetrahedral), with L_AB = L_BC = 1.0:
float d13 = compute_1_3_distance(1.0f, 1.0f, 109.5f);  // = 1.633

b3DistanceJointDef angleDef = b3DefaultDistanceJointDef();
angleDef.bodyIdA = atomA;
angleDef.bodyIdB = atomC;  // note: A and C, not B
angleDef.localFrameA.p = {0, 0, 0};
angleDef.localFrameB.p = {0, 0, 0};
angleDef.length = d13;
angleDef.collideConnected = false;  // 1-3 atoms shouldn't collide either
angleDef.hertz = 0.0f;
angleDef.dampingRatio = 0.0f;

b3JointId angleJointId = b3CreateDistanceJoint(worldId, angleDef);
```

This constrains the angle indirectly. Given fixed bond lengths (enforced by the 1-2 joints), the 1-3 distance uniquely determines the angle. The solver handles all of these simultaneously via sequential impulses.

### 4. Dihedral Angles (1-4) → Optional Distance Joints

For torsion angles A-B-C-D, you can add a 1-4 distance joint between A and D. The distance depends on the dihedral angle φ, the bond lengths, and the bond angles:

```
|AD| = f(L_AB, L_BC, L_CD, θ_ABC, θ_BCD, φ)
```

This is more complex to compute but follows the same principle. For a simple demo, you might skip dihedrals or use soft springs (enable spring on the distance joint with a low hertz) to allow torsion.

### 5. Non-Bonded Interactions → Automatic Sphere Collisions

Any pair of atoms not connected by a joint (or with `collideConnected = true`) will collide via sphere-sphere contact. This gives you:
- **Repulsion** (hard sphere) — the normal contact impulse
- **Friction** — tangential impulse (set low for atoms)

For **attractive** interactions (Lennard-Jones well, hydrogen bonds), Box3D has no built-in force field. You'd apply forces manually each step:

```c
void Step() {
    // Apply LJ attraction between non-bonded pairs
    for (int i = 0; i < atomCount; i++) {
        for (int j = i+1; j < atomCount; j++) {
            if (is_bonded(i, j) || is_1_3(i, j)) continue;
            
            b3Vec3 ri = b3Body_GetPosition(atoms[i]);
            b3Vec3 rj = b3Body_GetPosition(atoms[j]);
            b3Vec3 d = b3Sub(rj, ri);
            float r = b3Length(d);
            if (r > cutoff) continue;
            
            // F = -dU/dr * r_hat, U = 4ε[(σ/r)^12 - (σ/r)^6]
            float sr6 = powf(sigma/r, 6.0f);
            float f = 24.0f * epsilon * (2.0f*sr6*sr6 - sr6) / (r*r);
            b3Vec3 force = b3MulSV(f, d);
            b3Body_ApplyForceToCenter(atoms[i], force, false);
            b3Body_ApplyForceToCenter(atoms[j], b3Neg(force), false);
        }
    }
    
    b3World_Step(worldId, timeStep, subStepCount);
}
```

Note: this force application is **explicit** (MD-style, Jacobi), while the bond/angle constraints are **implicit** (SI). This hybrid works because the LJ forces are soft (low stiffness) and don't need tiny time steps, while the hard constraints (bonds, angles) are handled implicitly by the solver.

## What About "As Rigid as Possible" Within SI?

Your question about adding a proper angular constraint within the SI framework is interesting. There are two interpretations:

### Option A: 3-Body Angular Constraint (New Joint Type)

You'd add a new joint type that takes 3 bodies and constrains the angle A-B-C directly. The velocity constraint would be:

```
C = arccos(BA·BC / |BA||BC|) - θ₀ = 0
Ċ = dC/dt → J_A·v_A + J_B·v_B + J_C·v_C = 0
```

The Jacobian has terms for all three bodies. This is **mathematically cleaner** than the 1-3 distance trick, but requires:
- Modifying `b3JointSim` to support 3 bodies (currently hardcoded to 2)
- Modifying the constraint graph (currently assumes 2-body joints)
- Modifying the solver prepare/warm-start/solve functions
- Adding new prepare/solve functions

This is a **significant engine modification**. The joint system, constraint graph, and solver are all built around the 2-body assumption.

### Option B: 1-3 Distance (SHAKE-style, No Engine Changes)

As described above. The 1-3 distance constraint is **mathematically equivalent** to the angular constraint when bond lengths are fixed. The only difference is that the 1-3 distance also constrains the bond lengths (redundantly), which actually helps stability.

**My recommendation: Use Option B (1-3 distances).** It's the standard approach in constraint-based MD (SHAKE/LINCS), requires zero engine modifications, and the solver already handles it optimally. The "as rigid as possible" formulation would give the same result with much more implementation effort.

### Option C: Soft Angular Springs (MD-Style)

Instead of hard constraints, use the motor joint's angular spring to apply a torque proportional to angle deviation. This is the "harmonic angle potential" approach from MD:

```
U_angle = 1/2 * k_θ * (θ - θ₀)²
τ = -k_θ * (θ - θ₀) * dθ/d(orientation)
```

You'd apply this as a torque on body B each step. This is **explicit** (MD-style), so it needs `dt < 2/√(k_θ/I)` for stability. But it allows the angle to fluctuate (thermally), which is physically correct for MD. Use this if you want **soft, thermalizable** angles rather than rigid constraints.

## Summary: Which Approach?

| Approach | Rigidity | Engine Changes | Stability | Thermal Fluctuation |
|---|---|---|---|---|
| 1-3 distance joints (SHAKE) | Rigid | None | Excellent (implicit) | No |
| 3-body angular joint | Rigid | Major | Excellent (implicit) | No |
| Soft angular springs (MD) | Soft | None | Conditional (dt-limited) | Yes |
| Motor joint angular spring | Soft | None | Good (implicit soft) | Yes |

**For a demo**: Use 1-3 distance joints (SHAKE-style). It's the cleanest, requires no engine changes, and maps directly to how real constraint-based MD works.

**For realistic MD with thermal fluctuations**: Use soft angular springs applied as torques each step, combined with distance joints for bonds. The bonds are rigid (implicit), angles are soft (explicit spring). This is actually how many MD packages work — rigid bonds, flexible angles.

**For "as rigid as possible" within SI**: The 1-3 distance approach IS "as rigid as possible" within SI — it's a hard constraint solved by sequential impulses. A dedicated 3-body angular joint would give the same result with far more implementation work.

Want me to write a concrete sample implementation as a Box3D sample file?

---

# USER


No by "As Rigid as possible" I meant to make ports
like here
/home/prokophapala/git/FireCore/pyBall/RigidAtomFF/RRsp3/RRsp3.cl

can we implement this easily?

---

# GLM 5.2

Now I understand your code. This is a fundamentally different and more elegant approach than SHAKE/distance joints. Let me explain the mapping.

## What RRsp3 Does: Rigid Ports

Your code models each atom as a **rigid body with attachment points (ports)**. Each port is a local offset `r_local` that rotates with the atom's quaternion. The constraint is:

```
tip = pos_i + R(q_i) · r_local    ← world-space port tip position
constraint: tip = pos_j            ← tip must coincide with neighbor's center
```

This simultaneously constrains **distance** (the tip-to-center distance) and **angle** (the port direction in the body frame). When a node has multiple ports, it must rotate to align *all* ports toward their neighbors simultaneously — this is the "as rigid as possible" behavior. The constraint generates both a linear impulse (`dpos`) and a torque (`dtheta = cross(r_arm, P) * invI`) because the port is offset from the center of mass.

The XPBD impulse formula in your code (line 764-769):
```c
float3 rxn = cross(r_arm, n);
float w_ang = dot(rxn, rxn) * invI;        // angular effective mass
float alpha = 1.0f / (K * dt2);             // XPBD compliance
float w_total = w_i + w_j + w_ang + alpha;
float impulse_mag = dist / w_total;
```

This is exactly a **point-to-point constraint with a lever arm** — the constraint point is at the port tip, not at the body's center of mass.

## The Mapping: Spherical Joints = Ports

Box3D's **spherical joint** (ball-and-socket) does exactly this. It constrains a point on body A (defined in A's local frame) to coincide with a point on body B. The constraint generates a linear impulse at the anchor point, which — because the anchor is offset from the center of mass — produces both linear and angular corrections. This is the port constraint, verbatim.

```c
// For each port on a node atom:
b3SphericalJointDef portDef = b3DefaultSphericalJointDef();
portDef.bodyIdA = nodeAtomId;           // the node (has orientation)
portDef.bodyIdB = neighborAtomId;       // the neighbor (cap or node)
portDef.localFrameA.p = r_local;        // port position in node's local frame
portDef.localFrameB.p = {0, 0, 0};      // neighbor's center
portDef.collideConnected = false;        // bonded atoms don't collide

// For rigid ports (hard constraint):
// (spring disabled by default → rigid)

// For soft ports matching RRsp3 stiffness K:
portDef.enableSpring = true;
portDef.hertz = sqrtf(K / m) / (2.0f * B3_PI);  // convert stiffness to Hz
portDef.dampingRatio = 0.1f;

b3JointId portJointId = b3CreateSphericalJoint(worldId, portDef);
```

When the node has 4 ports, you create 4 spherical joints. The SI solver iterates over all 4, and the node rotates to find the best compromise alignment — **this is "as rigid as possible" within the sequential impulse framework**. No engine modifications needed.

## How Each RRsp3 Variant Maps to Box3D

| RRsp3 variant | What it does | Box3D equivalent |
|---|---|---|
| **Massfull XPBD** (lines 659-784) | XPBD with physical inertia, linear+angular recoil | **Spherical joint with spring** — same math, same physics. [b3MakeSoft](cci:1://file:///home/prokophapala/git_SW/box3d/src/solver.h:263:0-305:1) computes the same bias/mass/impulse scales. |
| **Massless Newton-Raphson** (lines 916-1084) | Solves `H·ω = torque` directly for optimal rotation | **No equivalent** — Box3D's SI iterates to the same solution but doesn't solve the rotation directly. Converges slower but handles dynamic topology. |
| **Massless Shape Match / Kabsch** (lines 1089-1255) | Polar decomposition for best-fit rotation | **No equivalent** — same as above. SI iteration is the iterative approximation. |
| **Massless Eigen / Davenport q-method** (lines 1260-1354) | 4×4 eigenproblem for optimal quaternion | **No equivalent** — but SI with enough sub-steps converges to the same answer. |

The key difference: RRsp3's massless variants solve for the optimal rotation **directly** (one-shot or few iterations), which is faster per iteration but assumes fixed topology (neighbor list doesn't change). Box3D's SI approaches the same solution **iteratively** through repeated impulse application, which is slower per iteration but handles dynamic topology (contacts appear/disappear) naturally.

## What Box3D Gives You For Free

Things RRsp3 implements as separate kernels, Box3D handles natively:

1. **Sphere-sphere collisions** — `compute_collision_cluster_rigid` kernel → Box3D's broad phase + narrow phase + contact solver. Automatic, with spatial acceleration (dynamic AABB tree), exclusion lists via `collideConnected = false`, and speculative contacts for fast-moving atoms.

2. **Parallelism** — workgroup clustering + ghost atoms → Box3D's graph coloring + per-block atomic CAS. Different approach, same goal: partition constraints into independent groups for multi-core parallelism. Box3D uses Gauss-Seidel (better convergence) instead of Jacobi.

3. **Integration** — `predict_dynamics` / `update_velocities_dynamics` → Box3D's `b3IntegrateVelocitiesTask` / `b3IntegratePositionsTask`. Box3D also handles gyroscopic torque (important for non-spherical atoms).

4. **Exclusion lists** — `excl1_global` / `excl2_global` → `collideConnected = false` on joints + collision filters on shapes.

5. **Heavy-ball momentum** — `dpos_mom * beta` in `apply_corrections_rigid_ports` → Box3D's warm starting serves a similar purpose (reusing previous impulses).

## What You'd Need to Implement Manually

1. **Lennard-Jones / attractive forces** — apply per-step via `b3Body_ApplyForceToCenter()` for non-bonded pairs. This is explicit (MD-style), same as your `predict_dynamics` force application.

2. **Damping** — `vel *= damp` in `update_velocities_dynamics` → Box3D's `linearDamping` / `angularDamping` on the body definition.

3. **Fixed atoms** — `fixmask` → set body type to `b3_staticBody` or use `b3Body_SetType()`.

4. **Custom port stiffness per bond** — set `hertz` per spherical joint. RRsp3 uses `stiffness_flat[idx]`; Box3D uses per-joint spring parameters.

## Concrete Implementation Sketch

```c
// 1. Create atoms as dynamic bodies with sphere shapes
for (int i = 0; i < natoms; i++) {
    b3BodyDef def = b3DefaultBodyDef();
    def.type = (fixmask[i] & 7) ? b3_staticBody : b3_dynamicBody;
    def.position = {pos[i].x, pos[i].y, pos[i].z};
    def.linearDamping = damp;
    def.angularDamping = damp;
    
    b3BodyId atomId = b3CreateBody(worldId, def);
    
    b3ShapeDef shapeDef = b3DefaultShapeDef();
    shapeDef.density = 1.0f / (pos[i].w * volume);  // invMass → mass
    shapeDef.material.friction = 0.1f;
    
    b3Sphere sphere = {.radius = radius[i]};
    b3CreateShapeSphere(atomId, sphere, shapeDef);
    
    atoms[i] = atomId;
}

// 2. Create port constraints as spherical joints
for (int inode = 0; inode < nnodes; inode++) {
    for (int k = 0; k < 4; k++) {
        int idx = inode * 4 + k;
        int j = neighbors[idx];
        if (j < 0 || stiffness[idx] <= 0.0f) continue;
        
        b3SphericalJointDef portDef = b3DefaultSphericalJointDef();
        portDef.bodyIdA = atoms[inode];           // node
        portDef.bodyIdB = atoms[j];               // neighbor
        portDef.localFrameA.p = {port_local[idx].x, port_local[idx].y, port_local[idx].z};
        portDef.localFrameB.p = {0, 0, 0};
        portDef.collideConnected = false;
        
        // Soft constraint matching RRsp3 stiffness K
        float K = stiffness[idx];
        float m = 1.0f / pos[inode].w;
        portDef.enableSpring = true;
        portDef.hertz = sqrtf(K / m) / (2.0f * B3_PI);
        portDef.dampingRatio = 0.1f;
        
        b3CreateSphericalJoint(worldId, portDef);
    }
}

// 3. Step with sub-stepping (replaces RRsp3 iteration loop)
b3World_Step(worldId, dt, subStepCount);
// subStepCount ≈ n_rot_substeps in your code
```

## The Fundamental Difference

Your RRsp3 code solves the rotation **per atom locally** — each workgroup processes one cluster, and the rotation solve (Newton-Raphson / Kabsch / eigen) finds the best orientation for that node given its neighbors' positions. This is a **local solve** — it doesn't propagate information to neighbors within the same iteration. Information propagates across iterations (Jacobi-style).

Box3D's SI solves constraints **globally across the island** — when a port constraint on atom A is solved, the impulse changes A's velocity, which is immediately visible to the next constraint involving A (Gauss-Seidel). Information propagates within a single iteration. This means SI typically needs **fewer iterations** than RRsp3's Jacobi-style approach, but each iteration is sequential (within a graph color).

For a molecule with ~100 atoms, the difference is:
- RRsp3: ~10-20 Jacobi iterations, each fully parallel, each doing a direct rotation solve
- Box3D: ~4 sub-steps × 1 SI iteration, each sub-step sequential within colors but graph-color parallel

Both converge to the same rigid configuration. Box3D's approach is simpler to implement (no custom rotation solver) and handles collisions natively, but doesn't have the per-atom direct rotation solve that makes RRsp3's massless variants efficient on GPU.

**Bottom line: Yes, this is easy to implement in Box3D using spherical joints. No engine modifications needed. The "as rigid as possible" behavior emerges from SI iteration over multiple port constraints on the same body.**

---

# USER

OK, so now implement this basic principle in this engine, make simple new demo where we can have n=10 molecules of nNonan.mol

---

# GLM 5.2

## Implementation: Done

A working demo has been created. It loads nNonane (C9H20, 29 atoms, 28 bonds) from an MDL V2000 MOL file and spawns 10 molecules in a grid above the ground. The molecules fall under gravity, collide with each other and the ground, and the rigid port constraints maintain molecular geometry.

### Key Insight: Two Joints Per Bond

The initial implementation used one spherical joint per bond (port on A toward B). This constrained the distance but **not the orientation** — atom B had `localFrameB.p = {0,0,0}`, so `rB = 0` and `cross(rB, impulse) = 0`, meaning B received no angular correction from that joint.

The fix: **create two spherical joints per bond**:
1. Port on A toward B: `localFrameA.p = B_init - A_init`, `localFrameB.p = {0,0,0}`
2. Port on B toward A: `localFrameA.p = A_init - B_init`, `localFrameB.p = {0,0,0}`

Now each atom has a port toward every neighbor. A carbon with 4 bonds has 4 port constraints. The SI solver iterates over all 4, and the torques from each port compete to find the orientation that best aligns all ports simultaneously — the "as rigid as possible" behavior.

What's still free: rotation around a single bond axis (torsion). An atom with only 1 bond (like H) can spin around that bond. But for atoms with 2+ bonds, the multiple port directions over-constrain the orientation, effectively locking it. This matches real molecular behavior.

### How It Maps to the Solver

Looking at `src/spherical_joint.c`, the point-to-point constraint (lines 604-645):

```c
// rA = rotated localFrameA.p (the port offset from A's center)
// rB = rotated localFrameB.p (= 0, since B's port is at its center)
b3Vec3 cdot = b3Sub(b3Sub(b3Add(vB, cross(wB, rB)), vA), cross(wA, rA));

// K matrix includes cross(rA, invIA, rA) — the angular effective mass from A's port
b3Matrix3 kA = b3MulMM(sA, b3MulMM(base->invIA, sA));

// Impulse applied:
vA -= mA * impulse;
wA -= iA * cross(rA, impulse);   // ← torque on A from port offset
vB += mB * impulse;
wB += iB * cross(rB, impulse);   // ← zero when rB = 0
```

With two joints per bond, both A and B get `cross(r, impulse)` torque terms. The SI solver iterates over all joints on a body, so multiple ports on the same atom produce competing torques that converge to the best-fit orientation.

### Demo Location

- **Sample file**: `samples/sample_molecules.cpp` — registered as **Molecules > Rigid Ports** in the sample picker
- **MOL parser**: `samples/mol_loader.h` — minimal MDL V2000 parser (atoms, bonds, element symbols)
- **Data file**: `data/mol/nNonan.mol` — copy of the nNonane molecule
- **Build**: added `sample_molecules.cpp` to `samples/CMakeLists.txt`

### Running

```bash
cd /home/prokophapala/git_SW/box3d
./build/bin/samples
```

Select **Molecules > Rigid Ports** from the sample picker (Ctrl+O or menu). 10 nNonane molecules spawn in a grid and fall under gravity. The UI panel shows atom/bond counts and has a Soft Spring toggle and Respawn button.

### What Box3D Handles For Free

- Sphere-sphere collisions between non-bonded atoms (excluded volume)
- Broad-phase spatial acceleration (dynamic AABB tree)
- Graph coloring for parallel constraint solving
- Gyroscopic torque for non-spherical inertia
- Warm starting (reuses previous impulses, similar to RRsp3's heavy-ball momentum)
- Sleep/wake for performance