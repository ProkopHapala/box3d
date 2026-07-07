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

---

# USER

now look here @overview.md I'm interested about the methods of Sequential impulses, if it can be used for relaxation of molecules. Currently I use projective dynamics or position based dynamics to solve problem of time-step being limited by hard degrees of freedom, while soft degrees of freedom then relax slow, so we wnat some spliting schemem, where we solve hard dofs separately (relax them) and then we propagate soft-degrees with long time step. (hard are bond lenghs, angles, hard sphere colisions; soft are dihedrals angles, atractive coulomb and vand der walls interactions; we consider all hard DOFs to be harmonic localy explresed by some stiffness matrix, althou collision are inequality constrains, they are still harmonic in penetration r<R)

now if sequantal impulses is more stable and more efficient than XPBD or Projective synamics for ridigs body system ,it may be also faster for molecular cynamics? I do not care too much about accuracy (conservation of energy etc.) I want model which is qualitatively fine (e.g. proprly conserve momentum both linear and rotational) and super stable and fast, to relax system as fast as possible close to minimum (have reasonable estimate of minimu menergy structure), I even replace Lenard-Jones by some simplified polynominak and hard spehre collisons to improve the speed, so do not care about accuracy too much f it does not break momentum conservation and other fundamental inveriants.

---

# GLM 5.2

Let me look at the Box3D solver implementation to understand how Sequential Impulses works in this codebase, then I can give you an informed comparison.

Now let me look at the actual constraint math — the distance joint (closest analog to bonds) and the core solve loop:



Let me also check the soft constraint formulation and the contact solver's inequality handling:




Now I have a thorough understanding of the Box3D solver. Here's my analysis:

---

## Box3D's Sequential Impulses vs XPBD/PD for Molecular Relaxation

### How Box3D's Solver Works

The solver pipeline per substep (see `@/home/prokop/git/box3d/src/solver.c:1228-1305`) is:

1. **Integrate velocities** — apply forces, gravity, damping, gyroscopic torque
2. **Warm start** — re-apply accumulated impulses from previous step
3. **Solve** (ITERATIONS=1) — iterate over all constraints by graph color, with bias (position correction)
4. **Integrate positions** — `x += h * v`
5. **Relax** (RELAX_ITERATIONS=1) — iterate again *without* bias, just velocity-level constraint solving

Each constraint (joint or contact) computes:
- **Effective mass** `K = J · M⁻¹ · Jᵀ` (scalar per DOF for distance, 3×3 for angular)
- **Impulse**: `λ = -K⁻¹ · (Cdot + bias) - impulseScale · λ_old`
- **Velocity update**: `v += M⁻¹ · Jᵀ · λ`

For **soft constraints** (springs), the [b3MakeSoft](cci:1://file:///home/prokop/git/box3d/src/solver.h:263:0-305:1) function (`@/home/prokop/git/box3d/src/solver.h:264-306`) maps hertz + damping ratio to a mass-scale/impulse-scale split:
- `massScale + impulseScale == 1` always
- High hertz → massScale→1, impulseScale→0 (rigid)
- Low hertz → massScale→0, impulseScale→1 (soft, previous impulse dominates)

### Key Architectural Features Relevant to Your Problem

**1. Graph coloring for parallelism** (`@/home/prokop/git/box3d/src/constraint_graph.c:216-275`): Constraints touching disjoint bodies are grouped into "colors" and solved in parallel. This is critical for molecular systems with thousands of bonds — bonds that don't share atoms can be solved simultaneously.

**2. Warm starting** (`@/home/prokop/git/box3d/src/distance_joint.c:310-347`): Accumulated impulses from the previous step are re-applied. This is the main reason SI converges fast — for a system near equilibrium, the first iteration is already close to the solution. This is analogous to carrying over Lagrange multipliers in XPBD.

**3. Soft vs rigid constraint unification**: The distance joint (`@/home/prokop/git/box3d/src/distance_joint.c:349-536`) shows both modes in one code path. When `enableSpring=true` with `hertz > 0`, it uses the soft formulation (massScale/impulseScale from [b3MakeSoft](cci:1://file:///home/prokop/git/box3d/src/solver.h:263:0-305:1)). When rigid, it uses bias-based position correction. This is exactly the hard/soft DOF splitting you need.

**4. Inequality constraints (contacts)** (`@/home/prokop/git/box3d/src/contact_solver.c:418-470`): Normal impulses are clamped to `≥ 0` (non-penetration). This is the same structure as your hard-sphere collisions: `λ = max(0, λ_accumulated + Δλ)`.

### SI vs XPBD/PD for Your Molecular Problem

| Aspect | Sequential Impulses (Box3D) | XPBD | Projective Dynamics |
|--------|---------------------------|------|---------------------|
| **Stability** | Excellent — velocity-level, implicit damping via softness params | Good — position-level, but stiffness limited by `dt` | Good for soft, bad for stiff (needs many iterations) |
| **Momentum conservation** | **Exact** per constraint — impulses are equal-and-opposite by construction (J·M⁻¹·Jᵀ formulation) | Exact — constraint forces are symmetric | Exact for linear constraints |
| **Hard/soft splitting** | Natural — just set hertz high for bonds, low for dihedrals | Natural — compliance parameter | Natural — alternating local/global |
| **Warm starting** | Built-in, critical for convergence | Built-in (Lagrange multiplier carry) | Not typical |
| **Parallelism** | Graph coloring — constraints on disjoint bodies solved in parallel | Same possible, but typically sequential | Global solve (Cholesky) or Jacobi |
| **Speed per iteration** | Very fast — scalar effective mass per 1D constraint, no global solve | Similar to SI | Global solve is expensive for large systems |
| **Convergence rate** | Slow for stiff systems without substepping (Gauss-Seidel) | Similar to SI | Better for soft constraints (Newton on global system) |

### My Assessment for Your Use Case

**SI is well-suited for your molecular relaxation problem**, and here's why:

**Advantages:**

- **Your hard DOFs (bonds, angles, hard-sphere collisions) map directly** to Box3D's distance joints, revolute/weld joints, and contact constraints. The impulse-clamping for contacts (`max(0, ...)`) is exactly hard-sphere non-penetration.

- **Your soft DOFs (dihedrals, Coulomb, van der Waals)** map to soft joints with low hertz. The [b3MakeSoft](cci:1://file:///home/prokop/git/box3d/src/solver.h:263:0-305:1) formulation gives you tunable stiffness + damping in one parameter pair (hertz, damping ratio), which is simpler than choosing XPBD compliance values.

- **Momentum conservation is exact** — each impulse update is `vA -= mA·P, vB += mB·P` with the same `P`. Linear and angular momentum are conserved per-constraint. This is a structural property of the J·M⁻¹·Jᵀ formulation, not a numerical accident.

- **Graph coloring gives you free parallelism** — molecular bond graphs are sparse and locally tree-like, so they color well (typically 4-8 colors for 3D molecular graphs). This means 4-8× speedup on multicore CPUs.

- **Warm starting** means that after the first few timesteps, the system is near equilibrium and each step converges in 1-2 iterations. This is the key to fast relaxation — you're not solving from scratch each step.

**The main limitation for your splitting scheme:**

- **Gauss-Seidel convergence for stiff systems is slow.** SI with ITERATIONS=1 per substep is a single Gauss-Seidel sweep. For very stiff bonds (high hertz), you need many substeps. Box3D handles this by substepping (`subStepCount` in `@/home/prokop/git/box3d/src/physics_world.c:1094`). Your "hard DOFs separately, then propagate soft" splitting is not directly supported — SI solves all constraints together in each iteration.

- **No native multi-rate splitting.** Box3D doesn't have a built-in mechanism to solve hard constraints with small dt and soft constraints with large dt. You'd need to implement this yourself: e.g., run N substeps with only bond/angle constraints active (high hertz), then one step with dihedral/Coulomb/VdW active (low hertz).

### Practical Recommendation

If you want to use SI for molecular relaxation, the most pragmatic approach is:

1. **Represent bonds as distance joints** with high hertz (or rigid, hertz=0 with bias)
2. **Represent angles as 3-DOF constraints** (similar to weld/spherical joints with angular limits)
3. **Represent hard-sphere collisions as contacts** (inequality, impulse clamped ≥ 0)
4. **Represent dihedrals as soft angular springs** (low hertz, damping ratio ~0.7)
5. **Represent Coulomb/VdW as external forces** applied during velocity integration (not as constraints — they're not stiff enough to need constraint treatment, and making them forces avoids the Gauss-Seidel convergence penalty)

The substep count becomes your splitting parameter: use more substeps when the system is far from equilibrium (large residual forces), fewer when it's relaxed. Warm starting makes each substep cheap after the first few.

**Bottom line:** SI is likely faster than XPBD or PD for your case because (a) warm starting eliminates most of the solve cost after the initial transient, (b) graph coloring parallelizes well on molecular graphs, and (c) the per-iteration cost is minimal (scalar effective mass, no global linear solve). The tradeoff is that you lose the global convergence guarantee of PD's Cholesky solve — but since you said accuracy doesn't matter and you want fast relaxation to a qualitative minimum, SI's bias-corrected Gauss-Seidel is adequate.

---

# Chat GPT 5.5

## Verdict

**Yes, sequential impulses can be very useful for your molecular relaxation problem, especially for hard-sphere contacts and rigid bond/angle constraints. But SI is not intrinsically more stable or faster than XPBD in a fundamental mathematical sense.**

The main advantage comes from the particular combination used by Box3D:

* projected Gauss–Seidel updates,
* accumulated impulses and warm starting,
* substeps rather than many iterations,
* explicit inequality handling for contacts,
* graph coloring, SIMD, and cache-friendly data layout,
* a second “relaxation” solve that removes artificial constraint energy.

Current Box3D calls this solver **Soft Step**. It is essentially Catto’s TGS-Soft scheme rather than just the old 2006 sequential-impulse algorithm. ([GitHub][1])

For your application, my recommendation is:

> **Use a custom SI/PGS solver as the hard-constraint projector inside a FIRE-like or heavy-ball molecular minimizer. Do not replace the complete molecular relaxation algorithm with an unmodified Box3D simulation loop.**

Your current PD collision plan introduces temporary harmonic constraints for close pairs and then adds a smoothing potential around the cutoff. SI lets you replace that whole repulsive inner region by a genuine inequality constraint. 

---

# 1. SI, XPBD, and PD are solving closely related problems

Suppose the generalized coordinates are (q), generalized velocities are (u), and a hard constraint is

[
C_a(q)=0
]

or, for a collision,

[
C_a(q)\ge 0.
]

Its Jacobian is

[
J_a = \frac{\partial C_a}{\partial q}.
]

At one velocity step, we seek the smallest mass-weighted velocity correction that satisfies the linearized constraint:

[
\min_{\Delta u}
\frac12 \Delta u^\mathrm{T}M\Delta u
]

subject to

[
J_a(u+\Delta u)+b_a=0.
]

The solution has the form

[
\Delta u=M^{-1}J_a^\mathrm{T}\Delta p_a,
]

where the scalar or small-vector impulse is

[
\Delta p_a
==========

-\left(J_aM^{-1}J_a^\mathrm{T}\right)^{-1}
\left(J_au+b_a\right).
]

Define the effective constraint mass

[
K_a=J_aM^{-1}J_a^\mathrm{T}.
]

Sequential impulses simply solve one row or small block at a time:

[
u\leftarrow u+M^{-1}J_a^\mathrm{T}\Delta p_a.
]

The next constraint immediately sees the corrected velocity. That is projected Gauss–Seidel.

For a contact, the accumulated normal impulse is projected:

[
p_a^{\mathrm{new}}
==================

\max\left(0,p_a^{\mathrm{old}}+\Delta p_a\right).
]

This is exactly the right mathematical structure for hard-sphere nonpenetration.

## Comparison with XPBD

XPBD performs an almost analogous update in position space:

[
\Delta\lambda_a
===============

\frac{-C_a-\tilde\alpha_a\lambda_a}
{J_aM^{-1}J_a^\mathrm{T}+\tilde\alpha_a},
\qquad
\tilde\alpha_a=\frac{\alpha_a}{h^2},
]

followed by

[
q\leftarrow q+M^{-1}J_a^\mathrm{T}\Delta\lambda_a.
]

XPBD was specifically introduced to remove the timestep- and iteration-dependent stiffness of ordinary PBD. ([matthias-research.github.io][2])

Therefore:

* SI is a velocity/impulse projection.
* XPBD is a position/multiplier projection.
* With similar Gauss–Seidel ordering, compliance, warm starting, and substepping, they are much more alike than their names suggest.

A well-implemented XPBD solver will not suddenly be much worse than SI. Box3D’s advantage is largely the complete solver design and engineering, not a magical property of impulses.

## Comparison with Projective Dynamics

PD performs local constraint projections but combines them with a global linear solve. For fixed topology, this propagates information globally and can converge much faster than local SI sweeps. Its great advantage is that its global matrix can often be prefactorized. ([ACM Digital Library][3])

But dynamic contacts are exactly where PD becomes awkward:

* contact rows appear and disappear,
* contacts are inequalities,
* the matrix topology changes,
* factorization reuse becomes difficult,
* one needs active-set or complementarity logic.

That is where SI is substantially more natural.

---

# 2. Mapping your molecular degrees of freedom

## Bond lengths

For atoms (i,j),

[
C_{ij}=|r_j-r_i|-l_{ij}=0.
]

With

[
n=\frac{r_j-r_i}{|r_j-r_i|},
]

the Jacobian is

[
J=[-n,,+n],
]

and for point masses

[
K=w_i+w_j,\qquad w_i=\frac1{m_i}.
]

The impulse is equal and opposite:

[
\Delta v_i=-w_i n\Delta p,
\qquad
\Delta v_j=+w_j n\Delta p.
]

This preserves total linear momentum exactly, apart from floating-point roundoff.

For relaxation, I would normally make covalent bond lengths genuinely hard, (\alpha=0), unless their elastic response under load matters for the approximate geometry.

## Bond angles

There are two efficient possibilities.

### 1–3 distance constraint

If (l_{ij}), (l_{jk}), and the target angle (\theta_0) are fixed, constrain

[
l_{ik}^2
========

l_{ij}^2+l_{jk}^2
-2l_{ij}l_{jk}\cos\theta_0.
]

This converts the angle into another ordinary distance constraint. It is very cheap and robust.

This is the best option when bond lengths are already nearly exact.

### Direct three-body angular constraint

Define

[
a=r_i-r_j,\qquad b=r_k-r_j,
]

[
\hat a=\frac{a}{|a|},\qquad
\hat b=\frac{b}{|b|},
]

and avoid the expensive and poorly conditioned (\arccos):

[
C=\hat a\cdot\hat b-\cos\theta_0.
]

Writing (c=\hat a\cdot\hat b), its gradients are

[
g_i=\frac{\hat b-c\hat a}{|a|},
]

[
g_k=\frac{\hat a-c\hat b}{|b|},
]

[
g_j=-g_i-g_k.
]

Then

[
K=w_i|g_i|^2+w_j|g_j|^2+w_k|g_k|^2,
]

and

[
\Delta v_s=w_sg_s\Delta p.
]

Because

[
g_i+g_j+g_k=0,
]

the constraint creates no net linear impulse. Since the constraint is invariant under global rotation, its gradients also have zero net torque analytically:

[
\sum_s r_s\times g_s=0.
]

Thus a properly implemented angle row preserves both linear and angular momentum.

## Hard-sphere collisions

For two atomic cores,

[
C_{ij}=|r_j-r_i|-R_{ij}\ge0.
]

The accumulated impulse is constrained by

[
p_{ij}\ge0.
]

This is better than your temporary harmonic collision pair because:

* it never attracts atoms after they separate;
* no inner harmonic cutoff is necessary;
* no smoothing function is needed at the inner/outer boundary;
* arbitrarily large penalty stiffness is unnecessary;
* existing penetration is removed using the bias or position-relaxation stage.

For a custom molecular kernel you can even avoid a square root during much of the calculation by using

[
C_2=\frac12\left(|d|^2-R^2\right)\ge0,
\qquad d=r_j-r_i,
]

with Jacobian

[
J=[-d,+d].
]

Near (r=R), this has the same local constraint surface. Constraint scaling affects compliance and bias parameters, so those must be calibrated consistently, but for hard relaxation constraints it is a useful optimization.

Set contact friction and rolling resistance to zero. Game-style Coulomb friction between atomic spheres would create artificial sticking and unwanted metastable structures.

## Rigid atomic ports

For your oriented (sp^3)-like atoms, a port constraint is a point-to-point constraint:

[
x_i+R(q_i)p_{ia}
================

x_j+R(q_j)p_{jb}.
]

An impulse at an offset port generates both translation and torque:

[
\Delta v_i=-m_i^{-1}P,
\qquad
\Delta\omega_i=-I_i^{-1}(r_i\times P),
]

[
\Delta v_j=+m_j^{-1}P,
\qquad
\Delta\omega_j=+I_j^{-1}(r_j\times P).
]

This is naturally a (3\times3) block constraint rather than three unrelated scalar constraints. Box3D’s point-anchor/spherical-joint machinery represents this type of constraint; your previous port analysis identified the same mapping. 

For your own implementation, I would solve all ports incident on one atom as a small (6\times6) translation–rotation block when possible. That will converge faster than scalar SI rows and approaches your local Newton/Kabsch rotation solvers.

---

# 3. The solver I would actually implement

I would use two nested levels.

## Outer level: soft-energy minimization

Evaluate only the relatively smooth forces:

[
F_i^{\mathrm{soft}}
===================

F_i^{\mathrm{dihedral}}
+
F_i^{\mathrm{Coulomb}}
+
F_i^{\mathrm{attractive\ vdW}}
+
F_i^{\mathrm{external}}.
]

These forces must be implemented symmetrically so that internal interactions satisfy

[
\sum_i F_i=0,
\qquad
\sum_i r_i\times F_i=0.
]

Use FIRE, damped heavy-ball, or a similar inertial minimizer for this outer step. FIRE was designed specifically for rapid atomic structural relaxation and is often much faster than plain damped MD. ([math.uni-bielefeld.de][4])

## Inner level: hard SI projection

Apply several cheap substeps or sweeps containing only:

* bond constraints,
* hard angle or port constraints,
* hard-sphere contacts.

A conceptual loop is:

```cpp
for (int outer = 0; outer < maxOuter; outer++) {

    update_neighbor_list(x);
    eval_soft_forces(x, fsoft);

    // FIRE or heavy-ball soft kick
    for (int i = 0; i < natom; i++) {
        v[i] += H * invMassOpt[i] * fsoft[i];
    }

    preserve_or_remove_rigid_body_modes(x, v);

    const double h = H / nHardSubsteps;

    for (int isub = 0; isub < nHardSubsteps; isub++) {

        // Reuse previous bond, angle, port, and persistent-contact impulses
        warm_start_constraints(v, lambda, h);

        // One graph-colored Gauss-Seidel sweep
        solve_bond_constraints(v, x, h);
        solve_angle_or_port_constraints(v, omega, x, q, h);
        solve_contact_constraints(v, omega, x, q, h);

        integrate_positions(x, q, v, omega, h);

        // Remove energy inserted by bias/penetration correction
        solve_constraints_without_bias(v, omega, x, q);
    }

    fire_update(v, fsoft, H);

    if (converged_projected_force_and_constraints()) break;
}
```

The current Box3D source similarly uses a primary solve and a relaxation solve around position integration, while Catto’s Soft Step formulation combines substepping, softness, warm starting, and relaxation. ([GitHub][5])

The observation that (N) small substeps with one constraint iteration can work better than one large step with (N) iterations is also supported by the “Small Steps” study. ([Macklin][6])

## Frozen soft forces

If Coulomb and attractive-force evaluation is expensive, freeze (F^{\mathrm{soft}}) across the inner substeps.

That gives exactly the splitting you want:

[
\text{one expensive soft evaluation}
\quad+\quad
N\text{ very cheap hard projections}.
]

This is not accurate molecular dynamics, but it is quite reasonable for structural relaxation. Add:

* a maximum displacement per macrostep;
* force clipping for badly initialized structures;
* adaptive (H);
* rejection or reduction of (H) when the soft energy rises too much.

---

# 4. SI alone does not solve the slow-soft-mode problem

This is the main limitation.

SI removes the instability caused by high-frequency bond and collision modes. It does **not** automatically accelerate:

* collective torsion motion,
* long-chain bending,
* rigid translation of molecular fragments,
* low-frequency rearrangement of a packed molecular cluster.

A local Gauss–Seidel solver propagates information through the constraint graph gradually. Long polymers and extended networks can still converge slowly.

Projective Dynamics can be better here because its global solve couples distant atoms in one iteration. Likewise, FIRE accelerates motion along persistent soft-force directions.

Therefore the best division of labour is:

[
\boxed{
\text{SI for feasibility and hard modes}
+
\text{FIRE or preconditioned descent for soft modes}
}
]

not

[
\text{SI for everything}.
]

For molecules containing recognizable rigid fragments—phenyls, cages, surface-bound groups—it may be even better to make each fragment a rigid body and optimize only fragment translations, rotations, and connecting dihedrals.

---

# 5. Momentum conservation during relaxation

## Internal SI constraints

For each correctly constructed constraint,

[
\Delta P_{\mathrm{total}}
=========================

\sum_i m_i\Delta v_i=0.
]

For rotationally invariant constraints,

[
\Delta L_{\mathrm{total}}
=========================

\sum_i r_i\times m_i\Delta v_i
+
\sum_i I_i\Delta\omega_i
=0.
]

This is a genuine advantage of velocity-impulse formulations: momentum conservation is structurally explicit.

However, an unmodified game engine may deliberately break these invariants through:

* per-body linear damping;
* per-body angular damping;
* velocity caps;
* sleeping;
* motors;
* externally imposed position corrections;
* friction against static geometry.

The uploaded architecture review correctly emphasizes that Box3D is a general rigid-body game engine rather than a specialized particle solver. 

For molecular relaxation:

* disable gravity;
* disable restitution;
* disable contact friction;
* disable sleeping until your own convergence test succeeds;
* avoid independent per-atom damping.

## Momentum-neutral damping

Separate the rigid-body motion of the complete isolated system:

[
V_{\mathrm{cm}}
===============

\frac{\sum_i m_iv_i}{\sum_i m_i},
]

and find (\Omega) from

[
I_{\mathrm{tot}}\Omega
======================

\sum_i (r_i-R_{\mathrm{cm}})
\times m_i(v_i-V_{\mathrm{cm}}).
]

Write

[
v_i
===

V_{\mathrm{cm}}
+
\Omega\times(r_i-R_{\mathrm{cm}})
+
v_i^{\mathrm{int}}.
]

Damp only the internal component:

[
v_i
\leftarrow
V_{\mathrm{cm}}
+
\Omega\times(r_i-R_{\mathrm{cm}})
+
\eta v_i^{\mathrm{int}}.
]

This removes vibrational energy without changing total linear or angular momentum.

For a minimization run, I would normally initialize

[
P_{\mathrm{total}}=0,\qquad L_{\mathrm{total}}=0
]

and preserve these zero modes.

---

# 6. Which method is likely fastest?

| Situation                       | Most promising method                                 |
| ------------------------------- | ----------------------------------------------------- |
| Dynamic hard-sphere contacts    | SI / projected Gauss–Seidel                           |
| Fixed covalent topology on CPU  | Colored SI, block-SI, or prefactored PD               |
| Fixed covalent topology on GPU  | XPBD/PD Jacobi or local block solves may win          |
| Long stiff chains/networks      | PD global solve usually propagates information faster |
| Pure point particles            | Custom SI is much cheaper than full Box3D             |
| Oriented rigid atoms with ports | Block SI is attractive                                |
| Fast local minimum search       | Projected FIRE + SI hard projector                    |
| Accurate finite bond stiffness  | XPBD compliance or explicit/implicit harmonic energy  |

Catto’s own Solver2D comparison does not establish that SI is universally cheaper than XPBD. At the compared settings, Soft Step and XPBD require similar numbers of constraint traversals; Catto also explicitly notes that Solver2D is not an optimized performance benchmark. ([box2d.org][7])

Box3D itself is highly optimized for large piles of general rigid bodies using data-oriented C17 code, multithreading, graph coloring, and SIMD. That is excellent engineering, but it carries machinery you do not need for atomic spheres: general shapes, contact manifolds, islands, events, sleeping, CCD, and full six-dimensional rigid bodies. ([GitHub][1])

So I would use Box3D as:

1. a reference implementation;
2. a prototype for testing molecular constraints;
3. a benchmark against your custom solver;

but probably not as the final high-performance molecular kernel.

---

# 7. Specific optimizations for your implementation

### Precompute colors for permanent constraints

Bond, angle, and port topology rarely changes. Color these constraints once. Within one color, no atom appears in two constraints, so the rows can be evaluated in parallel without atomics.

Dynamic contacts can be placed in temporary colors or solved with a Jacobi accumulation pass.

### Use block constraints

Avoid reducing everything to scalar rows.

* bond: (1\times1);
* angle: (1\times1);
* port coincidence: (3\times3);
* complete oriented atom update: potentially (6\times6).

Small direct block solves can greatly improve convergence for almost no cost.

### Warm-start permanent constraints

Store accumulated multipliers for:

* bonds,
* angles,
* ports,
* persistent contacts.

Warm starting is particularly effective near a converged structure because constraint reactions change only slightly between outer steps.

### Use optimization masses

You are not performing accurate dynamics, so physical atomic masses are not necessarily optimal.

Possible choices are:

[
m_i=1
]

for all atoms, or a diagonal stiffness estimate such as

[
m_i\propto
\sum_{a\ni i} k_a.
]

The latter acts as a simple preconditioner. It prevents highly connected atoms from moving much more than weakly connected atoms.

### Use continuation for difficult packing

Exact hard spheres can jam the structure in a poor local arrangement. A useful relaxation schedule is:

1. begin with slightly reduced hard-core radii;
2. use compliant contacts;
3. gradually increase radii;
4. gradually reduce compliance;
5. finish with exact nonpenetration.

The same can be done with angle and port stiffness.

### Test convergence using projected forces

Raw soft forces need not vanish at a constrained minimum because bond and contact reactions balance their normal components.

The relevant criterion is the force projected onto the feasible tangent space:

[
F_{\mathrm{tan}}
================

P_{\mathrm{constraints}}M^{-1}F.
]

In practice, apply your SI solver to a pseudo-velocity

[
v=M^{-1}F
]

with zero positional bias. The remaining velocity is an approximation to the projected descent direction. Stop when:

[
\max_i |v_i^{\mathrm{projected}}|<\epsilon_F,
]

together with small bond errors, angle errors, and penetration.

---

## Final recommendation

For your solver I would implement:

[
\boxed{
\begin{aligned}
&\text{soft compact-support attraction, Coulomb, dihedrals}\
&\qquad\downarrow\quad\text{FIRE/heavy-ball macrostep}\
&\text{colored block sequential impulses}\
&\quad\text{for bonds, angles/ports, and hard contacts}\
&\qquad\downarrow\
&\text{bias-free relaxation sweep}\
&\qquad\downarrow\
&\text{momentum-neutral damping}
\end{aligned}
}
]

This should be more natural and probably faster than your current iterative-PD treatment of **dynamic collisions**. It will not necessarily beat PD for the permanent covalent network, and it will not by itself accelerate the soft collective modes. The strongest design is therefore a hybrid: **SI as a hard-manifold/contact projector, FIRE or preconditioned descent as the actual minimizer.**

[1]: https://github.com/erincatto/box3d "GitHub - erincatto/box3d: Box3D is a 3D physics engine for games · GitHub"
[2]: https://matthias-research.github.io/pages/publications/XPBD.pdf?utm_source=chatgpt.com "XPBD: Position-Based Simulation of Compliant ..."
[3]: https://dl.acm.org/doi/10.1145/2601097.2601116?utm_source=chatgpt.com "Projective dynamics: fusing constraint projections for fast ..."
[4]: https://www.math.uni-bielefeld.de/~gaehler/papers/fire.pdf?utm_source=chatgpt.com "Structural Relaxation Made Simple"
[5]: https://raw.githubusercontent.com/erincatto/box3d/main/src/solver.c "raw.githubusercontent.com"
[6]: https://mmacklin.com/smallsteps.pdf?utm_source=chatgpt.com "Small Steps in Physics Simulation"
[7]: https://box2d.org/posts/2024/02/solver2d/?utm_source=chatgpt.com "Solver2D"

---

# USER

you say: "The solver I would actually implement, Outer level: soft-energy minimization, Inner level: hard SI projection"

yes, exactly that is what I menat by my analogy to projective dynamics or XPBD, the ideas of time splintin outer interial loop and inner relaxation loop with interial diagonal penalty M/dt^2 to allow seamless connection between outer and inner loop (initerial is explicit in outer loop, but only as regularization term in inner relaxation loop).

you say "substeps rather than many iterations", this puzzles me. For me sub-step and iteration is the same thing. I have outer loop (where I compute costly long-range non-covalent interaction O(n^2) but K_soft<<K_hard) and the inner lop where I compute cheap hard bonding and contact collision interation (O(n) but K_hard>>K_soft), that is the key disticition. But how this change sub-step vs iteration?


you say "SI will not automatically solve the problem of slow soft modes". But the think is that the more stable I can get the hard DOFs the longer time step I can allow for soft DOFs. My goal is to allow as long time step as possible (like SHAKE iused in MD from 1977 but without the cost of rigorous matrix solver for lagrange multiplies).

you say "PD may be better because it explicitly solve the hard constrain - e.g. using Cholesky/LU" but I do not wan to do costly linear solver. In fact I kinda like about SI that it works woth impulses which are more close to Forces used in normal molecualr dymanamics. 

As I said I do not relly care about accuracy here - my "Hard sheres" does not need to be really hard, in reality atoms are soft, hard sphere is just apporximation of the reall repulsive potential. I do not really care about exact shape of the potential, the only think I care about that the sphere sill avoid penetrating each other too much. The same with angles and both lenghs - harmonic constrains are prefectly fine, I I do not care abut ensuring exacly distances or angles. I just need fast solver relaxing to approximative distances and angles, which allow me to do as long time step as possible on the soft DOFs. 

I'm not sure if the warm starint is usefull as it is aditional degree for feedom which I have to read from memory (and on GPU I'm memory bound). 

I don;t know how to emphesise enoguhg: Yes, I solve similar problem as game engine - high stiffnes ratio. But the goal is not to get as accurate as possible contrain satisfaction at decent time step. But allow as long as possible time step with decent (approximative) constrain satistication. It is not about accuracy but about speed.  

---


Yes—your formulation is the right one, and my previous answer overemphasized **constraint accuracy**. Your actual objective is:

[
\boxed{\text{maximize the outer soft-force timestep }H}
]

while spending only a few cheap (O(N)) passes preventing bonds, angles, and atomic cores from becoming unreasonably distorted.

You are not trying to solve a rigorous constrained MD problem. Therefore exact Lagrange multipliers, exact complementarity, and highly converged inner solves are unnecessary.

## 1. Substep versus iteration: the precise difference

Suppose you evaluate the expensive soft force once per macrostep:

[
F_s^n=F_s(x^n),
]

and form the inertial predictor

[
y=x^n+Hv^n+H^2M^{-1}F_s^n.
]

Your PD-like inner problem is

[
x^{n+1}\approx
\arg\min_x
\left[
\frac{1}{2H^2}(x-y)^TM(x-y)
+
E_{\rm hard}(x)
\right].
]

### Several inner iterations

You keep the same:

* predictor (y),
* macrostep (H),
* inertial penalty (M/H^2),

and repeatedly improve the solution:

[
x^{(k+1)}=\mathcal R_H(x^{(k)};y).
]

No simulated time passes between iterations. You are solving one proximal minimization problem more accurately.

### Several substeps

Instead split

[
H=s h,\qquad h=\frac{H}{s}.
]

Then perform:

[
v\leftarrow v+hM^{-1}F_s^n,
]

[
v\leftarrow\operatorname{HardSolve}(x,v,h),
]

[
x\leftarrow x+hv,
]

and repeat (s) times, while still using the same frozen expensive (F_s^n).

Each substep:

* advances position,
* creates a new inertial predictor,
* recomputes bond directions, angle Jacobians and contact normals,
* activates or deactivates contacts,
* uses the microstep penalty (M/h^2), not (M/H^2).

Therefore it is a sequence of (s) different approximate implicit problems rather than (s) iterations on one problem.

This is what the “small steps” literature and Catto’s terminology mean: one large step with (s) solver iterations is compared with (s) smaller integrations containing one solver iteration each. The latter was found to behave better for strongly nonlinear constrained systems at similar constraint-traversal cost. ([Macklin][1])

### In your terminology

If your inner loop already does:

```cpp
solve_hard();
integrate_positions();
update_geometry();
```

on every pass, **your inner loop is already substepping**.

If instead it does:

```cpp
for (iter) relax_same_predicted_configuration();
integrate_positions_once_at_end();
```

then it is an iterative solve.

So you were not conceptually missing anything. The distinction is only about whether you integrate and relinearize between hard passes.

---

# 2. Why substeps help you increase the outer timestep

For a completely linear problem with fixed Jacobians, iterations and substeps can be rather similar. Your hard molecular interactions are nonlinear:

[
C_{\rm bond}=|r_{ij}|-l_0,
]

[
C_{\rm angle}
=============

\hat r_{ij}\cdot\hat r_{kj}-\cos\theta_0,
]

[
C_{\rm contact}=|r_{ij}|-R_{ij}.
]

With a very large predictor displacement, one macrostep can produce:

* strongly rotated bond directions;
* atoms jumping through each other;
* a collision normal evaluated at an irrelevant configuration;
* several new contacts appearing simultaneously;
* a poor linearization of angular constraints.

Doing several hard microsteps under the same frozen soft force means the hard manifold is followed gradually:

```text
evaluate expensive Fsoft once
        ↓
soft kick h
hard correction
move h
        ↓
update geometry
soft kick h
hard correction
move h
        ↓
...
```

The total soft impulse remains

[
\sum_{a=1}^{s}hF_s=HF_s,
]

so you have not increased the number of expensive soft-force evaluations. You have only spent more cheap (O(N)) hard passes.

That is exactly aligned with your objective.

Catto additionally uses smaller solver substeps to make soft constraints behave more rigidly and re-evaluates their positional bias during the substeps. ([Box2D][2]) Box3D exposes this as a substepping Soft Step solver. ([Box2D][3])

---

# 3. You do not need hard constraints at all

Given your stated priorities, I would not formulate bonds and contacts as exact constraints. I would use a **local implicit harmonic impulse**.

This has the force-like character you like, but removes the explicit spring stability limit.

Consider any locally scalar hard coordinate

[
C(x),
]

with generalized relative velocity

[
u=Jv,
\qquad
J=\frac{\partial C}{\partial x},
]

and inverse effective mass

[
w=JM^{-1}J^T.
]

For a bond between two point particles,

[
w=\frac1{m_i}+\frac1{m_j}.
]

Take a harmonic spring and damping force:

[
f_c=-kC-cu.
]

Instead of evaluating it explicitly at the beginning of the microstep, evaluate it implicitly at the end:

[
p=hf_c^{,\prime},
]

[
u'=u+wp,
]

[
C'\approx C+hu'.
]

Thus

[
p=-hkC'-hcu'
]

and consequently

[
\boxed{
p=
-\frac{
hkC+h(c+hk)u
}{
1+h(c+hk)w
}
}
]

This is one local implicit spring solve. It requires:

* one constraint value (C);
* one relative velocity (u);
* one effective inverse mass (w);
* one scalar impulse (p).

No global matrix, no Lagrange-multiplier solve and no convergence to an exact constraint.

## Dimensionless form

Define

[
a=h^2kw,
\qquad
d=hcw.
]

Then

[
\boxed{
p=
-\frac1w
\frac{
aC/h+(a+d)u
}{
1+a+d
}
}
]

and the post-correction constraint velocity is

[
u'
==

\frac{u-aC/h}{1+a+d}.
]

For an isolated undamped constraint,

[
C'
==

\frac{C+hu}{1+a}.
]

Therefore (a) has a very intuitive algorithmic meaning:

* (a=1): remove about half of the predicted error;
* (a=3): remove about 75%;
* (a=9): remove about 90%;
* (a\rightarrow\infty): hard projection.

You could avoid thinking in physical stiffness entirely and specify a desired correction fraction

[
\rho=\frac{1}{1+a}.
]

For example,

[
\rho=0.1
\quad\Longrightarrow\quad
a=9.
]

That seems much closer to what you actually want than specifying chemically realistic (k).

---

# 4. Application to your interactions

## Bonds

For

[
C=|r_j-r_i|-l_0,
\qquad
n=\frac{r_j-r_i}{|r_j-r_i|},
]

use

[
u=n\cdot(v_j-v_i),
]

[
w=w_i+w_j.
]

Compute (p) and apply equal and opposite impulses:

[
v_i\leftarrow v_i-w_i pn,
]

[
v_j\leftarrow v_j+w_j pn.
]

This exactly preserves total linear momentum.

You can choose (a_{\rm bond}\sim 3\ldots20) rather than a physical bond stiffness. Coupled bonds will correct less than the isolated estimate, but several microsteps propagate the corrections.

## Angles

Use either a direct three-body angular Jacobian or your rigid-port representation.

For an angle (i-j-k), an inexpensive constraint is

[
C=
\hat r_{ji}\cdot\hat r_{jk}
-\cos\theta_0.
]

Compute its three gradients (g_i,g_j,g_k), then

[
u=g_i\cdot v_i+g_j\cdot v_j+g_k\cdot v_k,
]

[
w=w_i|g_i|^2+w_j|g_j|^2+w_k|g_k|^2.
]

The same scalar implicit impulse formula applies:

[
v_s\leftarrow v_s+w_sg_sp.
]

Because the gradients come from a translation- and rotation-invariant scalar, the resulting internal impulse conserves both linear and angular momentum analytically.

## Soft hard-sphere collision

Define

[
C=r-R.
]

Negative (C) is penetration. Use the same implicit spring impulse but project it to be repulsive:

[
\boxed{
p\leftarrow\max(0,p)
}
]

and apply it along the contact normal.

You may also activate the contact speculatively when

[
C+hu<0,
]

meaning that the pair is predicted to penetrate during the microstep.

This produces exactly the behaviour you describe:

* some penetration is acceptable;
* large penetration is quickly suppressed;
* there is no singular Lennard-Jones wall;
* there is no need for exact complementarity;
* the interaction is unconditionally stable for an isolated linearized contact;
* no attraction is generated by the collision solver.

Your existing PD collision design converts close pairs into temporary harmonic constraints and adds a smooth outer potential.  The local implicit impulse lets you retain the same approximate harmonic-core idea without inserting the pair into a global PD matrix.

---

# 5. Substeps and iterations can still be combined

The most useful hierarchy for you is probably:

```text
soft macrostep H:
    evaluate Fsoft once

    for s hard microsteps, h = H/s:
        apply frozen Fsoft for h

        do one hard SI sweep:
            bonds
            angles / rigid ports
            collisions

        integrate positions by h
```

Possibly use two hard sweeps if the topology is highly coupled:

```text
for microstep:
    soft kick
    hard sweep 1
    hard sweep 2
    integrate
```

But I would first test:

[
\boxed{\text{many microsteps × one hard sweep}}
]

rather than

[
\boxed{\text{one macrostep × many hard sweeps}}.
]

Both cost approximately the same number of constraint evaluations. The first continually refreshes the nonlinear geometry.

For small independent molecules, an especially efficient GPU implementation is:

* one workgroup per molecule or molecular cluster;
* load atom positions and velocities into local memory;
* loop over several microsteps inside one kernel;
* process bonds, angles and internal collisions sequentially or by a few precomputed colors;
* write positions and velocities back only after all microsteps.

Then “sequential” impulses do not imply repeated global-memory traffic or repeated kernel launches.

---

# 6. Warm starting is completely optional for your objective

There are two distinct ideas that are often conflated.

## Accumulated impulse within one solve

For a rigorous unilateral contact solver, one stores

[
\lambda\leftarrow\max(0,\lambda+\Delta\lambda)
]

during multiple PGS iterations. This prevents later iterations from violating the one-sided nature of the contact.

## Persistent warm starting

At the next timestep, initialize (\lambda) from the value stored in the previous timestep.

That second operation is the extra global-memory read you are worried about.

You do **not** need persistent warm starting.

For your first implementation I would use:

* no persistent (\lambda) for bonds;
* no persistent (\lambda) for angles;
* no persistent (\lambda) for contacts;
* one repulsive implicit contact impulse per microstep;
* possibly one additional sweep.

For one contact visit per microstep, simply compute

[
p=\max(0,p_{\rm implicit})
]

from the current (C) and (u). There is nothing to warm-start.

Warm starting is valuable in Box3D because it must keep stacks of resting rigid bodies quiet with very few iterations. Catto’s soft-constraint formula explicitly includes the accumulated impulse. ([Box2D][2]) Your molecules have permanent local constraints, approximate satisfaction is acceptable, and repeated microsteps already provide temporal propagation.

On a bandwidth-bound GPU, your instinct is reasonable: do not add one float per constraint unless profiling shows that it eliminates enough complete constraint passes to pay for itself.

---

# 7. What ultimately limits the macrostep (H)

Once the hard modes are handled implicitly, the maximum useful (H) is no longer determined by

[
\omega_{\rm hard}\sim\sqrt{k_{\rm hard}/m}.
]

It is determined by three remaining things.

## Remaining soft curvature

For explicit soft forces, the largest remaining soft eigenfrequency still gives a scale

[
H\lesssim O\left(\frac1{\omega_{\rm soft}}\right).
]

That is exactly the benefit of eliminating hard modes: the limiting frequency becomes the highest **unprojected soft** frequency.

## Staleness of the frozen soft force

Even if the integrator does not numerically explode, a huge (H) evaluates

[
F_s(x^n)
]

and then moves to a position where that force is no longer a useful approximation.

Hard substeps cannot repair a stale Coulomb, vdW or torsional force direction.

For minimization, the practical control is therefore usually a displacement/trust-region limit:

[
\max_i |H v_i| < \Delta x_{\max}
]

or

[
\max_i
\left|
H^2M_i^{-1}F_i
\right|
<
\Delta x_{\max}.
]

## Contact topology

If a soft kick moves an atom through several radii before collision detection is refreshed, no local contact solver can reconstruct the missed collision reliably.

Your hard microstep size should therefore approximately satisfy

[
h,|v_{ij}| \lesssim \eta R_{ij},
]

with perhaps

[
\eta\sim0.2\ldots0.5.
]

This restriction affects only the cheap hard microstep (h), not the expensive soft-force macrostep (H). Increase the number of hard microsteps when velocities become large.

---

# 8. A concrete minimal algorithm

```cpp
// Expensive macrostep
evalSoftForces(pos, fsoft);          // Coulomb, attraction, torsion, etc.

const float h = H / nsub;

for (int isub = 0; isub < nsub; isub++) {

    // Frozen soft force, distributed over microsteps
    for (int i = 0; i < natom; i++) {
        vel[i] += h * invMass[i] * fsoft[i];
    }

    // One or optionally two colored hard sweeps
    for (int color = 0; color < nBondColors; color++) {
        solveImplicitBonds(color, h, aBond, dBond);
    }

    for (int color = 0; color < nAngleColors; color++) {
        solveImplicitAnglesOrPorts(color, h, aAngle, dAngle);
    }

    updateLocalCollisionPairs(pos);

    solveImplicitRepulsiveContacts(
        h,
        aCollision,
        dCollision
    );

    // Actual substep: this is what distinguishes it from an iteration
    for (int i = 0; i < natom; i++) {
        pos[i] += h * vel[i];
    }
}

// Strong relaxation damping or FIRE update
updateOptimizerVelocity(pos, vel, fsoft);
```

Reasonable starting values might be:

[
a_{\rm bond}=10,
\qquad
a_{\rm angle}=3\ldots10,
\qquad
a_{\rm collision}=5\ldots20,
]

with enough damping to prevent the constraint impulses from producing oscillations. These are algorithmic parameters, not physical force constants.

Then tune mainly:

* (H): outer motion per expensive force evaluation;
* (n_{\rm sub}): number of cheap hard corrections;
* (a): fraction of hard error removed per visit;
* optimizer damping/FIRE parameters.

## The corrected central recommendation

Your target is not:

> solve hard constraints accurately at a conventional timestep.

It is:

> make hard interactions sufficiently implicit and dissipative that they cease to determine the timestep, then push the expensive soft-force macrostep as far as force staleness and basin geometry permit.

For this objective, **one local implicit SI sweep per hard microstep, no warm starting, no exact constraints, and no global linear solve** is a very plausible optimum. It is essentially a cheap approximate proximal solver where the constraint reaction is represented as an impulse rather than as a converged position projection.

[1]: https://mmacklin.com/smallsteps.pdf "Small Steps in Physics Simulation"
[2]: https://box2d.org/posts/2024/02/solver2d/ "
        
            Solver2D :: Box2D 
        
    "
[3]: https://box2d.org/posts/2026/06/announcing-box3d/ "
        
            Announcing Box3D :: Box2D 
        
    "

