
# Box3D Architecture Review

## Introduction

Box3D is a **3D rigid body physics engine** written in portable C17. It is designed for real-time games and simulations, emphasizing determinism, cache efficiency, and multi-core CPU parallelism. The codebase in `src/` is approximately 80 files organized around a **data-oriented design** with opaque IDs, contiguous arrays, and a multithreaded solver.

CODEMAP:
Box3D Physics Engine: Core Simulation Pipeline
https://windsurf.com/codemaps/6381f6db-3340-4f66-af0f-3bc48bb53fed-fe86ab10a43f3d18

### What Box3D Is

Box3D simulates the motion of **rigid bodies** — solid objects that do not deform under forces. Each body has a position, orientation (quaternion), linear velocity, and angular velocity. Bodies interact through **contacts** (collision) and **joints** (constraints). The engine advances time in discrete steps, computing new velocities and positions while resolving all constraints so that bodies don't interpenetrate and joints stay connected.

### What Box3D Is Not

Box3D is **not** a soft body, cloth, fluid, or particle engine. It does not simulate deformable objects, mass-spring systems, finite element methods, or SPH fluids. These are fundamentally different problem domains requiring different solvers and data structures.

Box3D is also **not GPU-accelerated**. All physics computation runs on the CPU. The engine uses two CPU parallelism strategies:

- **SIMD (Single Instruction, Multiple Data)** — SSE2 on x86, NEON on ARM. Processes 4 contact constraints simultaneously by packing them into "wide" data structures. This is intra-core parallelism, exploiting data-level parallelism within a single CPU core.
- **Multithreading** — Multiple CPU cores work on the same time step in parallel. The solver partitions constraints into groups that touch disjoint sets of bodies (graph coloring), so different cores can solve different groups without data races. This is inter-core parallelism.

The only GPU usage in the repository is in the **samples application**, which uses sokol_gfx (OpenGL) for *rendering* — drawing the physics world to the screen. No physics computation is offloaded to the GPU.

This is a deliberate design choice. CPU physics engines offer better determinism (identical results across runs and platforms), easier debugging, predictable latency, and simpler integration with game logic. GPU physics is harder to make deterministic and introduces frame-latency pipeline stalls.

---

## What Can Be Simulated

### Rigid Body Dynamics

Full 6-DOF (degrees of freedom) rigid body simulation: 3 translational + 3 rotational. Three body types:

- **Static** — never moves, infinite mass. Used for walls, floors, terrain. No state is stored for static bodies in the solver, saving memory and cache pressure.
- **Kinematic** — moves according to user-specified velocity, but is not affected by forces or collisions. Used for moving platforms, doors, elevators.
- **Dynamic** — full simulation. Affected by gravity, forces, torques, impulses, and collisions. Has finite mass and inertia.

Forces, torques, and impulses can be applied to dynamic bodies at any point in world space. Gravity is configurable per-world. Linear and angular damping prevent perpetual motion.

### Contacts and Friction

When two shapes overlap, the engine creates a **contact** — a persistent constraint that prevents further interpenetration. Each contact has:

- **Normal impulse** — prevents shapes from pushing through each other along the contact normal direction. This is the "non-penetration" constraint.
- **Tangent (friction) impulse** — resists relative sliding motion along the contact surface. Uses the **Coulomb friction model**: the friction force is bounded by `μ × N`, where `μ` is the friction coefficient and `N` is the normal force. This means friction is proportional to how hard the surfaces are pressed together.
- **Rolling resistance** — for spheres and capsules, an additional torque resists rolling. This simulates the energy loss when a ball rolls across a rough surface.
- **Restitution (bounce)** — when shapes collide at significant speed, a portion of the approaching velocity is reflected. The coefficient of restitution `e` ranges from 0 (no bounce) to 1 (perfect bounce). To prevent jitter at low velocities, restitution is only applied when the approach velocity exceeds a configurable threshold.
- **Tangent velocity** — per-material conveyor belt support. A surface can be given a tangential velocity, causing objects resting on it to be carried along, simulating a conveyor belt.

### Collision Shapes

Six shape types, each with different collision algorithms:

- **Sphere** — simplest shape. Fast collision tests. Has rolling resistance.
- **Capsule** — a cylinder with hemispherical caps. Good for characters and rounded objects. Has rolling resistance.
- **Convex hull** — arbitrary convex polyhedron defined by vertices and faces. The most general convex shape. Collision detection uses **SAT** (Separating Axis Theorem) for hull-hull and **GJK** (Gilbert-Johnson-Keerthi) for distance queries. Identical hulls are deduplicated by content hash and shared via reference counting.
- **Triangle mesh** — arbitrary triangle soup for static geometry. Has an internal BVH (Bounding Volume Hierarchy) for fast queries. Supports per-triangle materials. Edge identification helps suppress "ghost collisions" — false contacts at shared mesh edges.
- **Height field** — regular grid of heights for terrain. Supports holes, compressed heights, and per-cell materials. Static only.
- **Compound** — aggregates multiple sub-shapes (capsules, hulls, meshes, spheres) into one baked collision shape. Designed for offline baking and world streaming. Materials are deduplicated.

### Terrain

Terrain is supported via **height fields** (regular grids) and **triangle meshes** (arbitrary geometry). Both are static-only shapes. Any static body with any shape acts as immovable ground.

### Joints

Nine joint types connect bodies with constraints that remove degrees of freedom:

| Joint | Description |
|---|---|
| **Distance** | Maintains a target distance between two anchor points. Supports springs (stiffness in Hz + damping ratio), limits (min/max length), and a motor. |
| **Revolute** | Hinge joint — allows rotation about one axis only. Supports angle limits, motor, and spring. Useful for doors, wheels on an axle, pendulums. |
| **Prismatic** | Slider joint — allows translation along one axis only. Supports limits, motor, and spring. Useful for pistons, sliding doors. |
| **Spherical** | Ball-and-socket — allows 3-DOF rotation, no translation. Supports cone limits and motor. Useful for ragdolls, chains. |
| **Weld** | Rigidly connects two bodies. In practice the solver is approximate, so weld joints can be slightly soft — chains of welds may flex. Supports spring parameters for intentional softness. |
| **Wheel** | Combines a suspension spring + steering axis + rotation motor. Designed for vehicle wheels. |
| **Motor** | Applies relative motion (linear or angular) between two bodies without restricting any DOFs. Useful for controlled movement without hard constraints. |
| **Parallel** | Keeps body axes parallel, like a universal (Cardan) joint. Allows rotation but constrains orientation. |
| **Filter** | No constraint at all — only disables collision between the two connected bodies. Useful for bodies connected by other joints that shouldn't also collide with each other. |

All joints support **springs** (specified as frequency in Hz and damping ratio, not raw spring constants), **limits**, and **motors**.

### Continuous Collision Detection (CCD)

Fast-moving objects can "tunnel" through thin obstacles in discrete time stepping — a bullet might pass completely through a wall in a single step. Box3D provides two mechanisms:

- **Speculative contacts** — when AABBs (bounding boxes) overlap but shapes haven't touched yet, a contact constraint is created preemptively. The solver resolves it as if the shapes are already touching, preventing tunneling for moderately fast objects.
- **Time of Impact (TOI)** — for very fast objects ("bullets"), the engine sweeps the shape from its previous position to its current position and computes the exact time fraction at which it first touches another shape. The body is then moved only to that point, preventing tunneling entirely.

### Sensors

Shapes can be marked as **sensors** — they detect overlap with other shapes but generate no contact constraints. Sensors maintain double-buffered visitor lists to detect enter/exit events. Useful for triggers, area detection, and game logic.

### Events

Box3D generates several event types after each step:

- **Body move events** — bodies that changed position (efficient alternative to polling all bodies)
- **Contact begin/end events** — when two shapes start or stop touching
- **Contact hit events** — contacts with significant approach velocity (for sound/particle effects)
- **Sensor begin/end events** — when a shape enters or leaves a sensor zone
- **Joint events** — joints that exceeded force/torque thresholds (for breakable joints)

---

## Solver Theory

This section explains the mathematical methods Box3D uses to simulate physics. Understanding these is essential for using the engine effectively — choosing sub-step counts, understanding why joints stretch, and knowing why the simulation behaves the way it does.

### The Constraint Problem

A physics engine's core job is to advance bodies forward in time while satisfying constraints. The two main constraint types are:

1. **Non-penetration constraints** — shapes must not overlap. When they do, a contact constraint pushes them apart.
2. **Joint constraints** — connected bodies must maintain specified relationships (fixed distance, shared axis, etc.).

Each constraint can be expressed as a **velocity constraint**: the relative velocity between two bodies at the constraint point must satisfy some condition. For a contact, the relative velocity along the normal must be non-negative (bodies must not be approaching). For a distance joint, the relative velocity along the connection axis must be zero.

Mathematically, for each constraint `i` we have:

```
J_i · v = 0     (equality constraint, e.g. distance joint)
J_i · v ≥ 0     (inequality constraint, e.g. non-penetration)
```

where `J_i` is the **Jacobian** (a row vector that extracts the relevant component of relative velocity) and `v` is the vector of all body velocities. The solver must find impulses `λ` such that applying them satisfies all constraints simultaneously.

This is a **Linear Complementarity Problem (LCP)** — a well-known problem in optimization. The inequality constraints (contacts) make it harder than a simple linear system.

### Sequential Impulses (SI)

Box3D uses the **Sequential Impulses** method, developed by Erin Catto (the same author behind Box2D). SI is an iterative method for solving the velocity-level LCP.

**How it works:**

1. **Integrate velocities** — apply external forces (gravity) to get preliminary velocities: `v = v₀ + h · f / m`
2. **Iterate** over all constraints N times (typically 1 iteration per sub-step):
   - For each constraint, compute the relative velocity at the constraint point
   - Compute the impulse needed to satisfy the constraint: `λ = -(J·v) / (J·M⁻¹·Jᵀ)`
   - Apply the impulse to both bodies: `v_A -= λ · M_A⁻¹ · J_A`, `v_B += λ · M_B⁻¹ · J_B`
   - Clamp `λ` for inequality constraints (contacts: accumulated normal impulse must be ≥ 0; friction: accumulated tangent impulse must be ≤ μ × normal impulse)
3. **Integrate positions** — advance positions using the corrected velocities: `x = x₀ + h · v`

The key idea: instead of solving the giant coupled system all at once (which would require matrix inversion), we solve **one constraint at a time**, updating velocities immediately. Each constraint's solution affects the velocities seen by the next constraint. By iterating multiple times, the constraints converge toward a mutually consistent solution.

**Relation to Gauss-Seidel:** SI is mathematically equivalent to **Projected Gauss-Seidel** (PGS) applied to the velocity LCP. "Gauss-Seidel" means we use the most recently updated values immediately (unlike Jacobi which waits until all constraints are processed before updating). "Projected" means we clamp the solution to satisfy inequality constraints (the complementarity conditions). The term "Sequential Impulses" is Catto's more intuitive framing: think of applying corrective impulses one at a time.

**Why not Jacobi?** In Jacobi iteration, all constraints are solved using velocities from the *previous* iteration, then all updates are applied at once. This converges slower than Gauss-Seidel because information propagates more slowly between constraints. However, Jacobi is trivially parallelizable (all constraints are independent within an iteration). Gauss-Seidel is inherently sequential — each constraint depends on the result of the previous one.

Box3D gets the best of both worlds through **graph coloring** (explained below): constraints are grouped into "colors" where no two constraints in the same color share a body. Within a color, constraints can be solved in any order (like Jacobi) because they touch independent bodies. Colors are solved sequentially (like Gauss-Seidel), so information propagates between colors. This allows multi-core parallelism while retaining Gauss-Seidel-like convergence.

### Warm Starting

**The problem:** At each time step, the solver starts from scratch and must converge to a solution. For resting contacts (a stack of boxes on the floor), the correct impulses are nearly the same every step. Re-solving from zero each time is wasteful and can cause visible jitter as the solver "rediscovers" the solution.

**The solution:** Warm starting. At the end of each step, the accumulated impulse for each contact is saved. At the start of the next step, this saved impulse is applied immediately, before the main solve iterations begin. This gives the solver a head start — it begins from a configuration that is already nearly correct.

**Why it matters:** Without warm starting, a stack of 10 boxes would need many iterations to settle each frame, and would likely jitter visibly. With warm starting, the stack remains stable with just 1-2 iterations per sub-step because the impulses from the previous frame are already approximately correct.

Warm starting is enabled by default in Box3D (`b3WorldDef::enableWarmStarting = true`).

### Sub-Stepping

**The problem:** A single solver iteration can only propagate constraint information across one constraint per iteration. If you have a chain of 10 links, information needs to travel from one end to the other — this takes at least 10 iterations. With few iterations, the chain stretches under gravity.

**The solution:** Sub-stepping. Instead of taking one large time step `dt`, the engine takes `N` smaller steps of size `dt/N`. Each sub-step runs the full solver pipeline (integrate, solve, integrate positions). This gives the solver more opportunities to propagate information.

```
Full step:     dt = 1/60s, 1 solve iteration → information travels 1 constraint
Sub-stepped:   4 × (1/240s), 1 solve iteration each → information travels 4 constraints
```

The recommended sub-step count is **4**. More sub-steps improve joint stiffness and reduce stretching, at the cost of more computation. Long joint chains benefit from higher sub-step counts.

### Position Correction (Relax / Stabilization)

**The problem:** The velocity solver prevents *future* penetration, but it doesn't fix *existing* penetration. If two shapes are already overlapping, the velocity constraint says "stop approaching" but doesn't push them apart. Over many frames, small numerical errors accumulate and objects slowly sink into each other.

**The solution:** A separate **position correction** stage (called "relax" in Box3D). After the velocity solve, the engine applies **bias impulses** that push overlapping shapes apart. The bias is proportional to the penetration depth:

```
λ_bias = -β · (penetration / h)
```

where `β` is a bias factor (0 = no correction, 1 = full correction in one step) and `h` is the sub-step time. This is **not** Position-Based Dynamics (PBD) — it works in velocity space, adding extra velocity to separate the shapes, then integrating that velocity into position.

The bias factor is derived from the **softness model** (explained next), which controls how quickly position errors are corrected.

### Soft Constraints (Spring-Damper Model)

Real-world contacts and joints are not perfectly rigid. A rubber ball deforms slightly. A hinge has some play. Box3D models this using **soft constraints** — instead of infinitely stiff constraints, each constraint has a configurable stiffness (frequency in Hz) and damping ratio.

The softness parameters are computed by `b3MakeSoft(hertz, zeta, h)`:

```
ω = 2π · f          (angular frequency from Hz)
a₁ = 2ζ + hω        (damping + stiffness term)
a₂ = hω · a₁        (combined term)
a₃ = 1 / (1 + a₂)

biasRate    = ω / a₁       (how fast position error is corrected)
massScale   = a₂ · a₃      (effective mass scaling)
impulseScale = a₃          (how much previous impulse is retained)
```

- **High frequency (stiff)** → fast correction, but can cause instability if too high relative to the time step. The constraint behaves like a stiff spring.
- **Low frequency (soft)** → slow, gentle correction. The constraint behaves like a soft spring.
- **Damping ratio ζ = 0** → no damping, the spring oscillates. ζ = 1 → critical damping, no oscillation. ζ between 0 and 1 → underdamped (some bounce).

This is the same formulation used by Box2D, derived from Erin Catto's research on soft constraints in sequential solvers. It allows game designers to tune the "feel" of contacts and joints without dealing with raw spring constants.

### What the Solver Is NOT

To avoid confusion, here are methods Box3D does **not** use:

- **Not Position-Based Dynamics (PBD)** — PBD solves constraints directly in position space (move points to satisfy constraints). Box3D works in velocity space: compute corrective impulses → update velocities → integrate positions. PBD is popular for cloth and soft body simulation (e.g., Unity's cloth, NVIDIA Flex). Box3D's approach is more stable for rigid body stacks but cannot simulate deformable objects.

- **Not Extended PBD (XPBD)** — XPBD adds compliance (inverse stiffness) to PBD to make it time-step independent. While conceptually similar to Box3D's soft constraint model, the algorithm is fundamentally different (position projection vs. velocity impulses).

- **Not Conjugate Gradient (CG)** — CG is a Krylov subspace method that solves sparse linear systems exactly in O(n) iterations for chain-like structures. It requires assembling a global system matrix. Box3D never assembles a matrix — it works constraint-by-constraint. CG would be more accurate for large systems but is harder to parallelize and less cache-friendly.

- **Not Jacobi** — Pure Jacobi solves all constraints independently using previous-iteration velocities, then averages. Box3D uses Gauss-Seidel (immediate updates) within graph colors. Pure Jacobi converges roughly twice as slowly as Gauss-Seidel for the same iteration count.

- **Not a direct solver** — Direct solvers (LU decomposition, etc.) factorize the constraint matrix and solve exactly in one pass. They are O(n³) in the number of constraints, impractical for real-time physics with thousands of contacts.

### Pros and Cons of the Sequential Impulse Approach

**Pros:**
- **Fast** — O(n × iterations) where n is the number of constraints. No matrix assembly or factorization.
- **Cache-friendly** — constraints are processed one at a time, accessing only the two bodies involved. Fits well with data-oriented design.
- **Parallelizable** — graph coloring allows multi-core scaling without changing the algorithm.
- **Warm starting** — naturally supports reusing previous solutions, dramatically improving stability for resting contacts.
- **Deterministic** — given the same inputs and iteration count, produces identical results. Critical for networking and replay.
- **Handles inequality constraints** — the projection step (clamping) naturally handles non-penetration and friction limits.

**Cons:**
- **Approximate** — the solution is not exact. With few iterations, joints stretch and stacks are soft. This is why the docs say "the Box3D solver is approximate so the joints can be soft in some cases."
- **Iteration-dependent quality** — more iterations = stiffer = more expensive. There is no "correct" answer; you trade accuracy for performance.
- **No global coupling** — because constraints are solved locally (one at a time), long chains of constraints propagate information slowly. A 20-link chain needs ~20 sub-steps for a signal to travel end-to-end.
- **Not suitable for soft body** — the velocity-impulse formulation assumes rigid bodies. Deformable objects require different solvers (PBD, FEM, mass-spring).
- **Position drift** — velocity-level solving doesn't directly fix position errors. The separate position correction stage is an approximation and can introduce energy.

### Gyroscopic Torque (3D-Specific)

A unique challenge in 3D (vs 2D) is the **gyroscopic effect**: a spinning body with non-spherical inertia experiences torque from its own rotation. This is why a spinning top stays upright and why a spinning rifle bullet follows a stable trajectory.

The equation is nonlinear:

```
I · (ω₂ - ω₁) + h · (ω₂ × I·ω₂) = 0
```

where `I` is the inertia tensor, `ω` is angular velocity, and `×` is the cross product. Box3D solves this using **Newton-Raphson iteration** in local coordinates (where the inertia tensor is constant), with a Jacobian derived by Erin Catto. This improves stability for long skinny bodies (rods, bars) that would otherwise tumble unrealistically.

This is one of the key differences between Box3D and 2D engines like Box2D — in 2D, angular velocity is a scalar and there is no gyroscopic effect.

---

## Connection to Molecular Dynamics and Other Solvers

If you come from molecular dynamics (MD), particle simulation, or position-based methods (PBD, Projective Dynamics, VBD), the relationship to Box3D's solver may not be immediately obvious. This section bridges that gap.

### The MD Approach: Explicit Force Integration

In molecular dynamics, you integrate Newton's equations of motion directly. For a system of particles with pairwise interactions (bonds, Lennard-Jones, etc.), a typical velocity Verlet / leapfrog step looks like:

```
1. Compute forces:    f_i = -dU/dr_i          (gradient of potential energy)
2. Update velocities: v_i += f_i * dt / m_i    (explicit Euler step)
3. Damping:           v_i *= clamp(0, 1 - dt*γ)  (velocity damping)
4. Update positions:  r_i += v_i * dt          (leapfrog / velocity Verlet)
```

For a harmonic bond between particles `i` and `j` with rest length `L` and stiffness `k`:

```
U = 1/2 * k * (|r_ij| - L)²
f_i = -k * (|r_ij| - L) * r̂_ij
```

The force pushes particles toward the rest length. The stiffer the bond (higher `k`), the faster it corrects deviations — but also the smaller the time step must be. The stability condition for explicit integration of a spring is:

```
dt < 2 / sqrt(k / m)     (i.e. dt < 2/ω where ω = sqrt(k/m) is the natural frequency)
```

For a "hard" constraint (rigid bond, non-penetration), you need `k → ∞`, which means `dt → 0`. This is the fundamental problem: **explicit integration cannot handle stiff constraints efficiently.** You need `dt ~ 1 fs` for bonds in MD, which is why MD simulations are computationally expensive.

### How Box3D Differs: Implicit Velocity Constraints

Box3D does not compute forces from potentials and integrate them explicitly. Instead, it works at the **velocity level** and solves constraints **implicitly**.

For the same bond (distance joint) between bodies `A` and `B`, Box3D formulates a velocity constraint:

```
J · v = 0
```

where `J = r̂_ijᵀ` extracts the relative velocity along the bond direction, and `v = [v_A, v_B]ᵀ` is the combined velocity vector. The constraint says: "the relative velocity along the bond must be zero" (the bond is not stretching or compressing).

Instead of computing a force and hoping it's large enough, Box3D computes the **exact impulse** needed to satisfy this constraint in one shot:

```
λ = -(J · v) / (J · M⁻¹ · Jᵀ)
```

where `M⁻¹` is the diagonal inverse mass matrix. The denominator `J · M⁻¹ · Jᵀ` is the **effective mass** of the constraint — it's the 1D equivalent of the Hessian / stiffness matrix `H` you mentioned. For a distance constraint between two equal-mass particles:

```
J · M⁻¹ · Jᵀ = 1/m_A + 1/m_B = 2/m    (for equal masses)
```

The impulse `λ` is then applied:

```
v_A -= (λ / m_A) * r̂_ij
v_B += (λ / m_B) * r̂_ij
```

This is **exactly the velocity change that a stiff spring would produce** in the limit `k → ∞` — but computed directly, without oscillation, without tiny time steps. The constraint is satisfied exactly for this pair, regardless of `dt`.

### The Key Insight: Force vs. Impulse

The relationship between MD and SI is best understood through the force-impulse connection:

| | MD (explicit spring) | Box3D (sequential impulse) |
|---|---|---|
| **What is computed** | Force `f = -k·(d - L)·r̂` | Impulse `λ = -(J·v) / (J·M⁻¹·Jᵀ)` |
| **How it's applied** | `v += f·dt/m` (explicit) | `v += λ/m` (implicit, exact) |
| **Stability** | Conditional: `dt < 2/ω` | Unconditional: stable for any `dt` |
| **Accuracy** | Exact for soft springs | Exact for single constraint, approximate for coupled constraints |
| **Stiffness** | Limited by `dt` | No limit — constraint is satisfied exactly |
| **Energy** | Conserved (symplectic) | May drift (dissipative due to clamping) |

The impulse `λ` is related to the MD force by `λ = f · dt`. But instead of computing `f` from a potential gradient (which requires small `dt` for stability), SI computes `λ` directly from the velocity constraint (which is stable for any `dt`).

In the linear approximation you mentioned — where the Hessian stiffness matrix `H = d²U/dr²` approximates the system — the SI denominator `J·M⁻¹·Jᵀ` is exactly the **projected stiffness** of the constraint. For a harmonic potential, `H = k` and `J·M⁻¹·Jᵀ = k_eff` where `k_eff` accounts for the masses of both bodies. SI solves the constraint using this effective stiffness, but does so one constraint at a time rather than assembling the full Hessian.

### Jacobi vs. Gauss-Seidel: The MD Connection

Your intuition is correct: **MD with explicit integration is essentially a Jacobi solver**.

In MD, all forces are computed from the current positions simultaneously, then all velocities are updated at once:

```
for each pair (i,j):  compute f_ij from current r    ← all using same positions
for each particle i:  v_i += (Σ_j f_ij) * dt / m_i   ← all updated simultaneously
for each particle i:  r_i += v_i * dt
```

This is **Jacobi**: every interaction uses the state from the *beginning* of the step. No interaction sees the result of another within the same step. Information propagates at one constraint per time step.

**Sequential Impulses is Gauss-Seidel**: constraints are solved one at a time, and each constraint sees the *updated* velocities from all previously solved constraints:

```
for each constraint (i,j):
    compute λ from current v_i, v_j     ← uses already-updated velocities
    apply λ to v_i, v_j immediately     ← next constraint sees this update
```

Information propagates faster: within one iteration, a change to body `A` is immediately visible to the next constraint involving `A`. This is why SI converges faster than Jacobi for the same iteration count.

**The tradeoff:** Jacobi (MD-style) is trivially parallel — all force computations are independent. Gauss-Seidel (SI) is inherently sequential — each constraint depends on the previous. Box3D's graph coloring is the bridge: within a color, constraints are independent (Jacobi-like, parallelizable); across colors, they're sequential (Gauss-Seidel-like, better convergence).

### Why MD Needs Tiny Steps but Box3D Doesn't

The fundamental difference is **explicit vs. implicit** treatment of constraints:

**MD (explicit):** The constraint force is computed from the *current* position error, then integrated forward. If the position error is large, the force is large, which causes a large velocity change, which causes overshoot, which causes a larger error next step → instability. The stability limit is `dt < 2/ω`.

**Box3D (implicit):** The constraint impulse is computed to satisfy the velocity constraint *exactly*. There is no overshoot because the impulse is not proportional to the error — it's the exact amount needed to zero the relative velocity. The constraint is satisfied regardless of `dt`.

This is the same reason why implicit methods in numerical ODE solvers (backward Euler, etc.) are unconditionally stable, while explicit methods (forward Euler, Runge-Kutta) have step size limits.

**However**, Box3D's accuracy still depends on `dt` and iteration count:
- Large `dt` with few iterations → constraints are satisfied approximately, joints stretch
- Small `dt` (or more sub-steps) with more iterations → constraints are satisfied more accurately
- The stability is unconditional, but **accuracy is not** — you trade computation for precision

### Connection to PBD, Projective Dynamics, and VBD

You mentioned implementing Projective Dynamics, Position-Based Dynamics (PBD), and Vertex Block Descent (VBD). These are all related approaches that share the same fundamental idea: **solve constraints iteratively rather than through explicit force integration**. The differences are in *what space* they work in and *how* they parallelize.

#### Position-Based Dynamics (PBD)

PBD works in **position space** instead of velocity space. Instead of computing impulses to fix velocities, it directly moves points to satisfy constraints:

```
for each constraint (i,j):
    compute position correction Δx = -C(x) / |∇C|² * ∇C
    move x_i and x_j by Δx (weighted by inverse masses)
```

where `C(x)` is the constraint function (e.g., `C = |r_ij| - L` for a distance constraint). After all constraints are solved, velocity is derived from the position change: `v = (x_new - x_old) / dt`.

**Relationship to Box3D:** Both are iterative constraint solvers. PBD projects positions; Box3D projects velocities. PBD's position projection is equivalent to Box3D's position correction (relax) stage, but PBD does it as the *primary* solve, while Box3D does velocity first, then position correction as a secondary pass.

**PBD's stiffness problem:** PBD's stiffness depends on iteration count and time step — more iterations = stiffer. This is similar to Box3D's sub-stepping. XPBD fixes this by adding compliance parameters, making stiffness independent of `dt` and iteration count — similar to how Box3D's `b3MakeSoft` provides time-step-independent softness.

#### Projective Dynamics

Projective Dynamics generalizes PBD with a **local-global alternating solve**:

- **Local step:** For each constraint, compute the "ideal" position that satisfies it locally (projection). This is like computing the impulse for one constraint in SI.
- **Global step:** Assemble all the local projections into a linear system `(M/dt² + H) · x = ...` and solve it. This finds the best compromise across all constraints simultaneously.

The matrix `H` is the Hessian / stiffness matrix — exactly the `H = d²U/dr²` you mentioned. `M/dt²` is the mass matrix scaled by time step.

**Relationship to Box3D:** Projective Dynamics' global step is what Box3D *doesn't* do. Box3D never assembles a global system — it applies constraints one at a time (Gauss-Seidel). Projective Dynamics assembles and solves the full system, which is more accurate but more expensive (matrix factorization, though the matrix is constant for fixed topology and can be prefactorized).

**The tradeoff:**
- Projective Dynamics: more accurate per iteration (global solve), but each iteration is expensive (matrix solve). Good for cloth and soft bodies with fixed connectivity.
- Box3D (SI): less accurate per iteration (local only), but each iteration is cheap (no matrix). Good for rigid bodies with changing contact topology (contacts appear and disappear every frame).

This is why Box3D doesn't use Projective Dynamics: **contact topology changes every frame**. The constraint graph is different each step as objects collide and separate. Assembling and factorizing a new matrix every frame would be prohibitively expensive. SI's constraint-by-constraint approach handles dynamic topology naturally.

#### Vertex Block Descent (VBD)

VBD is a parallel variant of PBD that addresses the sequential nature of Gauss-Seidel. It processes vertices in **blocks** — groups of vertices that can be updated independently because they don't share constraints. Within a block, all vertices are updated simultaneously (Jacobi-like). Across blocks, updates are sequential (Gauss-Seidel-like).

**Relationship to Box3D:** VBD's block decomposition is conceptually identical to Box3D's **graph coloring**. Both partition constraints into independent groups that can be processed in parallel:

| | VBD | Box3D graph coloring |
|---|---|---|
| **Space** | Position space | Velocity space |
| **Partition unit** | Vertex blocks | Graph colors |
| **Within partition** | Parallel (Jacobi) | Parallel (Jacobi) |
| **Across partitions** | Sequential (Gauss-Seidel) | Sequential (Gauss-Seidel) |
| **Target** | Cloth / deformable | Rigid body contacts |

The core idea is the same: get Gauss-Seidel-like convergence with Jacobi-like parallelism by partitioning into independent groups.

### Summary: The Unified View

All these methods are iterative solvers for constrained dynamics. They differ in three axes:

**1. What space do they work in?**
- **MD (explicit):** Force space → integrate to velocity → integrate to position. Stable only for soft constraints.
- **Box3D (SI):** Velocity space → compute impulses → integrate to position. Stable for hard constraints, approximate for coupled systems.
- **PBD / XPBD:** Position space → project positions → derive velocities. Stable for hard constraints, stiffness depends on iterations (XPBD fixes this).
- **Projective Dynamics:** Position space → local-global alternating solve. More accurate per iteration, but needs matrix solve.

**2. How do they parallelize?**
- **MD:** Trivially parallel (all forces independent) — pure Jacobi.
- **Box3D:** Graph coloring — Jacobi within colors, Gauss-Seidel across colors. Multi-core via atomic CAS on solver blocks.
- **PBD:** Sequential (Gauss-Seidel) or parallel (Jacobi, less stable).
- **VBD:** Block decomposition — same idea as graph coloring.
- **Projective Dynamics:** Global solve is parallelizable (prefactorized matrix), local step is embarrassingly parallel.

**3. How do they handle dynamic topology?**
- **MD:** Topology is usually fixed (bonds don't break). Can precompute neighbor lists.
- **Box3D:** Topology changes every frame (contacts created/destroyed). Must handle dynamically. SI handles this naturally — no matrix to reassemble.
- **PBD / Projective Dynamics:** Typically used with fixed topology (cloth, soft body). Projective Dynamics' prefactorized matrix assumes fixed connectivity.
- **VBD:** Can handle dynamic topology but is designed for cloth/meshes.

**The fundamental equation connecting all approaches:**

For a constraint `C(x) = 0` with stiffness `k`, the correction can be written as:

```
MD:          f = -k · C(x) · ∇C          → v += f·dt/m    (explicit, unstable for large k)
SI:          λ = -Ċ(v) / (J·M⁻¹·Jᵀ)      → v += λ·Jᵀ/M   (implicit velocity, stable)
PBD:         Δx = -C(x) / |∇C|² · ∇C      → x += Δx        (implicit position, stable)
ProjDyn:     (M/dt² + H) · x = ...         → solve globally  (implicit position, exact per iteration)
```

All four are approximations of the same physics. The difference is *what* you solve for (force, impulse, position, or global position) and *how* you iterate (Jacobi, Gauss-Seidel, graph-colored, or global solve). Box3D chose velocity-space impulses with graph coloring because it handles hard contacts with dynamic topology efficiently on multi-core CPUs.

---

## Component Architecture

Below are the key software components, what they do, and how they interconnect.

---

## 1. Core & Platform Layer

**Files:** `@/home/prokophapala/git_SW/box3d/src/core.h`, `@/home/prokophapala/git_SW/box3d/src/core.c`, `@/home/prokophapala/git_SW/box3d/src/platform.h`

The foundation: platform detection (Windows/Linux/macOS/WASM), SIMD selection (SSE2/NEON), atomic primitives, mutex/semaphore/thread abstractions, and global allocation hooks.

Key types:
```c
// core.h:101-109
typedef struct b3AtomicInt { int value; } b3AtomicInt;
typedef struct b3AtomicU32 { uint32_t value; } b3AtomicU32;
```

Global allocators are hookable:
```c
// core.h:127-130
void* b3Alloc( size_t size );
void* b3AllocZeroed( size_t size );
void b3Free( void* mem, size_t size );
```

---

## 2. Container & ID System

**Files:** `@/home/prokophapala/git_SW/box3d/src/container.h`, `@/home/prokophapala/git_SW/box3d/src/id_pool.h`, `@/home/prokophapala/git_SW/box3d/src/id_pool.c`

**`b3Array(T)`** is a macro-based generic dynamic array (like `std::vector<T>`). Every major data structure in the engine uses it:

```c
// container.h:12-18
#define b3DeclareArray( T ) \
	typedef struct b3DynamicArray_##T { struct T* data; int count; int capacity; } b3DynamicArray_##T
```

**`b3IdPool`** provides stable integer IDs with a free list. When you destroy a body/shape/joint, its ID is recycled:

```c
// id_pool.h:8-12
typedef struct b3IdPool {
	b3Array( int ) freeArray;
	int nextIndex;
} b3IdPool;
```

The world maintains separate pools for bodies, shapes, joints, contacts, and islands. User-facing IDs (`b3BodyId`, `b3ShapeId`, etc.) wrap these raw indices with a generation counter for dangling-pointer safety.

---

## 3. Math & SIMD

**Files:** `@/home/prokophapala/git_SW/box3d/src/math_internal.h`, `@/home/prokophapala/git_SW/box3d/src/math_functions.c`, `@/home/prokophapala/git_SW/box3d/src/simd.h`, `@/home/prokophapala/git_SW/box3d/src/simd.c`

Public math types (`b3Vec3`, `b3Quat`, `b3Transform`, `b3Matrix3`, `b3AABB`) live in [include/box3d/math_functions.h](cci:7://file:///home/prokophapala/git_SW/box3d/include/box3d/math_functions.h:0:0-0:0). Internal math ([math_internal.h](cci:7://file:///home/prokophapala/git_SW/box3d/src/math_internal.h:0:0-0:0)) adds helpers for solver use.

SIMD is auto-selected at compile time: **SSE2** on x86, **NEON** on ARM, with a scalar fallback. The SIMD layer processes 4 contacts in parallel in the solver's "wide" contact constraints.

---

## 4. World (`b3World`)

**Files:** `@/home/prokophapala/git_SW/box3d/src/physics_world.h`, `@/home/prokophapala/git_SW/box3d/src/physics_world.c`

The **central hub**. Everything lives inside `b3World`:

```c
// physics_world.h:126-286
typedef struct b3World {
	b3Stack stack;                    // scratch memory for step
	b3BroadPhase broadPhase;          // AABB tree + pair finding
	b3ConstraintGraph constraintGraph; // graph-colored constraints

	b3IdPool bodyIdPool;              // stable body IDs
	b3Array( b3Body ) bodies;         // sparse: ID -> body metadata

	b3IdPool solverSetIdPool;
	b3Array( b3SolverSet ) solverSets; // static/awake/disabled/sleeping

	b3IdPool jointIdPool;
	b3Array( b3Joint ) joints;

	b3IdPool contactIdPool;
	b3Array( b3Contact ) contacts;

	b3IdPool islandIdPool;
	b3Array( b3Island ) islands;

	b3IdPool shapeIdPool;
	b3Array( b3Shape ) shapes;

	void* hullDatabase;               // ref-counted shared hull data

	b3Array( b3Sensor ) sensors;
	b3Array( b3TaskContext ) taskContexts; // per-thread storage

	// Events
	b3Array( b3BodyMoveEvent ) bodyMoveEvents;
	b3Array( b3ContactBeginTouchEvent ) contactBeginEvents;
	// ...end events are double-buffered

	b3Vec3 gravity;
	uint16_t generation;
	// ...scheduler, callbacks, recording, etc.
} b3World;
```

The main entry point is **[b3World_Step](cci:1://file:///home/prokophapala/git_SW/box3d/src/physics_world.c:1023:0-1194:1)**:

```c
// physics_world.c:1024-1195
void b3World_Step( b3WorldId worldId, float timeStep, int subStepCount )
{
	// 1. Update broad-phase pairs -> create new contacts
	b3UpdateBroadPhasePairs( world );

	// 2. Narrow phase: collide all awake contacts in parallel
	b3Collide( &context );

	// 3. Solve: integrate velocities, solve constraints, integrate positions
	b3Solve( world, &context );

	// 4. Update sensors
	b3OverlapSensors( world );

	// 5. Swap end-event double buffer, record state hash
}
```

---

## 5. Broad Phase (`b3BroadPhase`)

**Files:** `@/home/prokophapala/git_SW/box3d/src/broad_phase.h`, `@/home/prokophapala/git_SW/box3d/src/broad_phase.c`, `@/home/prokophapala/git_SW/box3d/src/dynamic_tree.c`

Uses a **`b3DynamicTree`** (bounding volume hierarchy / BVH) per body type:

```c
// broad_phase.h:27-47
typedef struct b3BroadPhase {
	b3DynamicTree trees[b3_bodyTypeCount];  // separate trees: static, kinematic, dynamic
	b3BitSet movedProxies[b3_bodyTypeCount];
	b3Array( int ) moveArray;
	b3MoveResult* moveResults;
	b3MovePair* movePairs;
	b3HashSet pairSet;  // tracks existing shape pairs with contacts
} b3BroadPhase;
```

When a body moves, [b3BufferMove](cci:1://file:///home/prokophapala/git_SW/box3d/src/broad_phase.h:66:0-78:1) marks its proxy. Each step, `b3UpdateBroadPhasePairs` queries moved proxies against all trees to find new overlapping pairs, which triggers contact creation.

The **`b3DynamicTree`** is a binary AABB tree with category-bit filtering for fast ray casts, region queries, and box casts. It's also usable standalone for game spatial queries.

---

## 6. Bodies (`b3Body` / `b3BodySim` / `b3BodyState`)

**Files:** `@/home/prokophapala/git_SW/box3d/src/body.h`, `@/home/prokophapala/git_SW/box3d/src/body.c`

Bodies are **split into three structs** for cache efficiency:

- **`b3Body`** — organizational metadata: ID, set index, linked-list heads for shapes/contacts/joints, mass, inertia, sleep state, flags. This is the "slow path" struct.

```c
// body.h:75-131
typedef struct b3Body {
	int setIndex;      // which solver set (static/awake/sleeping)
	int localIndex;    // position within that set
	int headContactKey; // linked list of contacts
	int headShapeId;    // linked list of shapes
	int headJointKey;   // linked list of joints
	int islandId;
	float mass;
	b3Matrix3 inertia;
	b3BodyType type;    // static, kinematic, dynamic
	uint16_t generation;
	// ...
} b3Body;
```

- **`b3BodySim`** — simulation data used for integration: transform, center of mass, force, torque, inverse mass/inertia. Lives in the solver set's contiguous `bodySims` array.

- **`b3BodyState`** — the hot solver data: velocities and delta position/rotation. Only exists for awake bodies. Designed for SIMD scatter-gather (56 bytes).

```c
// body.h:161-176
typedef struct b3BodyState {
	b3Vec3 linearVelocity;   // 12
	b3Vec3 angularVelocity;  // 12
	b3Vec3 deltaPosition;    // 12  (delta to reduce float error far from origin)
	b3Quat deltaRotation;    // 16  (delta rotation, identity for static)
	uint32_t flags;          // 4
} b3BodyState;
```

---

## 7. Shapes (`b3Shape`)

**Files:** `@/home/prokophapala/git_SW/box3d/src/shape.h`, `@/home/prokophapala/git_SW/box3d/src/shape.c`

Shapes bind collision geometry to a body with material properties. They form a doubly-linked list per body (`prevShapeId`/`nextShapeId`).

```c
// shape.h:25-64
typedef struct b3Shape {
	int id;
	int bodyId;
	int prevShapeId, nextShapeId;  // body's shape linked list
	int proxyKey;                   // broad-phase proxy
	b3ShapeType type;               // sphere, capsule, hull, mesh, heightField, compound
	float density;
	b3AABB aabb, fatAABB;           // current and inflated AABB
	int materialCount;
	b3SurfaceMaterial material;     // inline single material
	b3SurfaceMaterial* materials;   // heap array for multi-material
	b3Filter filter;
	union {
		b3Capsule capsule;
		b3Sphere sphere;
		const b3HullData* hull;
		b3Mesh mesh;
		const b3HeightFieldData* heightField;
		const b3CompoundData* compound;
	};
	// ...
} b3Shape;
```

Shape types: **sphere**, **capsule**, **convex hull**, **triangle mesh** (static only, has internal BVH), **height field** (static only, regular grid), **compound** (static only, baked aggregation).

The hull database (`world->hullDatabase`) deduplicates identical hull data by content hash, sharing one copy with reference counting.

---

## 8. Contacts (`b3Contact`)

**Files:** `@/home/prokophapala/git_SW/box3d/src/contact.h`, `@/home/prokophapala/git_SW/box3d/src/contact.c`

Contacts represent the **persistent interaction between two shapes**. They're created by the broad phase and updated each step by the narrow phase.

```c
// contact.h:96-158
typedef struct b3Contact {
	int setIndex;      // which solver set
	int colorIndex;    // graph color (for parallel solving)
	int localIndex;    // index within color
	b3ContactEdge edges[2];  // doubly-linked list per body
	int shapeIdA, shapeIdB;
	int islandId;
	b3Manifold* manifolds;
	int manifoldCount;
	b3Quat cachedRotationA, cachedRotationB;  // for contact recycling
	float friction, restitution, rollingResistance;
	union {
		b3ConvexContact convexContact;  // SAT/simplex cache
		b3MeshContact meshContact;      // per-triangle cache
	};
	// ...
} b3Contact;
```

Contact lifecycle per step:
1. **Broad phase** finds overlapping AABB pairs → `b3CreateContact`
2. **Narrow phase** ([b3CollideTask](cci:1://file:///home/prokophapala/git_SW/box3d/src/physics_world.c:549:0-778:1)) runs the appropriate collide function (e.g., `b3CollideHulls`, `b3CollideSpheres`) to compute manifolds
3. State changes are processed serially: started touching → link to island + add to constraint graph; stopped touching → unlink + remove from graph; disjoint → destroy

---

## 9. Constraint Graph (`b3ConstraintGraph`)

**Files:** `@/home/prokophapala/git_SW/box3d/src/constraint_graph.h`, `@/home/prokophapala/git_SW/box3d/src/constraint_graph.c`

**Graph coloring** for parallel constraint solving. Constraints (contacts + joints) sharing a body get different colors. Constraints of the same color can be solved simultaneously because they touch disjoint bodies.

```c
// constraint_graph.h:39-67
typedef struct b3GraphColor {
	b3BitSet bodySet;                    // which bodies are in this color
	b3Array( b3JointSim ) jointSims;
	b3Array( int ) convexContacts;
	b3Array( b3ContactSpec ) contacts;
	struct b3ContactConstraintWide* wideConstraints;  // SIMD-packed
	int wideConstraintCount;
	struct b3ManifoldConstraint* manifoldConstraints; // mesh/overflow
	// ...
} b3GraphColor;

typedef struct b3ConstraintGraph {
	b3GraphColor colors[B3_GRAPH_COLOR_COUNT];  // ~32 colors + overflow
} b3ConstraintGraph;
```

The last color is **overflow** for when a single body touches too many others (exceeds color count). Overflow constraints are solved serially.

---

## 10. Solver (`b3Solve` / `b3StepContext`)

**Files:** `@/home/prokophapala/git_SW/box3d/src/solver.h`, `@/home/prokophapala/git_SW/box3d/src/solver.c`, `@/home/prokophapala/git_SW/box3d/src/contact_solver.h`, `@/home/prokophapala/git_SW/box3d/src/contact_solver.c`

The solver runs **sequential impulses** with sub-stepping. Per time step:

```c
// solver.h:69-82  — solver stages
typedef enum b3SolverStageType {
	b3_stagePrepareJoints,
	b3_stagePrepareWideContacts,
	b3_stagePrepareContacts,
	b3_stageIntegrateVelocities,
	b3_stageWarmStart,
	b3_stageSolve,           // iterative: repeated per sub-step
	b3_stageIntegratePositions,
	b3_stageRelax,           // position stabilization
	b3_stageRestitution,
	b3_stageStoreWideImpulses,
	b3_stageStoreImpulses,
} b3SolverStageType;
```

Work is partitioned into **`b3SolverBlock`s** claimed by workers via atomic CAS on a per-block `syncIndex`. This avoids the cache-line stampede of a single shared counter:

```c
// solver.h:112-127
typedef struct b3SolverBlock {
	int startIndex;
	uint16_t count;
	uint8_t blockType;    // body/joint/wideContact/contact/graphJoint/...
	uint8_t colorIndex;
} b3SolverBlock;

typedef struct b3SyncBlock {
	b3SolverBlock block;
	b3AtomicInt syncIndex;  // monotonic across sub-step iterations
} b3SyncBlock;
```

**Contact constraints** come in two flavors:
- **`b3ContactConstraintWide`** — SIMD-packed, 4 contacts per wide block, used for convex-convex contacts
- **`b3ManifoldConstraint`** — scalar, used for mesh contacts and overflow

The **softness** model uses spring-damper parameters (Hertz + damping ratio) converted to bias/mass/impulse scales:

```c
// solver.h:264-306
static inline b3Softness b3MakeSoft( float hertz, float zeta, float h ) {
	float omega = 2.0f * B3_PI * hertz;
	// ... returns biasRate, massScale, impulseScale
}
```

---

## 11. Solver Sets (`b3SolverSet`)

**Files:** `@/home/prokophapala/git_SW/box3d/src/solver_set.h`, `@/home/prokophapala/git_SW/box3d/src/solver_set.c`

Bodies and constraints are grouped into sets for **memory locality**:

```c
// solver_set.h:30-55
typedef struct b3SolverSet {
	b3Array( b3BodySim ) bodySims;
	b3Array( b3BodyState ) bodyStates;  // only in awake set
	b3Array( b3JointSim ) jointSims;    // sleeping/disabled joints
	b3Array( int ) contactIndices;      // non-touching contacts in awake set
	b3Array( b3IslandSim ) islandSims;
	int setIndex;
} b3SolverSet;
```

Set indices (from `physics_world.h:47-53`):
```c
enum b3SetType {
	b3_staticSet = 0,       // all static bodies
	b3_disabledSet = 1,     // disabled bodies
	b3_awakeSet = 2,        // active simulation
	b3_firstSleepingSet = 3, // sleeping islands, one set each
};
```

When an island sleeps, its data is copied from the awake set into a new sleeping set. When woken, it's merged back. This keeps the awake set small and cache-friendly.

---

## 12. Islands (`b3Island`)

**Files:** `@/home/prokophapala/git_SW/box3d/src/island.h`, `@/home/prokophapala/git_SW/box3d/src/island.c`

Islands are **connected components** of the body-contact-joint graph. Each awake island can be solved independently. When bodies come to rest, islands can be **split** to allow smaller pieces to sleep:

```c
// island.h:48-75
typedef struct b3Island {
	int setIndex;
	int localIndex;
	int islandId;
	int constraintRemoveCount;  // triggers split consideration
	b3Array( int ) bodies;
	b3Array( b3ContactLink ) contacts;  // cached for fast union-find
	b3Array( b3JointLink ) joints;
} b3Island;
```

`b3SplitIsland` uses **union-find** on the contact/joint links to find new connected components after a constraint is removed.

---

## 13. Joints

**Files:** `@/home/prokophapala/git_SW/box3d/src/joint.h`, `@/home/prokophapala/git_SW/box3d/src/joint.c`, plus per-type files: [distance_joint.c](cci:7://file:///home/prokophapala/git_SW/box3d/src/distance_joint.c:0:0-0:0), [revolute_joint.c](cci:7://file:///home/prokophapala/git_SW/box3d/src/revolute_joint.c:0:0-0:0), [prismatic_joint.c](cci:7://file:///home/prokophapala/git_SW/box3d/src/prismatic_joint.c:0:0-0:0), [spherical_joint.c](cci:7://file:///home/prokophapala/git_SW/box3d/src/spherical_joint.c:0:0-0:0), [weld_joint.c](cci:7://file:///home/prokophapala/git_SW/box3d/src/weld_joint.c:0:0-0:0), [wheel_joint.c](cci:7://file:///home/prokophapala/git_SW/box3d/src/wheel_joint.c:0:0-0:0), [motor_joint.c](cci:7://file:///home/prokophapala/git_SW/box3d/src/motor_joint.c:0:0-0:0), [parallel_joint.c](cci:7://file:///home/prokophapala/git_SW/box3d/src/parallel_joint.c:0:0-0:0)

Joints connect two bodies with constraints that remove degrees of freedom. Each joint type has its own `b3JointSim` variant with accumulated impulses. Joints are added to the constraint graph and solved in the same stage pipeline as contacts.

```c
// joint.h:28-62
typedef struct b3Joint {
	int setIndex;
	int colorIndex;     // graph color
	int localIndex;
	b3JointEdge edges[2];  // doubly-linked list per body
	b3JointType type;      // distance, revolute, prismatic, spherical, weld, wheel, motor, parallel, filter
	bool collideConnected;
	uint16_t generation;
} b3Joint;
```

---

## 14. Narrow Phase Collision

**Files:** `@/home/prokophapala/git_SW/box3d/src/distance.c` (GJK), `@/home/prokophapala/git_SW/box3d/src/convex_manifold.c` (SAT for hulls), `@/home/prokophapala/git_SW/box3d/src/triangle_manifold.c` (mesh contacts), `@/home/prokophapala/git_SW/box3d/src/mesh_contact.c`, `@/home/prokophapala/git_SW/box3d/src/capsule.c`, `@/home/prokophapala/git_SW/box3d/src/sphere.c`

These implement the actual collision detection algorithms:
- **GJK** (`b3ShapeDistance`) — closest points between convex shapes via simplex iteration
- **SAT** (`b3CollideHulls`) — separating axis theorem for hull-hull, with warm-started `b3SATCache`
- **TOI** (`b3TimeOfImpact`) — continuous collision detection via sweep advancement
- **Mesh contacts** — per-triangle collision with cached triangle indices
- **Mover** ([mover.c](cci:7://file:///home/prokophapala/git_SW/box3d/src/mover.c:0:0-0:0)) — capsule-based character controller collision

---

## 15. Compound Shapes

**Files:** `@/home/prokophapala/git_SW/box3d/src/compound.c`, `@/home/prokophapala/git_SW/box3d/src/compound.h`

Compounds aggregate multiple sub-shapes (capsules, hulls, meshes, spheres) into a single static collision shape. Designed for **offline baking** and world streaming. Material deduplication uses a hash map (`b3MaterialMap`), and hull/mesh sharing uses content-addressed databases.

---

## 16. Sensors

**Files:** `@/home/prokophapala/git_SW/box3d/src/sensor.h`, `@/home/prokophapala/git_SW/box3d/src/sensor.c`

Sensors are shapes that detect overlap without generating contact constraints. Each sensor maintains double-buffered visitor lists (`overlaps1`/`overlaps2`) to detect enter/exit events:

```c
// sensor.h:27-33
typedef struct b3Sensor {
	b3Array( b3Visitor ) hits;       // TOI-based hits
	b3Array( b3Visitor ) overlaps1;  // current step overlaps
	b3Array( b3Visitor ) overlaps2;  // previous step overlaps
	int shapeId;
} b3Sensor;
```

---

## 17. Recording & Replay

**Files:** `@/home/prokophapala/git_SW/box3d/src/recording.h`, `@/home/prokophapala/git_SW/box3d/src/recording.c`, `@/home/prokophapala/git_SW/box3d/src/recording_replay.h`, `@/home/prokophapala/git_SW/box3d/src/recording_replay.c`, `@/home/prokophapala/git_SW/box3d/src/recording_ops.inl`

A deterministic recording system. Every world mutation (body creation, shape addition, joint creation, step) is logged as an op stream. The recording includes a registry of geometry (hulls, meshes) and per-step state hashes for validation. This enables **deterministic replay** and cross-platform determinism testing.

---

## 18. Multithreading

**Files:** `@/home/prokophapala/git_SW/box3d/src/scheduler.h`, `@/home/prokophapala/git_SW/box3d/src/scheduler.c`, `@/home/prokophapala/git_SW/box3d/src/parallel_for.h`, `@/home/prokophapala/git_SW/box3d/src/parallel_for.c`, `@/home/prokophapala/git_SW/box3d/src/parallel_joint.c`

Box3D uses **data parallelism** — multiple cores work on the same world step, not task parallelism. The internal scheduler creates worker threads. `b3ParallelFor` splits work into ranges and dispatches to workers. The solver's block-based CAS design (see §10) distributes contention across per-block atomics.

Per-thread storage is in `b3TaskContext`:

```c
// physics_world.h:71-120
typedef struct b3TaskContext {
	b3Arena arena;                    // scratch allocator
	b3Array( b3SensorHit ) sensorHits;
	b3BitSet contactStateBitSet;      // contact state changes
	b3BitSet jointStateBitSet;
	b3BitSet hitEventBitSet;
	b3BitSet enlargedSimBitSet;
	b3BitSet awakeIslandBitSet;
	// profiling counters...
	char cacheLine[64];  // prevent false sharing
} b3TaskContext;
```

---

## 19. Memory Management

**Files:** `@/home/prokophapala/git_SW/box3d/src/arena_allocator.h`, `@/home/prokophapala/git_SW/box3d/src/arena_allocator.c`, `@/home/prokophapala/git_SW/box3d/src/block_allocator.h`, `@/home/prokophapala/git_SW/box3d/src/block_allocator.c`

- **`b3Arena`** — bump allocator for per-thread scratch memory during narrow phase. Reset after each step.
- **`b3BlockAllocator`** — slab allocator for fixed-size blocks (used for manifolds). Pools are per-manifold-count.
- **`b3Stack`** — world-level stack allocator for temporary arrays during a step.

The design minimizes per-frame allocations: after warmup, memory is recycled via free lists and array reuse.

---

## Interconnection Summary

```
User creates World
  └── b3World holds everything
       ├── BroadPhase (3 DynamicTrees)
       │    └── finds overlapping AABB pairs → creates Contacts
       │
       ├── Bodies (sparse array via IdPool)
       │    └── each Body points to SolverSet
       │         └── SolverSet holds contiguous BodySim + BodyState arrays
       │
       ├── Shapes (sparse array)
       │    └── each Shape has a proxy in BroadPhase
       │    └── geometry: sphere/capsule/hull/mesh/heightField/compound
       │
       ├── Contacts (sparse array)
       │    └── Narrow phase updates manifolds each step
       │    └── Touching contacts → added to ConstraintGraph (graph coloring)
       │    └── Non-touching → stored in awake SolverSet
       │
       ├── Joints (sparse array)
       │    └── Added to ConstraintGraph by color
       │
       ├── ConstraintGraph (B3_GRAPH_COLOR_COUNT colors)
       │    └── Each color: disjoint body sets → parallel solving
       │    └── Wide (SIMD) constraints for convex contacts
       │    └── Manifold constraints for mesh/overflow
       │
       ├── Islands (connected components)
       │    └── Awake island → solved together
       │    └── Sleep → moved to sleeping SolverSet
       │    └── Split → union-find on contact/joint links
       │
       └── Step pipeline:
            1. b3UpdateBroadPhasePairs → new contacts
            2. b3Collide (parallel) → update manifolds, state changes
            3. b3Solve (parallel, sub-stepped):
               prepare → integrate velocities → warm start → solve → integrate positions → relax → restitution
            4. b3OverlapSensors → sensor events
            5. Events available to user (begin/end touch, hit, body move)
```