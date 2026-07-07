# USER

look in 
/home/prokop/git/FireCore/pyBall/RigidAtomFF/RRsp3/
/home/prokop/git/FireCore/pyBall/RigidAtomFF/RRsp3/RRsp3.cl
/home/prokop/git/FireCore/pyBall/RigidAtomFF/RRsp3/RRsp3_momentum_design.md

and figure out if you would be able to implement it in pure c here (no OpenCL) and then integrate it with visualization (samples)
@overview.md @visualization.md 


---

# GLM 5.2

Let me read all the referenced files in parallel.

Let me read the rest of the `.cl` file and explore the box3d source structure.

Let me check the existing molecular sample and the box3d math types to understand integration points.

Now I have a complete picture. Here's my assessment:

---

# Feasibility Analysis: Porting RRsp3 to Pure C in Box3D

## Short Answer

**Yes, this is very feasible.** The RRsp3 OpenCL code is straightforward numerical computation that translates cleanly to C. The box3d visualization system already has the infrastructure to display it, and there's even a precedent ([sample_molecules.cpp](cci:7://file:///home/prokop/git/box3d/samples/sample_molecules.cpp:0:0-0:0)).

## What RRsp3 Does (from the `.cl` and design doc)

The system is an **XPBD rigid-body molecular dynamics solver** with:

- **Atoms** as rigid bodies: position (`float4` with invMass in `.w`), quaternion, linear velocity, angular velocity
- **Port constraints**: each atom has up to 4 "ports" (local-space lever arms). A bond is enforced by constraining the rotated port tip to the neighbor atom's center
- **Collisions**: sphere-sphere Jacobi positional correction with exclusion lists
- **Multiple solver variants**: massfull (physical rotational inertia), massless (geometric alignment via Newton/Kabsch/eigen)
- **Heavy-ball momentum**: accelerates Jacobi convergence within a timestep
- **Cluster-sorted topology with ghost atoms**: a GPU-specific optimization for workgroup-local memory

## Translation: OpenCL → C

The port is **mechanical** because the code has no GPU-specific algorithms, only GPU-specific memory management:

| OpenCL construct | C equivalent | Difficulty |
|---|---|---|
| `__kernel` functions | Regular C functions called in a loop | Trivial |
| `get_global_id(0)` | Loop index `i` | Trivial |
| `__local` memory + `barrier()` | **Eliminated** — on CPU all atoms are directly accessible | Simplification |
| Ghost atoms / cluster sorting | **Eliminated** — no workgroup memory walls | Major simplification |
| `float4` / `float3` / `cross` / `dot` | `b3Vec4` / `b3Vec3` / `b3Cross` / `b3Dot` (already in box3d) | Trivial |
| `atomic_inc` | Not needed (single-threaded) | Eliminated |

**The ghost atom machinery (~200 lines of `build_local_topology_rigid`) disappears entirely.** On CPU, neighbor lookups are just `pos[j]` — no mapping from global to local indices.

## What Needs to Be Implemented

### 1. Core solver (~600 lines of C)
Translate these kernels to C functions:
- `predict_dynamics` — symplectic Euler position/rotation prediction
- `compute_collision_cluster_rigid` — O(n²) or grid-based sphere collision (no ghosts needed)
- `compute_ports_cluster_rigid` — the core XPBD port constraint (pick one variant; I'd start with `orig`/`current` massfull)
- `apply_corrections_rigid_ports` — accumulate corrections with heavy-ball momentum
- `update_velocities_dynamics` — derive velocities from position deltas

### 2. Data structures (~100 lines)
```c
typedef struct {
    b3Vec4 pos;       // xyz = position, w = invMass
    b3Quat quat;
    b3Vec3 vel;
    b3Vec3 omega;
    b3Vec3 pos_prev;
    b3Quat quat_prev;
    float  radius;
    int    fixmask;
} RRsp3_Atom;

typedef struct {
    int    neighbor[4];   // global atom indices (-1 = unused)
    b3Vec3 port_local[4]; // local-space lever arms
    float  stiffness[4];  // per-port spring stiffness
} RRsp3_Node;
```

### 3. Math helpers (~100 lines)
Most already exist in box3d's `math_functions.h` (`b3Quat`, `b3Vec3`, `b3Cross`, `b3Dot`, etc.). A few need adding:
- `quat_rotate(q, v)` — rotate vector by quaternion
- `quat_from_axis_angle(axis, angle)`
- `quat_delta_rotvec(q_new, q_old)` — extract rotation vector from quaternion delta
- `solve_3x3` — 3×3 linear solve via adjugate
- `mat3_to_quat` — rotation matrix to quaternion
- `apply_K_matrix` — for the eigen variant

### 4. Visualization integration (~200 lines)
Following the [sample_molecules.cpp](cci:7://file:///home/prokop/git/box3d/samples/sample_molecules.cpp:0:0-0:0) pattern:
- New `sample_rrsp3.cpp` inheriting from [Sample](cci:2://file:///home/prokop/git/box3d/samples/sample.h:106:0-229:1)
- [Step()](cci:1://file:///home/prokop/git/box3d/samples/sample_molecules.cpp:223:1-231:2) calls the RRsp3 solver instead of `b3WorldStep`
- [Render()](cci:1://file:///home/prokop/git/box3d/samples/sample.h:121:1-121:25) draws atoms as spheres (`DrawSphere`), bonds as lines (`DrawLine`), port tips as points (`DrawPoint`)
- ImGui controls for: solver variant, iteration count, dt, stiffness, momentum β, collision on/off

## Key Design Decision: Standalone vs. Integrated

There are two approaches:

- **Standalone solver** (recommended): The RRsp3 solver runs its own physics loop independent of box3d's `b3World`. The sample uses box3d **only for visualization** (Draw* calls). This is the cleanest approach — no need to reconcile XPBD with box3d's Sequential Impulses solver. The existing [sample_molecules.cpp](cci:7://file:///home/prokop/git/box3d/samples/sample_molecules.cpp:0:0-0:0) already shows this is possible (it uses box3d's world, but a standalone solver would just skip `b3WorldStep`).

- **Integrated into box3d's solver**: Would require implementing a new joint type or solver set. Much more complex, and the XPBD position-level approach conflicts with box3d's velocity-level SI solver.

## Estimated Effort

| Component | Lines | Effort |
|---|---|---|
| Math helpers | ~100 | 1 hour |
| Data structures + init | ~100 | 1 hour |
| Core solver (5 functions) | ~600 | 4-6 hours |
| Sample + visualization | ~200 | 2 hours |
| CMakeLists integration | ~20 | 30 min |
| Testing/debugging | — | 2-4 hours |
| **Total** | **~1000** | **~1-2 days** |

## Recommendation

Start with the **massfull `orig` variant** — it's the most physically complete and the simplest to implement (single kernel, no iterative substeps or eigen solves). Add the massless variants later as options.

The visualization integration is straightforward since [sample_molecules.cpp](cci:7://file:///home/prokop/git/box3d/samples/sample_molecules.cpp:0:0-0:0) already demonstrates the pattern: load a `.mol` file, create atoms, render spheres + bond lines.

Want me to start implementing?

---

## USER

OK, look in the OpenCL kernles I uploaded, I want to make clear cookbook how to weite this "ports" based aprroach which collisions into efficient scheme using all thre methods projective-dynamics, XPBD, and Sequantial-impulses with soft-steps. MAybe 4th one is vertex-block-descent or augumented vertex block descent. 

I mean the general idea shoul be common - we have quaternion which have some rotation dynamics with torque and angular velocity, we have position of centeral atom, within the local rotation frame given by quaternion we have ~4 ports where this atom want his negibhors to be (this is symmetric other atoms have this also - it si like "as rigid as posible" paper https://igl.ethz.ch/projects/ARAP/arap_web.pdf  Olga Sorkine and Marc Alexa 2007 )

We can have some AABB boudnig boxes from broad space collision and split Leanrd jones or Morse  to appoximative harmonic repulsive part ( Kdr^2 for dr<0, 0 for dr>0, dr=r-R ). This is the common part you do not have to care too much about. Bute care about what would differe in the formuls in the kernels for the 3-4 cases , PD,XPBD,SI,(A)VBD.

Give me discussion and pseudocode, finally try to skerhe maby minimalustic C/C++ code (you do not have to care about OpenCL paralleization local memory and such ocmplication in this conpetucal description) 

also you say you do not need hard constrains or linearization. So you can use sub-steps instead of linear solver iterations. The reason why I wanted to use lienar solver instead of force/impulse substesps on hard-DOFs wha because Projective dynamics provides physically rigorous formulas how to connect the our interial loop and linearized relaxation inner loop via the inertial diagonal regularization term M/dt^2, position based dynamics have something kimilar with 1/Kdt^2 or something like that. 

If you know how to connect sub-steps using hard-DOF forces (O(n)) with outer loop of soft-long range forces O(n^2) in physically rigorous way (without predictor corrector - like extrapolating position in future), I would be very happyly using it instead of linearization. 

---

# Chat GPT 5.5

Yes. Looking at `RRsp3.cl`, the cleanest way to write the cookbook is to separate three things:

1. **The molecular model** — rigid atom centers, quaternions, ports, collisions.
2. **The time-scale splitting** — expensive soft forces once per macrostep, cheap hard interactions in microsteps.
3. **The local solver used for hard interactions** — PD, XPBD, soft-step SI, or VBD/AVBD.

The most important correction to your current kernels is this:

> A bond should preferably be represented as one symmetric **port-to-port constraint**, processed once, rather than two port-to-neighbor-center constraints with heuristic halving.

Your current `compute_ports_cluster_rigid()` uses

[
x_j-\left(x_i+R_i a_{ij}\right),
]

and then guesses the neighbor’s rotational contribution by adding another copy of `w_ang` when `j` is also a node. That is only approximate. The proper symmetric constraint is

[
\boxed{
e_{ij}
======

## \left(x_j+R_j a_{ji}\right)

\left(x_i+R_i a_{ij}\right)
}
]

where (a_{ij}) and (a_{ji}) are the two local ports.

This automatically generates:

* equal and opposite linear reactions;
* torque on atom (i);
* torque on atom (j);
* exact conservation of internal linear and angular momentum;
* no double-counting or `0.5` heuristic.

That symmetric port energy is closely related to an ARAP local-rigidity energy: each atom carries a local frame, and neighboring frames try to map their local reference vectors into mutually compatible world-space geometry. ARAP similarly alternates between local rotations and globally compatible positions. ([Interactive Geometry Lab][1])

---

# 1. Common rigid-port molecular model

Each oriented atom has

```cpp
struct RigidAtom {
    Vec3  x;       // center
    Quat  q;       // orientation
    Vec3  v;       // linear velocity
    Vec3  omega;   // world-space angular velocity

    float invMass;
    Mat3  invInertiaWorld;
};
```

A bond stores two local ports:

```cpp
struct PortBond {
    int ia;
    int ib;

    Vec3 portA_local;
    Vec3 portB_local;

    float stiffness;
    float damping;
    float compliance;

    Vec3 lambda;       // only needed by XPBD/AVBD or warm-started SI
};
```

World-space port arms and tips are

[
r_A=R(q_A)a_A,\qquad
r_B=R(q_B)a_B,
]

[
p_A=x_A+r_A,\qquad
p_B=x_B+r_B.
]

The vector error is

[
e=p_B-p_A.
]

The port velocities are

[
u_A=v_A+\omega_A\times r_A,
]

[
u_B=v_B+\omega_B\times r_B,
]

so the relative port velocity is

[
u=u_B-u_A.
]

Define the skew matrix (S(r)) by

[
S(r)y=r\times y.
]

The (3\times3) inverse effective mass of the point-to-point constraint is

[
\boxed{
K
=

(w_A+w_B)I
+
S(r_A)I_A^{-1}S(r_A)^T
+
S(r_B)I_B^{-1}S(r_B)^T
}
]

where

[
w_i=m_i^{-1}.
]

If the inertia is isotropic,

[
I_i^{-1}=\iota_i I,
]

then

[
S(r)I^{-1}S(r)^T
================

\iota\left(|r|^2I-r r^T\right).
]

This is exactly the matrix version of your present scalar expression

```cpp
dot(cross(r_arm,n), cross(r_arm,n))*invI
```

but it handles all three error directions simultaneously and includes both atoms correctly.

Given a port impulse (P), apply

[
v_A\mathrel{-}=w_AP,
\qquad
\omega_A\mathrel{-}=I_A^{-1}(r_A\times P),
]

[
v_B\mathrel{+}=w_BP,
\qquad
\omega_B\mathrel{+}=I_B^{-1}(r_B\times P).
]

For position-level methods, replace velocity changes by (x) and rotation-vector changes:

[
x_A\mathrel{-}=w_AP,
\qquad
\theta_A\mathrel{-}=I_A^{-1}(r_A\times P),
]

and similarly for (B).

---

# 2. Collision model in the same notation

For atoms (A,B),

[
d=x_B-x_A,\qquad
r=|d|,\qquad
n=\frac{d}{r},
]

[
C=r-R_{AB}.
]

Penetration means

[
C<0.
]

The scalar inverse effective mass is simply

[
k_n=w_A+w_B
]

for spherical atoms with impulses through their centers.

The repulsive energy may be

[
E_{\mathrm{coll}}
=================

\frac12 k_c\min(C,0)^2.
]

Because the force passes through both centers, it generates no torque. All four methods differ only in how they use this same (C), gradient (n), and stiffness (k_c).

---

# 3. Physically principled outer/inner time splitting

You asked for a way to combine cheap hard microsteps with expensive soft forces **without extrapolating future positions**.

The proper framework is multiple-time-step operator splitting, essentially RESPA/Strang splitting.

Let

[
F(x)=F_s(x)+F_h(x),
]

where:

* (F_s): expensive Coulomb, attractive vdW, long-range interactions, torsions;
* (F_h): cheap ports, bonds and collisions.

Use a macrostep (H), divided into

[
s
]

microsteps

[
h=\frac{H}{s}.
]

The symmetric scheme is:

[
v\leftarrow v+\frac{H}{2}M^{-1}F_s(x),
]

then perform (s) hard microsteps, and finally evaluate the soft force at the new position:

[
v\leftarrow v+\frac{H}{2}M^{-1}F_s(x_{\mathrm{new}}).
]

The ending soft force is cached and reused as the starting force of the next macrostep. Thus, after initialization, you still need only **one new expensive soft-force evaluation per macrostep**.

```cpp
// Fs is already evaluated at x_n
kickVelocity(atoms, Fs, 0.5f*H);

for (int isub=0; isub<nsub; isub++) {
    hardMicrostep(atoms, bonds, contacts, h);
}

evaluateSoftForces(atoms, FsNew);       // one new expensive evaluation
kickVelocity(atoms, FsNew, 0.5f*H);

swap(Fs, FsNew);
```

This does not predict or extrapolate (F_s). It evaluates the force only at actual configurations.

For conservative explicit hard forces, this becomes ordinary reversible RESPA. In your case the hard propagator is replaced by a dissipative implicit/constraint solver, so the complete method is no longer symplectic, but the coupling between scales remains consistent.

For relaxation, that loss of reversibility is desirable.

---

# 4. Method A: Projective Dynamics

PD derives the hard step from the implicit Euler variational objective

[
\boxed{
\Phi(q)
=======

\frac{1}{2h^2}
(q-y)^TM(q-y)
+
E_h(q)
}
]

with inertial predictor

[
y=q^n+h\dot q^n.
]

If external hard-step forces are absent, this predictor is not an arbitrary extrapolation. It follows directly from the implicit Euler discretization.

Projective Dynamics introduces auxiliary projected configurations (p_c) and uses energies of the form

[
E_c(q,p_c)
==========

\frac{k_c}{2}
\left|
A_cq-p_c
\right|^2.
]

It alternates:

1. local projection (p_c);
2. global minimization in (q).

The global system is

[
\left(
\frac{M}{h^2}
+
\sum_c k_cA_c^TA_c
\right)q
========

\frac{M}{h^2}y
+
\sum_c k_cA_c^Tp_c.
]

That local-global structure is the defining feature of PD. ([users.cs.utah.edu][2])

## Port constraints in PD

For fixed current rotations, the port bond energy is already quadratic in positions:

[
E_{ij}
======

\frac{k}{2}
\left|
(x_B+r_B)-(x_A+r_A)
\right|^2.
]

The difficult part is rotation. A practical alternating PD/ARAP solve is:

### Local rotational step

For each atom (i), solve

[
R_i^*
=====

\arg\min_{R\in SO(3)}
\sum_{j\in N(i)}
k_{ij}
\left|
R a_{ij}-
t_{ij}
\right|^2,
]

where

[
t_{ij}
======

x_j+R_j a_{ji}-x_i.
]

This is your Kabsch, Davenport eigenvalue, or Newton rotation kernel:

* `compute_ports_cluster_rigid_substep_optimized`
* `compute_ports_cluster_rigid_shapematch`
* `compute_optimal_rotation_eigen`

These are already ARAP-like local rotation projections.

### Global translation step

With (R_i) fixed, minimize

[
\frac{1}{2h^2}\sum_i m_i|x_i-y_i|^2
+
\frac12\sum_{ij}
k_{ij}
\left|
x_j-x_i+r_j-r_i
\right|^2.
]

If you refuse a global Cholesky solve, use matrix-free Jacobi or colored Gauss–Seidel:

[
x_i
\leftarrow
\frac{
m_i y_i/h^2
+
\sum_{j}k_{ij}
\left(x_j+r_j-r_i\right)
}{
m_i/h^2+\sum_j k_{ij}
}.
]

This is an important formula for your cookbook:

[
\boxed{
x_i^{\mathrm{PD}}
=================

\frac{
(m_i/h^2)y_i+
\sum_j k_{ij}(x_j+r_j-r_i)
}{
m_i/h^2+\sum_jk_{ij}
}
}
]

It shows explicitly how the inertial diagonal term connects the outer dynamics to the inner relaxation.

## Collision in PD

For each penetrating pair, project the pair separation to the collision surface:

[
p_{ij}
======

R_{ij}n.
]

The local projected relative displacement is (p_{ij}), and the quadratic energy is

[
E_{ij}
======

\frac{k_c}{2}
\left|
(x_j-x_i)-p_{ij}
\right|^2.
]

For nonpenetrating pairs, set

[
k_c=0.
]

Your current collision kernel performs essentially a Jacobi projection but omits the explicit inertial diagonal and physical stiffness. It computes

```cpp
dl = penetration/(invMi + invMj);
```

which is the infinite-stiffness PBD limit.

## PD pseudocode

```cpp
void hardStepPD(System& s, float h, int niter)
{
    // Inertial target
    for (Atom& a : s.atoms) {
        a.xPred = a.x + h*a.v;
        a.qPred = integrateQuat(a.q, a.omega, h);
    }

    for (int iter=0; iter<niter; iter++) {

        // Local ARAP projection of orientations
        for (Atom& a : coloredAtoms) {
            a.q = solveBestPortRotation(a, s);
        }

        // Rebuild active collision projections
        buildProjectedContacts(s);

        // Matrix-free global/Jacobi translation step
        for (Atom& a : s.atoms) {
            Vec3 rhs = (a.mass/(h*h))*a.xPred;
            float diag = a.mass/(h*h);

            for (PortBond& b : incidentBonds(a)) {
                rhs  += b.k * desiredCenterFromNeighbor(a,b);
                diag += b.k;
            }

            for (Contact& c : incidentContacts(a)) {
                rhs  += c.k * projectedCenterFromNeighbor(a,c);
                diag += c.k;
            }

            a.xNew = rhs/diag;
        }

        swapPositions();
    }

    for (Atom& a : s.atoms) {
        a.v     = (a.x-a.xOld)/h;
        a.omega = quatDifference(a.q,a.qOld)/h;
    }
}
```

## PD assessment

PD is attractive if:

* you already like your Kabsch/eigen rotational solve;
* topology is mostly fixed;
* you can tolerate several matrix-free Jacobi iterations.

Its weakness is that the translation step is still a global coupled relaxation. Without a good linear solver, it becomes similar to PBD but with a well-derived inertial diagonal.

---

# 5. Method B: XPBD

XPBD works directly with constraints and compliance. It was introduced specifically to make effective stiffness less dependent on timestep and iteration count. ([matthias-research.github.io][3])

For a vector port constraint

[
C=e=p_B-p_A=0,
]

use compliance

[
\alpha=\frac1{k}.
]

The timestep-scaled compliance is

[
\tilde\alpha=\frac{\alpha}{h^2}.
]

Let (\lambda\in\mathbb R^3) be the accumulated constraint multiplier. The vector XPBD update is

[
\boxed{
\Delta\lambda
=============

*

\left(
K+\tilde\alpha I
\right)^{-1}
\left(
e+\tilde\alpha\lambda
\right)
}
]

Then apply

[
x_A\mathrel{-}=w_A\Delta\lambda,
]

[
\theta_A\mathrel{-}=I_A^{-1}(r_A\times\Delta\lambda),
]

[
x_B\mathrel{+}=w_B\Delta\lambda,
]

[
\theta_B\mathrel{+}=I_B^{-1}(r_B\times\Delta\lambda),
]

and

[
\lambda\mathrel{+}=\Delta\lambda.
]

This should replace the scalar-normal treatment in your current kernel.

Your current code computes

```cpp
dist   = length(diff);
n      = diff/dist;
alpha  = 1/(K*dt*dt);
denom  = wi + wj + wang + alpha;
dlambda = dist/denom;
```

That solves only the scalar component along the instantaneous error direction. It is valid as a scalar distance constraint but loses tangential port-error coupling. For a point-to-point port attachment, the (3\times3) vector solve is better and costs only a small symmetric matrix inversion.

## Collision XPBD

For

[
C=r-R\ge0,
]

the update is scalar:

[
\Delta\lambda
=============

\frac{
-C-\tilde\alpha\lambda
}{
w_A+w_B+\tilde\alpha
}.
]

Projection enforces unilateral repulsion:

[
\lambda_{\mathrm{new}}
======================

\max(0,\lambda+\Delta\lambda),
]

[
\Delta\lambda
=============

\lambda_{\mathrm{new}}-\lambda.
]

Then

[
x_A\mathrel{-}=w_A n\Delta\lambda,
\qquad
x_B\mathrel{+}=w_B n\Delta\lambda.
]

If you do only one visit per microstep and do not need accurate complementarity, you can reset (\lambda=0) every microstep.

## XPBD pseudocode

```cpp
void solvePortXPBD(Atom& A, Atom& B, PortBond& b, float h)
{
    Vec3 rA = rotate(A.q, b.portA_local);
    Vec3 rB = rotate(B.q, b.portB_local);

    Vec3 e = (B.x+rB) - (A.x+rA);

    Mat3 K =
        identity()*(A.invMass+B.invMass)
      + skew(rA)*A.invInertiaWorld*transpose(skew(rA))
      + skew(rB)*B.invInertiaWorld*transpose(skew(rB));

    float alphaTilde = b.compliance/(h*h);

    Vec3 dlambda =
        -inverse(K + identity()*alphaTilde)
         *(e + b.lambda*alphaTilde);

    b.lambda += dlambda;

    A.x -= dlambda*A.invMass;
    B.x += dlambda*B.invMass;

    A.q = applyRotation(
        A.q,
        -A.invInertiaWorld*cross(rA,dlambda)
    );

    B.q = applyRotation(
        B.q,
         B.invInertiaWorld*cross(rB,dlambda)
    );
}
```

## XPBD integration loop

```cpp
void hardStepXPBD(System& s, float h, int iterations)
{
    for (Atom& a : s.atoms) {
        a.xOld = a.x;
        a.qOld = a.q;

        a.x += h*a.v;
        a.q  = integrateQuat(a.q,a.omega,h);
    }

    // Optional: clear multipliers every microstep
    clearLambdas(s);

    for (int iter=0; iter<iterations; iter++) {
        for (Color color : bondColors)
            parallelFor(bonds[color], solvePortXPBD);

        buildContacts(s);

        for (Color color : contactColors)
            parallelFor(contacts[color], solveContactXPBD);
    }

    for (Atom& a : s.atoms) {
        a.v     = (a.x-a.xOld)/h;
        a.omega = quatDifference(a.q,a.qOld)/h;
    }
}
```

## XPBD assessment

XPBD is the closest to your current `RRsp3.cl`.

Required changes:

1. use one symmetric bond record;
2. include both port arms;
3. use a (3\times3) effective mass;
4. optionally store `lambda`;
5. remove heuristic halving and guessed neighbor inertia.

---

# 6. Method C: sequential impulses with soft steps

SI works at the velocity level.

The same port constraint has:

[
e=p_B-p_A,
]

[
u=
v_B+\omega_B\times r_B
----------------------

v_A-\omega_A\times r_A.
]

For an infinitely rigid velocity constraint one would solve

[
P=-K^{-1}u.
]

But you want a soft, strongly damped local spring rather than an exact constraint.

Consider

[
f=-k e-c u.
]

Using backward Euler over the microstep gives a local implicit impulse

[
\boxed{
P=
--

\left[
I+(hc+h^2k)K
\right]^{-1}
\left[
hk,e+h(c+hk)u
\right]
}
]

This is the vector generalization of the scalar local implicit-spring formula.

It has precisely the connection you asked for:

* it uses no global linear solve;
* it contains the timestep (h) rigorously;
* it is an implicit discretization of the local harmonic hard energy;
* it does not need exact constraint satisfaction;
* arbitrarily large (k) does not create the explicit spring instability.

For very strong damping, you may simplify it to a positional-bias SI form:

[
P
=

*

K_{\mathrm{soft}}^{-1}
\left(
u+\beta\frac{e}{h}
\right),
]

where

[
K_{\mathrm{soft}}
=================

\gamma I+K.
]

Here:

* (\beta) controls the fraction of positional error removed per microstep;
* (\gamma) acts like compliance/softness.

Catto’s soft-step solver packages these terms into bias, effective-mass scaling and previous-impulse scaling. The central idea is still an implicit spring-damper impulse followed by position integration.

## Applying the impulse

```cpp
void applyPortImpulse(
    Atom& A, Atom& B,
    Vec3 rA, Vec3 rB,
    Vec3 P
){
    A.v     -= A.invMass*P;
    A.omega -= A.invInertiaWorld*cross(rA,P);

    B.v     += B.invMass*P;
    B.omega += B.invInertiaWorld*cross(rB,P);
}
```

## SI port solver

```cpp
void solvePortSI(
    Atom& A,
    Atom& B,
    const PortBond& b,
    float h
){
    Vec3 rA = rotate(A.q,b.portA_local);
    Vec3 rB = rotate(B.q,b.portB_local);

    Vec3 pA = A.x+rA;
    Vec3 pB = B.x+rB;
    Vec3 e  = pB-pA;

    Vec3 u =
        B.v + cross(B.omega,rB)
      - A.v - cross(A.omega,rA);

    Mat3 K =
        identity()*(A.invMass+B.invMass)
      + skew(rA)*A.invInertiaWorld*transpose(skew(rA))
      + skew(rB)*B.invInertiaWorld*transpose(skew(rB));

    float a = h*b.damping + h*h*b.stiffness;

    Mat3 lhs = identity() + K*a;

    Vec3 rhs =
          e*(h*b.stiffness)
        + u*(h*(b.damping+h*b.stiffness));

    Vec3 P = -inverse(lhs)*rhs;

    applyPortImpulse(A,B,rA,rB,P);
}
```

## SI collision

For collision, use the scalar version.

Let

[
C=r-R,
\qquad
u_n=n\cdot(v_B-v_A),
\qquad
w=w_A+w_B.
]

For (C<0), compute

[
p
=

*

\frac{
hkC+h(c+hk)u_n
}{
1+(hc+h^2k)w
}.
]

Then clamp:

[
p=\max(0,p).
]

```cpp
void solveContactSI(
    Atom& A,
    Atom& B,
    float R,
    float k,
    float c,
    float h
){
    Vec3 d = B.x-A.x;
    float r2 = dot(d,d);
    if (r2 >= R*R) return;

    float r = sqrt(r2);
    Vec3 n = d/r;

    float C  = r-R;
    float un = dot(n,B.v-A.v);
    float w  = A.invMass+B.invMass;

    float p =
        -(h*k*C + h*(c+h*k)*un)
        /(1.0f + (h*c+h*h*k)*w);

    p = max(0.0f,p);

    Vec3 P = n*p;

    A.v -= A.invMass*P;
    B.v += B.invMass*P;
}
```

## Soft-step microstep

```cpp
void hardStepSI(System& s, float h, int sweeps)
{
    for (int sweep=0; sweep<sweeps; sweep++) {

        for (Color color : bondColors)
            parallelFor(bonds[color], solvePortSI);

        buildContacts(s);

        for (Color color : contactColors)
            parallelFor(contacts[color], solveContactSI);
    }

    for (Atom& a : s.atoms) {
        a.x += h*a.v;
        a.q  = integrateQuat(a.q,a.omega,h);
    }
}
```

This is the cleanest answer to your question about replacing linear-solver iterations by hard-force substeps.

The physical connection is not the PD penalty

[
M/h^2
]

written explicitly in a global objective. Instead, that same inertia appears locally through the effective mass (K^{-1}), while (h^2k) determines how strongly the implicit spring competes with inertia.

For one scalar translational spring,

[
a=h^2kw.
]

That is exactly the ratio

[
\frac{\text{spring stiffness}}
{\text{inertial stiffness }m/h^2}.
]

So the SI formula contains the same balance as PD:

[
\boxed{
a
=

# \frac{k}{m_{\mathrm{eff}}/h^2}

h^2k,m_{\mathrm{eff}}^{-1}
}
]

It is merely solved locally in impulse form rather than globally in position form.

---

# 7. Method D: Vertex Block Descent

VBD minimizes the same implicit Euler objective as PD:

[
\Phi(q)
=======

\frac{1}{2h^2}(q-y)^TM(q-y)
+
E_h(q),
]

but instead of alternating constraint projections and a global solve, it minimizes (\Phi) with respect to one vertex/body block at a time.

VBD is therefore a block-coordinate Gauss–Seidel solver for implicit Euler. It was introduced as an unconditionally stable local-block method that can be stopped after a fixed iteration budget. ([arXiv][4])

For your oriented atom, one block should be six-dimensional:

[
z_i=
\begin{pmatrix}
x_i\
\theta_i
\end{pmatrix}.
]

Hold all neighbors fixed, then solve

[
H_i\Delta z_i=-g_i,
]

where (g_i) and (H_i) contain:

* inertial translation;
* inertial rotation;
* all incident ports;
* active collisions.

## Inertial contribution

[
g_{x,\mathrm{inertia}}
======================

\frac{m_i}{h^2}(x_i-y_i),
]

[
H_{xx,\mathrm{inertia}}
=======================

\frac{m_i}{h^2}I.
]

For rotation-vector increments,

[
g_{\theta,\mathrm{inertia}}
\approx
\frac{I_i}{h^2}\theta_{\mathrm{pred-error}},
]

[
H_{\theta\theta,\mathrm{inertia}}
\approx
\frac{I_i}{h^2}.
]

## Port contribution

For one incident bond with neighbor fixed,

[
e=p_j-p_i.
]

The Jacobian with respect to the local six-dimensional block is

[
J_i=
\begin{bmatrix}
-I & S(r_i)
\end{bmatrix}
]

up to the chosen sign convention.

Then

[
g_i
===

kJ_i^Te,
]

[
H_i
\approx
kJ_i^TJ_i.
]

This is a Gauss–Newton approximation.

Your existing `compute_ports_cluster_rigid_substep_optimized()` already computes the rotational (3\times3) sub-block

[
\sum_j k_j
\left(
|r_j|^2I-r_jr_j^T
\right)
]

and solves it by Newton iteration.

VBD simply extends that idea to the complete (6\times6) body block, including translation-rotation coupling.

## VBD pseudocode

```cpp
void updateBodyVBD(System& s, int i, float h)
{
    Atom& A = s.atoms[i];

    Vec6 g = zero6();
    Mat6 H = zero6x6();

    // Inertial regularization
    g.x = A.mass*(A.x-A.xPred)/(h*h);
    H.xx += identity()*A.mass/(h*h);

    Vec3 rotError = rotationDifference(A.q,A.qPred);
    g.r = A.inertiaWorld*rotError/(h*h);
    H.rr += A.inertiaWorld/(h*h);

    for (PortBond& b : incidentBonds(i)) {
        PortJacobian J;
        Vec3 e = evaluatePortErrorAndLocalJacobian(s,b,i,J);

        g += transpose(J)*(b.stiffness*e);
        H += transpose(J)*(b.stiffness*J);
    }

    for (Contact& c : incidentContacts(i)) {
        if (c.C >= 0.0f) continue;

        Vec6 J = contactJacobianForBody(c,i);

        g += J*(c.stiffness*c.C);
        H += outer(J,J)*c.stiffness;
    }

    H += diagonalRegularization();

    Vec6 dz = solve6x6(H,-g);

    A.x += dz.x;
    A.q  = applyRotation(A.q,dz.r);
}
```

Run atom colors sequentially and atoms within one color in parallel:

```cpp
for (int iter=0; iter<niter; iter++) {
    for (Color color : atomColors) {
        parallelFor(atoms[color], updateBodyVBD);
    }
}
```

## Why VBD may fit your ports extremely well

Compared with XPBD or SI:

* XPBD solves one bond constraint at a time;
* SI solves one bond impulse at a time;
* VBD solves all four ports of one atom together.

Thus VBD directly computes the best combined translational and rotational correction for the central atom.

That is almost exactly what your massless Newton/Kabsch kernels already try to do, except they currently:

1. solve rotation separately;
2. subsequently compute linear recoil;
3. do not include the inertial (6\times6) block consistently.

A proper rigid-body VBD kernel would unify those operations.

---

# 8. Augmented VBD

AVBD adds augmented-Lagrangian multipliers to VBD. It was developed to improve VBD for very high stiffness ratios and exact hard constraints. ([graphics.cs.utah.edu][5])

For a constraint (C(q)=0), augment the energy by

[
E_{\mathrm{aug}}
================

\lambda^TC(q)
+
\frac{\rho}{2}|C(q)|^2.
]

The body-local gradient uses

[
J^T(\lambda+\rho C),
]

and the approximate Hessian uses

[
\rho J^TJ.
]

After the body sweep, update

[
\lambda\leftarrow\lambda+\rho C.
]

For collisions, project the multiplier:

[
\lambda_n
\leftarrow
\max(0,\lambda_n+\rho C).
]

This allows very stiff or formally exact constraints without making the local Hessian arbitrarily ill-conditioned.

However, for your stated goal, plain VBD is probably enough. AVBD introduces multiplier memory and extra bandwidth, just like warm-started XPBD/SI.

Use AVBD only if you find that:

* bonds gradually stretch too much;
* a large stiffness ratio prevents VBD convergence;
* increasing the local penalty makes the (6\times6) solve poorly conditioned.

---

# 9. Unified comparison

All four methods use the same model:

[
e_{ij}
======

## x_j+R_ja_{ji}

x_i-R_ia_{ij}.
]

They differ only in what variable is corrected.

| Method       | Unknown corrected locally     | Inertia representation                   | Stored multiplier |
| ------------ | ----------------------------- | ---------------------------------------- | ----------------- |
| PD           | positions/orientations        | global (M/h^2) diagonal                  | no                |
| XPBD         | positions/orientations        | effective inverse mass plus (\alpha/h^2) | optional/yes      |
| SI soft step | velocities/angular velocities | effective mass in impulse solve          | optional          |
| VBD          | one body’s (x,\theta) block   | local (6\times6) (M/h^2) Hessian         | no                |
| AVBD         | one body block + multipliers  | local Hessian + augmented penalty        | yes               |

The formulas can be viewed as different approximations to the same implicit objective.

## PD

[
\left(
\frac{M}{h^2}+H
\right)\Delta q=-g.
]

Solve globally or approximately.

## VBD

Take one body block of the same equation:

[
\left(
\frac{M_i}{h^2}+H_{ii}
\right)\Delta q_i=-g_i.
]

## XPBD

Eliminate generalized-coordinate correction in favor of constraint multiplier:

[
\Delta\lambda
=============

*

\left(
JM^{-1}J^T+\alpha/h^2
\right)^{-1}
\left(C+\alpha\lambda/h^2\right).
]

## Soft SI

Discretize the spring-damper in velocity space:

[
P
=

*

\left[
I+(hc+h^2k)JM^{-1}J^T
\right]^{-1}
\left[
hkC+h(c+hk)\dot C
\right].
]

So yes: the inertial coupling that you value in PD is fully present in local soft SI. It appears through

[
JM^{-1}J^T
]

and the dimensionless combination

[
h^2kJM^{-1}J^T.
]

---

# 10. Recommended cookbook architecture

I would keep one common outer driver and swap only the hard solver.

```cpp
void macroStep(System& s, float H, int nsub, SolverKind method)
{
    // Cached Fsoft(x_n)
    applySoftKick(s, 0.5f*H);

    float h = H/nsub;

    for (int isub=0; isub<nsub; isub++) {

        updateBroadphaseIfNeeded(s);
        buildActiveContacts(s);

        switch(method) {
            case PD:
                hardStepPD(s,h,2);
                break;

            case XPBD:
                hardStepXPBD(s,h,1);
                break;

            case SOFT_SI:
                hardStepSI(s,h,1);
                break;

            case VBD:
                hardStepVBD(s,h,1);
                break;
        }
    }

    evaluateSoftForces(s);     // Fsoft(x_{n+1})
    applySoftKick(s,0.5f*H);

    applyRelaxationOptimizer(s);   // damping, FIRE, etc.
}
```

For an initial comparison, use the same budget:

* one expensive soft evaluation per macrostep;
* eight hard microsteps;
* one hard sweep per microstep;
* no warm starting;
* identical contact list and port topology.

Measure:

1. wall time per macrostep;
2. maximum stable (H);
3. average port error;
4. maximum penetration;
5. energy reached after a fixed wall time;
6. number of expensive soft evaluations to convergence.

That last quantity is probably the most meaningful metric for your application.

---

# 11. What I would implement first

My expected ranking for your particular rigid-port molecular problem is:

### 1. Rigid-body VBD

Likely best convergence per memory access because it solves all four ports and rotation of one atom together. It is also closest to your existing Newton/Kabsch kernels.

Implement a (6\times6) local Gauss–Newton solve with inertial diagonal:

[
H_i
===

\frac{M_i}{h^2}
+
\sum_{b\ni i}k_bJ_b^TJ_b
+
\sum_{c\ni i}k_cJ_c^TJ_c.
]

### 2. Soft-step SI

Likely simplest and cheapest. It preserves momentum cleanly and needs no multiplier memory. It may require more microsteps because each bond is processed separately.

### 3. XPBD

Very close to your existing implementation and useful as a reference. The multiplier can be reset each microstep initially.

### 4. Matrix-free PD

Useful because it gives the clearest variational reference, but without a global solver its advantage over VBD/XPBD may be limited.

---

# 12. Minimal common C++ skeleton

```cpp
struct Atom {
    Vec3 x, v;
    Quat q;
    Vec3 omega;

    float invMass;
    Mat3 invInertia;
};

struct Bond {
    int ia, ib;
    Vec3 portA;
    Vec3 portB;

    float k;
    float damping;
    float compliance;

    Vec3 lambda;
};

struct PortState {
    Vec3 rA, rB;
    Vec3 error;
    Vec3 relVelocity;
    Mat3 invEffectiveMass;
};

PortState evaluatePort(
    const Atom& A,
    const Atom& B,
    const Bond& b
){
    PortState s;

    s.rA = rotate(A.q,b.portA);
    s.rB = rotate(B.q,b.portB);

    s.error =
        (B.x+s.rB)
      - (A.x+s.rA);

    s.relVelocity =
        B.v + cross(B.omega,s.rB)
      - A.v - cross(A.omega,s.rA);

    Mat3 K =
        identity()*(A.invMass+B.invMass)
      + skew(s.rA)*A.invInertia*transpose(skew(s.rA))
      + skew(s.rB)*B.invInertia*transpose(skew(s.rB));

    s.invEffectiveMass = inverse(K);

    return s;
}

void applyPositionCorrection(
    Atom& A,
    Atom& B,
    const PortState& s,
    Vec3 dlambda
){
    A.x -= A.invMass*dlambda;
    B.x += B.invMass*dlambda;

    A.q = applyRotation(
        A.q,
        -A.invInertia*cross(s.rA,dlambda)
    );

    B.q = applyRotation(
        B.q,
         B.invInertia*cross(s.rB,dlambda)
    );
}

void applyVelocityImpulse(
    Atom& A,
    Atom& B,
    const PortState& s,
    Vec3 impulse
){
    A.v -= A.invMass*impulse;
    B.v += B.invMass*impulse;

    A.omega -= A.invInertia*cross(s.rA,impulse);
    B.omega += B.invInertia*cross(s.rB,impulse);
}
```

Then the three simplest variants become:

```cpp
void solveXPBD(Atom& A, Atom& B, Bond& b, float h)
{
    PortState s = evaluatePort(A,B,b);

    float at = b.compliance/(h*h);

    Mat3 W = inverse(s.invEffectiveMass) + identity()*at;

    Vec3 dlambda =
        -inverse(W)*(s.error + b.lambda*at);

    b.lambda += dlambda;

    applyPositionCorrection(A,B,s,dlambda);
}
```

```cpp
void solveSoftSI(Atom& A, Atom& B, Bond& b, float h)
{
    PortState s = evaluatePort(A,B,b);

    Mat3 K = inverse(s.invEffectiveMass);

    float a = h*b.damping + h*h*b.k;

    Vec3 rhs =
          s.error*(h*b.k)
        + s.relVelocity*(h*(b.damping+h*b.k));

    Vec3 impulse =
        -inverse(identity()+K*a)*rhs;

    applyVelocityImpulse(A,B,s,impulse);
}
```

```cpp
void solveInfinitePBD(Atom& A, Atom& B, Bond& b)
{
    PortState s = evaluatePort(A,B,b);

    Vec3 dlambda =
        -s.invEffectiveMass*s.error;

    applyPositionCorrection(A,B,s,dlambda);
}
```

PD and VBD use the same error and Jacobians but accumulate them into position-system or body-block normal equations.

---

The main conceptual simplification is:

[
\boxed{
\text{one symmetric port-to-port residual}
+
\text{one common rigid-body Jacobian}
}
]

Everything else—PD, XPBD, SI, VBD—is just a different strategy for deciding and applying the correction generated by that residual.

[Uploaded `RRsp3.cl`](sandbox:/mnt/data/RRsp3.cl)

[1]: https://igl.ethz.ch/projects/ARAP/arap_web.pdf?utm_source=chatgpt.com "As-Rigid-As-Possible Surface Modeling"
[2]: https://users.cs.utah.edu/~ladislav/bouaziz14projective/bouaziz14projective.pdf?utm_source=chatgpt.com "Fusing Constraint Projections for Fast Simulation"
[3]: https://matthias-research.github.io/pages/publications/XPBD.pdf?utm_source=chatgpt.com "XPBD: Position-Based Simulation of Compliant ..."
[4]: https://arxiv.org/abs/2403.06321?utm_source=chatgpt.com "Vertex Block Descent"
[5]: https://graphics.cs.utah.edu/research/projects/avbd/Augmented_VBD-SIGGRAPH25.pdf?utm_source=chatgpt.com "Augmented Vertex Block Descent"

---

# USER

OK, I think you did  not get the representation very well.

I'm my represetnation I have everything atom-centered. I do not solve bonds, I solve atoms position and rotation - this should be common for all methods SI,PD,XPBD,(A)VBD
my dynamical variables which are store in global mmeory are only position and rotation of atom (pos+normalized quaternion) and associated angular/linar vecloities (forces and torques are computed on the fly no need to store them). 

The other variables (like position of ports) are also computed n the fly, they do not need to be stored. For efficiency it may be usefull store it in some global axuliary buffer but that is just optimization not formulation. 

In you data type I see that bond contains  Bond{ Vec3 portA; Vec3 portB; }, I do not want that, bond should be only indexes of atoms/port and eventually stiffness maybe compliace or lambda if you want to store the lagrange multipliers between the steps.

you say "Why VBD may fit your ports extremely well ... XPBD solves one bond constraint at a time ... SI solves one bond impulse at a time ... VBD solves all four ports of one atom together."

this is in contrast to my view - in my view we always gather all  forces and torques from all ports (contstrains) and collisions on given atom, and only then we move the atom. Look again carefully how our kernels in RRsp3.cl works, I think you missed that aspect completely. Respect it and fomulate you distinction between the metods XPBD,PD,SI,VBD around this general atom-centered formulation.


---

You are correct. I misidentified your representation as an edge-centered constraint system. The uploaded kernel is very clearly an **atom-centered Jacobi gather–apply architecture**:

* `compute_ports_cluster_rigid()` assigns one work-item to a central atom, loops over its four ports, and accumulates one `sum_dpos` and one `sum_dtheta`.
* The reaction acting on each neighboring atom is written to a temporary `dpos_neigh[port_slot]`.
* `apply_corrections_rigid_ports()` gathers:

  * the atom’s own central-port correction,
  * all incoming recoils through `bkSlots`,
  * its accumulated collision correction,
  * and only then updates the atom once.
* `compute_collision_cluster_rigid()` similarly gathers all overlaps acting on one atom before writing one `dpos_coll[i]`.
* `compute_ports_cluster_rigid_substep_optimized()` already gathers all four ports into one torque and one (3\times3) rotational Hessian before updating the quaternion.

So the invariant architecture should indeed be:

[
\boxed{
\text{one atom thread}
\rightarrow
\text{gather all incident port/collision effects}
\rightarrow
\text{compute one translational/rotational update}
}
]

The distinction between PD, XPBD, SI, and VBD is **what is accumulated and what local equation produces that one atom update**, not whether constraints are processed one bond at a time.

---

# 1. Correct atom-centered representation

For each atom (i), the dynamic global state is only

[
x_i,\quad q_i,\quad v_i,\quad \omega_i.
]

Static topology and geometry are something like:

```cpp
struct AtomTopology {
    int   neigh[4];       // neigh[i][k] = atom wanted at port k
    Vec3  portLocal[4];   // fixed local port vectors
    float stiffness[4];

    // Optional only:
    float lambda[4];      // XPBD/AVBD or warm starting
};
```

There is no object containing `portA` and `portB`.

A directed port is simply

[
(i,k)\longrightarrow j.
]

Its world-space arm and tip are

[
r_{ik}=R(q_i)a_{ik},
]

[
t_{ik}=x_i+r_{ik}.
]

The atom-centered residual is

[
\boxed{
e_{ik}=t_{ik}-x_j
=================

x_i+R(q_i)a_{ik}-x_j.
}
]

The energy represented by all ports is

[
\boxed{
E_{\rm port}
============

\frac12\sum_i\sum_{k=1}^{n_i}
K_{ik}
\left|
x_i+R(q_i)a_{ik}-x_{j(i,k)}
\right|^2.
}
]

This is exactly the ARAP-like interpretation:

[
R_i a_{ik}
\approx
x_j-x_i.
]

Each atom has its own local rotation (R_i) that tries to align its local port pattern with the current neighbor vectors.

If atom (j) also has a port pointing toward (i), there is another directed residual

[
e_{jl}
======

x_j+R_j a_{jl}-x_i.
]

These are **two atom-centered terms**, not one bond with two stored ports.

---

# 2. One directed port: common local quantities

Define

[
e=x_i+r-x_j,
\qquad
C=|e|,
\qquad
n=\frac{e}{|e|}.
]

An infinitesimal atom rotation (\delta\theta_i) changes the port arm by

[
\delta r
========

\delta\theta_i\times r.
]

Therefore

[
\delta C
========

n\cdot\delta x_i
+
(r\times n)\cdot\delta\theta_i
------------------------------

n\cdot\delta x_j.
]

Define

[
s=r\times n.
]

The scalar constraint Jacobians are then

[
J_i=
\begin{bmatrix}
n^T & s^T
\end{bmatrix},
\qquad
J_j=-n^T.
]

The effective inverse mass is

[
\boxed{
w_{\rm eff}
===========

w_i+w_j+s^TI_i^{-1}s,
}
]

where

[
w_i=m_i^{-1}.
]

Crucially, there is **no angular term for atom (j)** in this directed port residual. Its orientation does not occur in

[
x_i+R_i a_{ik}-x_j.
]

The rotation of (j) is controlled by its own outgoing ports.

This means the current heuristic in `compute_ports_cluster_rigid()` that:

1. adds another copy of `w_ang` when the neighbor is a node, and
2. then halves the impulse,

does not follow naturally from the atom-centered directed-port energy.

If both directed port terms exist, their stiffnesses should be chosen consistently. For example, use (K/2) for each directed term if their combination should correspond approximately to one bond stiffness (K). That should be a topology/model decision, not a runtime `j_isnode` correction.

---

# 3. Common gather–apply pipeline

All methods can use the same architecture:

```cpp
for each solver iteration or microstep:

    // Evaluation/gather stage
    parallel_for atom i:
        clear atom accumulator

        for each outgoing port (i,k)->j:
            evaluate current port geometry
            compute method-specific contribution

            accumulate central contribution for i
            write/scatter recoil contribution for j

        gather collision contributions acting on i

    // Apply stage
    parallel_for atom i:
        gather:
            own outgoing-port contribution
            incoming port recoils
            collisions

        update atom i exactly once
```

Your existing buffers map naturally:

| Current buffer     | Meaning                                             |
| ------------------ | --------------------------------------------------- |
| `dpos_node[i]`     | sum of effects from atom (i)’s outgoing ports       |
| `drot_node[i]`     | sum/solve of rotational effects from outgoing ports |
| `dpos_neigh[slot]` | recoil sent to the target neighbor                  |
| `bkSlots[i]`       | incoming recoil slots for atom (i)                  |
| `dpos_coll[i]`     | all collision effects gathered on atom (i)          |

These are transient solver buffers, not additional molecular degrees of freedom.

---

# 4. XPBD in the atom-centered architecture

This is closest to the current `compute_ports_cluster_rigid()`.

For one directed port, XPBD gives

[
\tilde\alpha
============

# \frac{\alpha}{h^2}

\frac{1}{Kh^2},
]

[
\boxed{
\Delta\lambda_{ik}
==================

\frac{
-C_{ik}-\tilde\alpha\lambda_{ik}
}{
w_{\rm eff}+\tilde\alpha
}.
}
]

If you do not store multipliers, set

[
\lambda_{ik}=0
]

at each iteration or microstep:

[
\Delta\lambda_{ik}
==================

-\frac{C_{ik}}
{w_{\rm eff}+\tilde\alpha}.
]

The contributions are

[
\Delta x_i
==========

w_i n,\Delta\lambda,
]

[
\Delta\theta_i
==============

I_i^{-1}s,\Delta\lambda,
]

[
\Delta x_j
==========

-w_j n,\Delta\lambda.
]

Your atom thread gathers

[
\Delta x_i^{\rm own}
====================

\sum_{k\in{\rm out}(i)}
w_i n_{ik}\Delta\lambda_{ik},
]

[
\Delta\theta_i
==============

\sum_{k\in{\rm out}(i)}
I_i^{-1}s_{ik}\Delta\lambda_{ik}.
]

Each outgoing port also writes the target recoil

[
\Delta x_j^{(ik)}
=================

-w_jn_{ik}\Delta\lambda_{ik}.
]

The apply kernel later gathers all such incoming terms for (j).

## Atom-centered XPBD pseudocode

```cpp
void computeAtomXPBD(int i, float h)
{
    Vec3 dxOwn = {0,0,0};
    Vec3 dtheta = {0,0,0};

    for (int k=0; k<nport[i]; k++) {
        int j = neigh[i][k];
        if (j < 0) continue;

        Vec3 r = rotate(q[i], portLocal[i][k]);
        Vec3 e = x[i] + r - x[j];

        float C = length(e);
        if (C < 1e-8f) continue;

        Vec3 n = e/C;
        Vec3 s = cross(r,n);

        float weff =
              invMass[i]
            + invMass[j]
            + dot(s, invInertia[i]*s);

        float alphaTilde =
            compliance[i][k]/(h*h);

        float dlambda =
            (-C - alphaTilde*lambda[i][k])
            /(weff + alphaTilde);

        // Optional:
        lambda[i][k] += dlambda;

        dxOwn += invMass[i]*n*dlambda;
        dtheta += invInertia[i]*s*dlambda;

        recoil[i][k] =
            -invMass[j]*n*dlambda;
    }

    dxNode[i] = dxOwn;
    dthetaNode[i] = dtheta;
}
```

Then:

```cpp
void applyAtomXPBD(int i)
{
    Vec3 dx =
          dxNode[i]
        + collisionCorrection[i];

    for (int s=0; s<nIncoming[i]; s++) {
        dx += recoil[incomingSlot[i][s]];
    }

    x[i] += relaxation*dx;
    q[i] = applyRotation(
        q[i],
        relaxation*dthetaNode[i]
    );
}
```

This is XPBD-Jacobi: every atom gathers all effects, and all atoms are updated simultaneously.

---

# 5. Soft-step impulses in the same atom-centered architecture

For the directed port,

[
C=|x_i+r-x_j|,
]

and its velocity is

[
\boxed{
\dot C
======

n\cdot
\left(
v_i+\omega_i\times r-v_j
\right).
}
]

Again, (\omega_j) does not occur because the target of this directed port is the center of atom (j).

For a local harmonic spring-damper,

[
f=-KC-D\dot C,
]

a backward-Euler local solve gives the scalar impulse

[
\boxed{
p
=

*

\frac{
hKC+h(D+hK)\dot C
}{
1+(hD+h^2K)w_{\rm eff}
}.
}
]

Then gather velocity changes:

[
\Delta v_i
==========

w_i n p,
]

[
\Delta\omega_i
==============

I_i^{-1}s,p,
]

[
\Delta v_j
==========

-w_jnp.
]

## Atom-centered soft-impulse pseudocode

```cpp
void computeAtomSoftImpulse(int i, float h)
{
    Vec3 dvOwn = {0,0,0};
    Vec3 domega = {0,0,0};

    for (int k=0; k<nport[i]; k++) {
        int j = neigh[i][k];
        if (j < 0) continue;

        Vec3 r = rotate(q[i], portLocal[i][k]);
        Vec3 e = x[i] + r - x[j];

        float C = length(e);
        if (C < 1e-8f) continue;

        Vec3 n = e/C;
        Vec3 s = cross(r,n);

        float Cdot = dot(
            n,
            v[i] + cross(omega[i],r) - v[j]
        );

        float weff =
              invMass[i]
            + invMass[j]
            + dot(s,invInertia[i]*s);

        float K = stiffness[i][k];
        float D = damping[i][k];

        float impulse =
            -(h*K*C + h*(D+h*K)*Cdot)
            /(1.0f + (h*D+h*h*K)*weff);

        dvOwn += invMass[i]*n*impulse;
        domega += invInertia[i]*s*impulse;

        velocityRecoil[i][k] =
            -invMass[j]*n*impulse;
    }

    dvNode[i] = dvOwn;
    domegaNode[i] = domega;
}
```

Apply once:

```cpp
void applyAtomSoftImpulse(int i, float h)
{
    Vec3 dv =
          dvNode[i]
        + dvCollision[i];

    for (int s=0; s<nIncoming[i]; s++) {
        dv += velocityRecoil[incomingSlot[i][s]];
    }

    v[i] += dv;
    omega[i] += domegaNode[i];

    x[i] += h*v[i];
    q[i] = integrateQuaternion(q[i],omega[i],h);
}
```

## Important naming correction

If every impulse is computed from the same old (x,q,v,\omega), gathered, and then applied simultaneously, this is mathematically a:

[
\boxed{\text{Jacobi soft-impulse solver}}
]

not strictly “sequential impulses.”

Sequential impulses means later constraints see velocities changed by earlier constraints. Catto’s Soft Step combines soft constraints with substepping, but its constraint solve is still Gauss–Seidel-like within the available ordering/coloring. ([box2d.org][1])

You can retain your atom-centered gather while introducing some sequentiality by processing atom colors:

```cpp
for color in colors:
    parallel_for atom i in color:
        gather all port/collision impulses on i
        update i
```

That is an **atom-block Gauss–Seidel soft-impulse method**. But your current two-kernel gather/apply design is the Jacobi version, which may be preferable on the GPU.

---

# 6. Projective Dynamics in the atom-centered representation

PD minimizes an implicit objective

[
\boxed{
\Phi(x,q)
=========

\sum_i
\frac{m_i}{2h^2}|x_i-y_i|^2
+
\sum_i
\frac{1}{2h^2}
\left|
\operatorname{Log}(q_iq_{i,\rm pred}^{-1})
\right|*{I_i}^2
+
E*{\rm port}(x,q)
+
E_{\rm coll}(x).
}
]

Here

[
y_i=x_i^n+h v_i^n
]

is the inertial target.

In your atom-centered ARAP formulation, the natural PD splitting is:

1. **Per-atom local rotation projection.**
2. **Global translation compromise**, approximated by atom-centered Jacobi updates rather than Cholesky.

## PD local rotation step

For each atom (i), gather all its outgoing ports and solve

[
\boxed{
R_i^\star
=========

\arg\min_{R\in SO(3)}
\sum_k
K_{ik}
\left|
R a_{ik}-(x_j-x_i)
\right|^2.
}
]

This is exactly what your:

* shape-matching kernel,
* Davenport eigen kernel,
* Newton rotational kernel,

are doing.

`compute_ports_cluster_rigid_shapematch()` is therefore naturally interpreted as an ARAP/PD local projection.

Your Newton kernel gathers

[
\tau_i
======

\sum_k
K_{ik},r_{ik}\times
\left[x_j-(x_i+r_{ik})\right],
]

and

[
H_i^\theta
\approx
\sum_k K_{ik}
\left(
|r_{ik}|^2I-r_{ik}r_{ik}^T
\right).
]

That is already an atom-centered rotational block solve.

For dynamic PD, add rotational inertia:

[
H_i^\theta
\leftarrow
H_i^\theta+\frac{I_i}{h^2},
]

and add the predicted-orientation error to the rotational gradient.

## PD translation step

With rotations fixed, the directed port term is quadratic:

[
E_{ik}
======

\frac{K_{ik}}2
\left|
x_i+r_{ik}-x_j
\right|^2.
]

For atom (i), gather:

* outgoing targets (x_j-r_{ik});
* incoming targets (x_l+r_{lm});
* projected collision targets;
* inertial target (y_i).

The atom-centered Jacobi update is

[
\boxed{
x_i^{\rm new}
=============

\frac{
\displaystyle
\frac{m_i}{h^2}y_i
+
\sum_{(i,k)\to j}
K_{ik}(x_j-r_{ik})
+
\sum_{(l,m)\to i}
K_{lm}(x_l+r_{lm})
+
b_i^{\rm coll}
}{
\displaystyle
\frac{m_i}{h^2}
+
\sum_{{\rm outgoing}}K
+
\sum_{{\rm incoming}}K
+
D_i^{\rm coll}
}.
}
]

This is the proper PD analogue of your current gather/apply kernel.

The distinguishing feature is:

> PD gathers **projected target positions and diagonal coefficients**, then solves a local Jacobi approximation to the global quadratic system.

It does not gather per-port XPBD corrections or velocity impulses.

---

# 7. VBD in the atom-centered representation

VBD is even closer to your intended view.

VBD minimizes the implicit-Euler objective by gathering all energies incident on one vertex and solving one local block. The original method is explicitly formulated as vertex-level block coordinate descent with Gauss–Seidel iterations. ([arXiv][2])

For your rigid atom, the block is six-dimensional:

[
z_i=
\begin{pmatrix}
x_i\
\theta_i
\end{pmatrix}.
]

For one outgoing vector residual

[
e_{ik}=x_i+r_{ik}-x_j,
]

its Jacobian with respect to the central atom block is

[
\boxed{
J_{ik,i}
========

\begin{bmatrix}
I & -[r_{ik}]_\times
\end{bmatrix}.
}
]

An incoming directed port contributes only to translation of atom (i), because (i) appears only as the neighbor center:

[
J_{{\rm incoming},i}
====================

\begin{bmatrix}
-I & 0
\end{bmatrix}.
]

For each atom, gather the six-dimensional gradient

[
g_i
===

\frac{M_i}{h^2}(z_i-z_{i,\rm pred})
+
\sum_{\rm incident}
K_aJ_{a,i}^Te_a
+
g_i^{\rm coll},
]

and the local Gauss–Newton Hessian

[
\boxed{
H_i
===

\frac{M_i}{h^2}
+
\sum_{\rm incident}
K_aJ_{a,i}^TJ_{a,i}
+
H_i^{\rm coll}.
}
]

Then solve

[
\boxed{
H_i\Delta z_i=-g_i.
}
]

Only after gathering every outgoing port, every incoming port recoil, and every collision do you update

[
x_i\leftarrow x_i+\Delta x_i,
]

[
q_i\leftarrow
\exp(\Delta\theta_i),q_i.
]

## Atom-centered VBD pseudocode

```cpp
void computeAtomVBD(int i, float h)
{
    Vec6 g = {};
    Mat6 H = {};

    // Inertial diagonal regularization
    g.translation =
        mass[i]*(x[i]-xPred[i])/(h*h);

    g.rotation =
        inertia[i]*rotationError(q[i],qPred[i])/(h*h);

    H.tt += identity3()*mass[i]/(h*h);
    H.rr += inertia[i]/(h*h);

    // Outgoing ports: translation + rotation
    for (int k=0; k<nport[i]; k++) {
        int j = neigh[i][k];

        Vec3 r = rotate(q[i],portLocal[i][k]);
        Vec3 e = x[i]+r-x[j];

        Mat3 S = skew(r);
        Mat3x6 J = concatenate(identity3(),-S);

        float K = stiffness[i][k];

        g += transpose(J)*(K*e);
        H += transpose(J)*(K*J);
    }

    // Incoming ports: translation only
    for (int s=0; s<nIncoming[i]; s++) {
        PortRef p = incomingPort[i][s];

        int l = p.sourceAtom;
        int k = p.sourcePort;

        Vec3 r = rotate(q[l],portLocal[l][k]);
        Vec3 e = x[l]+r-x[i];

        // d e / d x_i = -I
        g.translation -= stiffness[l][k]*e;
        H.tt += identity3()*stiffness[l][k];
    }

    addCollisionGradientHessian(i,g,H);

    Vec6 dz = solve6x6(H,-g);

    dzAtom[i] = dz;
}
```

Then apply all atoms simultaneously for Jacobi VBD, or color atoms and apply after each color for standard Gauss–Seidel VBD.

## Relation to your Newton rotational kernel

Your `compute_ports_cluster_rigid_substep_optimized()` is already approximately:

[
H_{\theta\theta}\Delta\theta=-g_\theta.
]

So it is essentially the rotational (3\times3) sub-block of atom-centered VBD, but currently:

* translation is handled separately;
* translation–rotation coupling is omitted;
* the rotational inertial term (I/h^2) is omitted;
* incoming ports do not affect rotation, correctly;
* it is operating in a massless/geometric mode.

A full VBD kernel would extend this from (3\times3) rotation to one (6\times6) atom block.

---

# 8. AVBD in the same architecture

AVBD keeps the same atom-centered VBD gather, but augments each port constraint with a multiplier.

For a residual (e_a), use

[
E_a^{\rm aug}
=============

\lambda_a^Te_a
+
\frac{\rho_a}{2}|e_a|^2.
]

Its contribution to the atom block is

[
g_i
\mathrel{+}=
J_{a,i}^T
\left(
\lambda_a+\rho_ae_a
\right),
]

[
H_i
\mathrel{+}=
\rho_aJ_{a,i}^TJ_{a,i}.
]

After a sweep,

[
\lambda_a
\leftarrow
\lambda_a+\rho_ae_a.
]

For collisions, project the normal multiplier so that it cannot become attractive.

AVBD was introduced to improve VBD for hard constraints and large stiffness ratios using an augmented-Lagrangian formulation. ([graphics.cs.utah.edu][3])

But it does require persistent per-port or per-contact multiplier memory. Thus:

* plain VBD respects your minimal dynamic-state preference;
* AVBD pays extra bandwidth for better high-stiffness convergence;
* if multipliers are reset every time, the method falls back toward penalty VBD.

---

# 9. Correct distinction between the four methods

All four can use precisely your atom-centered gather–apply representation.

| Method       | What each atom gathers                                                       | Local atom result                                                            | When atom moves                             |
| ------------ | ---------------------------------------------------------------------------- | ---------------------------------------------------------------------------- | ------------------------------------------- |
| PD           | projected port targets, collision projections, inertial RHS and diagonal     | Jacobi solution of quadratic global step; separate local rotation projection | once after gathering                        |
| XPBD         | per-port compliant position corrections (\Delta\lambda) and neighbor recoils | (\Delta x_i,\Delta\theta_i)                                                  | once after gathering                        |
| Soft impulse | per-port velocity impulses based on (C,\dot C), plus collision impulses      | (\Delta v_i,\Delta\omega_i)                                                  | velocity updated once, then one integration |
| VBD          | full gradient and (6\times6) Hessian from every incident energy              | (\Delta x_i,\Delta\theta_i=-H_i^{-1}g_i)                                     | once after gathering                        |
| AVBD         | VBD gradient/Hessian plus augmented multipliers                              | same (6\times6) block solve                                                  | once after gathering                        |

The essential conceptual differences are:

### PD

[
\text{project each local energy}
\rightarrow
\text{solve global quadratic approximately}.
]

### XPBD

[
\text{solve each constraint multiplier approximately}
\rightarrow
\text{sum all coordinate corrections}.
]

### Soft impulses

[
\text{solve each local spring in velocity/impulse form}
\rightarrow
\text{sum all velocity changes}.
]

### VBD

[
\text{sum all energies into one atom Hessian}
\rightarrow
\text{solve the atom block directly}.
]

---

# 10. What should change in the current kernel

The corrected interpretation suggests these concrete changes.

## Keep

* `pos`, `quat`, `vel`, `omega` as the only dynamic atom state.
* `port_local[i*4+k]`.
* `neighs[i]`.
* `stiffness_flat[i*4+k]`.
* `bkSlots` or another reverse-adjacency mechanism.
* one atom work-item gathering all four ports.
* one apply kernel moving each atom once.
* collision gather per atom.

## Reconsider

### Neighbor angular effective mass

For the directed residual

[
x_i+R_i a_{ik}-x_j,
]

remove the guessed neighbor angular term. The denominator is

[
w_i+w_j+
(r_{ik}\times n)^TI_i^{-1}(r_{ik}\times n)
+\tilde\alpha.
]

### Automatic halving for node neighbors

Do not halve solely because `j_isnode`.

If reciprocal directed ports exist, they are two terms in the ARAP-like energy. Their relative weighting should be encoded in `stiffness_flat`, not inferred at runtime.

### Classify current kernels more accurately

* `compute_ports_cluster_rigid()`:
  atom-centered Jacobi XPBD/PBD-like gather.
* `compute_ports_cluster_rigid_shapematch()`:
  atom-centered ARAP/PD local rotation projection followed by XPBD-like translation.
* `compute_ports_cluster_rigid_substep_optimized()`:
  atom-centered rotational Newton/VBD-like block solve followed by XPBD-like translation.
* `apply_corrections_rigid_ports()`:
  common atom-centered gather/apply backend.
* `compute_collision_cluster_rigid()`:
  atom-centered Jacobi collision projection.

The current code is therefore already a **hybrid family**, not one pure solver.

---

# 11. The most natural fair comparison

Keep exactly the same topology, buffers, broad phase, and apply stage. Replace only the per-atom compute kernel:

```cpp
switch(method) {

case XPBD:
    gather scalar compliant corrections
    -> dx_node, dtheta_node, recoil;
    break;

case SOFT_IMPULSE:
    gather implicit spring impulses
    -> dv_node, domega_node, velocity_recoil;
    break;

case PD:
    local ARAP rotation projection;
    gather PD translation RHS and diagonal;
    -> x_new, q_new;
    break;

case VBD:
    gather 6D gradient and 6x6 Hessian;
    solve one body block;
    -> dx_node, dtheta_node;
    break;

case AVBD:
    same as VBD with lambda + rho*C;
    break;
}
```

My expectation is now slightly different from my previous answer:

* **VBD** is not special because it alone gathers all four ports. Your methods all gather all four ports.
* VBD is special because it combines all those ports into **one coupled atom-level Hessian solve** rather than summing independently calculated per-port corrections.
* Your current Newton rotation kernel already captures much of that advantage for the rotational part.
* The most informative next implementation is probably a full atom-centered (6\times6) VBD block and an atom-centered soft-impulse Jacobi kernel, compared against the current XPBD-like kernel.

The corrected source reference is here: [RRsp3(2).cl](sandbox:/mnt/data/RRsp3%282%29.cl).

[1]: https://box2d.org/posts/2024/02/solver2d/?utm_source=chatgpt.com "Solver2D"
[2]: https://arxiv.org/html/2403.06321v3?utm_source=chatgpt.com "Vertex Block Descent"
[3]: https://graphics.cs.utah.edu/research/projects/avbd/Augmented_VBD-SIGGRAPH25.pdf?utm_source=chatgpt.com "Augmented Vertex Block Descent"

---

# USER

OK, now we are on the same page, now try to focus again on key aspects
1) interaction between outer (soft O(n^2) steps) and inner loop (hard O(n) substreps), the intrial tem M/dt^2 or its analogs
2) momentum conservtion
3) you currently consider only the ports (bonds) but I also want to consider colisions with are outside the ports (unlike Electrostatics which is long range, the collison is contact and can be efficiently accelerated by AABB, do not focus on the accelration of braod-space collisions, that is another topic, jsut assume that each atom can easily find all atoms which are within interaction cuoff radius (~4A) using some efficient spatial map), but we need to agreagate the collision displacements or impulese with the bond displacement/impulese (and torques) before we update position and rotation (and associate linear and angular velocity) of each atom. That is crucial to make right to ensure stability and conservation of linear and anglualr momentum (realistic recoils Newton 3rd law).


---

# Chat GPT 5.5

Yes. The clean formulation is a **three-level structure**:

[
\boxed{
\text{soft macrostep}
;\supset;
\text{hard microsteps}
;\supset;
\text{one atom-centered gather and update}
}
]

At every hard microstep, **ports and collisions must contribute to the same translational and rotational accumulators before any atom is moved**.

---

# 1. Outer soft step versus inner hard microsteps

Let

* (H): expensive soft-force macrostep,
* (s): number of cheap hard microsteps,
* (h=H/s): hard microstep,
* (F_s): long-range electrostatic, attractive vdW/Morse, soft torsions,
* (U_h): port and collision energy.

A physically consistent multiple-time-step scheme is the symmetric kick–hard–kick split:

[
v_i
\leftarrow
v_i+\frac{H}{2m_i}F_{s,i}(x^n),
]

[
\omega_i
\leftarrow
\omega_i+\frac{H}{2}I_i^{-1}\tau_{s,i}(x^n),
]

then perform (s) hard microsteps of duration (h), followed by

[
v_i
\leftarrow
v_i+\frac{H}{2m_i}F_{s,i}(x^{n+1}),
]

[
\omega_i
\leftarrow
\omega_i+\frac{H}{2}I_i^{-1}\tau_{s,i}(x^{n+1}).
]

The force at the end of one macrostep is cached and reused as the first half-kick of the next one. Therefore, after initialization, this needs only **one new expensive soft-force evaluation per macrostep**. This is the same multiple-time-scale operator-splitting principle underlying reversible RESPA. ([Columbia University][1])

```cpp
void macroStep(float H, int nsub)
{
    applySoftKick(0.5f*H, Fsoft_cached, Tsoft_cached);

    const float h = H/nsub;

    for(int isub=0; isub<nsub; isub++){
        hardMicrostep(h);
    }

    evaluateSoftForces(Fsoft_new, Tsoft_new);

    applySoftKick(0.5f*H, Fsoft_new, Tsoft_new);

    swap(Fsoft_cached, Fsoft_new);
    swap(Tsoft_cached, Tsoft_new);
}
```

There is no extrapolated soft force. Both evaluations occur at actual configurations.

## What happens inside a hard microstep

This depends on the method.

### PD or VBD

Construct the free inertial target

[
x_i^\star=x_i+h v_i,
]

[
q_i^\star=\exp\left(\frac h2\omega_i\right)q_i.
]

Then approximately minimize

[
\boxed{
\Phi_h(x,q)
===========

\sum_i
\frac{m_i}{2h^2}|x_i-x_i^\star|^2
+
\sum_i
\frac{1}{2h^2}
\theta_i^T I_i\theta_i
+
U_{\rm port}(x,q)
+
U_{\rm coll}(x)
}
]

where

[
\theta_i=\operatorname{Log}(q_i(q_i^\star)^{-1}).
]

The inertial regularization is therefore

[
\boxed{
M/h^2,
}
]

not (M/H^2).

Projective Dynamics solves this implicit variational problem by local projections plus a global quadratic solve, while VBD minimizes the implicit-Euler objective by local vertex/body block updates. ([ACM Digital Library][2])

### XPBD

The same balance between inertia and stiffness appears through

[
\boxed{
\tilde\alpha=\frac{\alpha}{h^2}
=\frac{1}{K h^2}.
}
]

Thus a constraint denominator contains

[
w_{\rm eff}+\frac{1}{K h^2}.
]

This is XPBD’s analogue of (M/h^2+K). ([matthias-research.github.io][3])

### Soft-step impulse method

The same ratio appears as

[
\boxed{
h^2K,w_{\rm eff}
================

\frac{K}{m_{\rm eff}/h^2}.
}
]

The denominator of a local implicit spring impulse is

[
1+(hD+h^2K)w_{\rm eff}.
]

So the inertial regularization has not disappeared. It is encoded locally through the effective mass and (h^2K), rather than being assembled as a global matrix.

Catto’s Soft Step similarly combines soft spring–damper constraints with substepping; smaller substeps allow larger effective constraint frequencies. ([box2d.org][4])

## Crucial distinction

If the inner loop solves one fixed target (x^\star) repeatedly, it is performing **iterations of one macrostep problem** and must use (M/H^2).

If each inner pass:

1. advances time by (h),
2. updates (x,q,v,\omega),
3. recomputes port and collision geometry,

then it is a genuine microstep and must use

[
M/h^2.
]

---

# 2. Common atom-centered port and collision representation

For each atom (i), global dynamical state is only

[
x_i,\quad q_i,\quad v_i,\quad\omega_i.
]

Every hard microstep creates temporary accumulators:

[
J_i \quad\text{linear reaction},
]

[
L_i \quad\text{angular reaction}.
]

For SI these have units of physical impulse and angular impulse.

For XPBD they are position-level pseudo-impulses, but the gather structure is identical.

---

# 3. Directed port reaction

For port (k) of atom (i), pointing toward neighbor-center (j),

[
r_{ik}=R(q_i)a_{ik},
]

[
t_{ik}=x_i+r_{ik},
]

[
e_{ik}=x_j-t_{ik}.
]

Define

[
C=|e|,
\qquad
n=\frac{e}{|e|}.
]

A scalar reaction magnitude (p) produces a central reaction

[
P=p,n.
]

It acts:

* on atom (i), at its port tip;
* on atom (j), at its center.

Therefore:

[
J_i^{\rm port}\mathrel{+}=P,
]

[
L_i^{\rm port}\mathrel{+}=r_{ik}\times P,
]

[
J_j^{\rm recoil}\mathrel{-}=P.
]

There is no direct torque on (j), because this directed interaction acts on its center. Its orientation is controlled by its own outgoing ports.

This corresponds exactly to your current atom-centered design:

```text
dpos_node[i]        own outgoing port reactions
drot_node[i]        own outgoing port torques
dpos_neigh[slot]    recoil sent to target neighbor
bkSlots[i]          gather incoming recoils
```

---

# 4. Port momentum conservation

For one directed port interaction,

[
J_i+J_j=P-P=0.
]

Therefore linear momentum is exactly conserved.

The total angular impulse about an arbitrary origin is

[
\Delta L_{\rm total}
====================

x_i\times P
+
r_{ik}\times P
--------------

x_j\times P.
]

Thus

[
\Delta L_{\rm total}
====================

(t_{ik}-x_j)\times P.
]

But

[
t_{ik}-x_j=-e
]

and

[
P\parallel e.
]

Therefore

[
\boxed{
\Delta L_{\rm total}=0.
}
]

So the directed atom-port representation conserves angular momentum perfectly provided that:

1. the neighbor receives exactly (-P);
2. the central atom receives torque (r\times P);
3. (P) is parallel to the tip-to-center error (e).

This last condition matters.

A scalar radial port reaction preserves angular momentum. An arbitrary (3\times3) vector correction can produce

[
e\times P\ne0
]

and inject net angular momentum at finite port error.

Therefore, given your priorities, I recommend retaining a **scalar port constraint along (n)** rather than using an unconstrained vector impulse obtained from a full anisotropic (3\times3) inversion.

---

# 5. Sphere collision reaction

For a collision pair (i,j),

[
d=x_i-x_j,
\qquad
r=|d|,
\qquad
n=\frac{d}{r},
]

[
\delta=R_{ij}-r.
]

The collision is active when

[
\delta>0.
]

A repulsive scalar reaction (p_c\ge0) gives

[
P_c=p_c n.
]

Accumulate

[
J_i^{\rm coll}\mathrel{+}=P_c,
]

[
J_j^{\rm coll}\mathrel{-}=P_c.
]

Because the collision force passes through both centers,

[
L_i^{\rm coll}=L_j^{\rm coll}=0.
]

Linear momentum conservation is immediate:

[
P_c-P_c=0.
]

Angular momentum conservation is also immediate:

[
(x_i-x_j)\times P_c
===================

# d\times(p_cn)

0.

]

Thus ports and sphere collisions have exactly the same reaction structure:

| Interaction           |              Linear reaction |                  Torque |
| --------------------- | ---------------------------: | ----------------------: |
| Port (i\rightarrow j) |     (+P) on (i), (-P) on (j) | (r_{ik}\times P) on (i) |
| Sphere collision      | (+P_c) on (i), (-P_c) on (j) |                    zero |

---

# 6. Aggregate everything before updating the atom

For each atom:

[
\boxed{
J_i
===

J_i^{\rm own\ ports}
+
J_i^{\rm incoming\ recoils}
+
J_i^{\rm collisions}
}
]

and

[
\boxed{
L_i
===

L_i^{\rm own\ port\ torques}
+
L_i^{\rm collision\ torques}.
}
]

For spherical collisions, the second collision term is zero.

Only after all contributions have been gathered do you update the atom.

## Sequential/soft impulses

[
v_i\leftarrow v_i+m_i^{-1}J_i,
]

[
\omega_i
\leftarrow
\omega_i+I_i^{-1}L_i.
]

Then integrate once:

[
x_i\leftarrow x_i+h v_i,
]

[
q_i\leftarrow
\operatorname{normalize}
\left[
\exp\left(\frac h2\omega_i\right)q_i
\right].
]

## XPBD

Treat (J_i,L_i) as position pseudo-impulses:

[
x_i\leftarrow x_i+m_i^{-1}J_i,
]

[
q_i\leftarrow
\exp\left(\frac12I_i^{-1}L_i\right)q_i.
]

Then update velocities using the accumulated changes:

[
v_i^{n+1}
=========

v_i^\star+\frac{\Delta x_i}{h},
]

[
\omega_i^{n+1}
==============

\omega_i^\star+\frac{\Delta\theta_i}{h}.
]

Using accumulated deltas rather than subtracting two large absolute coordinates is also numerically preferable in single precision. Catto describes the same improvement for XPBD-style velocity reconstruction. ([box2d.org][4])

---

# 7. Soft-step SI formulas for both interactions

## Port impulse

The port separation rate is

[
u
=

n\cdot
\left[
v_j-
\left(
v_i+\omega_i\times r_{ik}
\right)
\right].
]

Positive (u) means the port and neighbor center are separating.

The effective inverse mass is

[
\boxed{
w_{\rm port}
============

m_i^{-1}
+
m_j^{-1}
+
(r_{ik}\times n)^T
I_i^{-1}
(r_{ik}\times n).
}
]

For stiffness (K_p), radial damping (D_p), and microstep (h), the local backward-Euler impulse is

[
\boxed{
p_p
===

\frac{
hK_p C
+
h(D_p+hK_p)u
}{
1+
(hD_p+h^2K_p)w_{\rm port}
}.
}
]

Then

[
P_p=p_p n.
]

This impulse can become negative when the port is closing too quickly. That is appropriate for a bilateral damped bond.

Use **radial damping only**, through (u=n\cdot v_{\rm relative}). Full vector damping proportional to the complete tip velocity difference can generate a noncentral force and therefore a net angular-momentum change.

## Collision impulse

Define separating velocity

[
u_c=n\cdot(v_i-v_j).
]

The effective inverse mass is

[
w_c=m_i^{-1}+m_j^{-1}.
]

For penetration (\delta>0),

[
\boxed{
p_c
===

\max\left[
0,
\frac{
hK_c\delta
----------

h(D_c+hK_c)u_c
}{
1+
(hD_c+h^2K_c)w_c
}
\right].
}
]

Then

[
P_c=p_c n.
]

This is a soft unilateral collision:

* penetration creates repulsion;
* approaching velocity increases repulsion;
* rapid separation can reduce the impulse to zero;
* no attractive contact impulse is allowed.

---

# 8. XPBD formulas for both interactions

For a port constraint

[
C_p=|e|,
]

the standard XPBD multiplier update is

[
\Delta\lambda
=============

\frac{
-C_p-\tilde\alpha_p\lambda
}{
w_{\rm port}+\tilde\alpha_p
},
]

with

[
\tilde\alpha_p=\frac{1}{K_ph^2}.
]

Define

[
p_p=-\Delta\lambda.
]

Then use the same central pseudo-impulse

[
P_p=p_p n.
]

Without persistent multipliers,

[
\boxed{
p_p
===

\frac{C_p}
{w_{\rm port}+1/(K_ph^2)}.
}
]

For collision penetration,

[
\boxed{
p_c
===

\frac{\delta}
{w_c+1/(K_ch^2)}.
}
]

Clamp it to be nonnegative.

The only difference from SI is that XPBD accumulates **position corrections**, while SI accumulates **velocity impulses**.

The port and collision contributions still go into the same per-atom (J_i,L_i) buffers before applying the update.

---

# 9. PD and VBD accumulation

Use the hard energy

[
U_h
===

U_{\rm port}+U_{\rm coll},
]

with

[
U_{\rm port}
============

\frac12
\sum_{i,k}
K_{ik}
\left|
x_j-x_i-R_i a_{ik}
\right|^2,
]

and

[
U_{\rm coll}
============

\frac12
\sum_{(i,j)}
K_{c,ij}
\max(0,R_{ij}-r_{ij})^2.
]

For each atom, gather both kinds of contributions before solving its block.

## VBD atom block

Let

[
z_i=
\begin{pmatrix}
x_i\
\theta_i
\end{pmatrix}.
]

Gather

[
g_i
===

g_i^{\rm inertia}
+
g_i^{\rm ports}
+
g_i^{\rm incoming}
+
g_i^{\rm collisions},
]

and

[
H_i
===

H_i^{\rm inertia}
+
H_i^{\rm ports}
+
H_i^{\rm incoming}
+
H_i^{\rm collisions}.
]

The inertial block is

[
\boxed{
H_i^{\rm inertia}
=================

\begin{pmatrix}
m_iI/h^2 & 0\
0 & I_i/h^2
\end{pmatrix}.
}
]

Then solve one (6\times6) system:

[
H_i\Delta z_i=-g_i.
]

Only then update (x_i,q_i).

Collisions add only translation blocks for spherical atoms:

[
g_{x,i}^{\rm coll}=-K_c\delta n,
]

[
H_{xx,i}^{\rm coll}\approx K_c nn^T.
]

## Momentum caveat for PD/VBD

The underlying global implicit equations are derived from translation- and rotation-invariant energies. A fully converged global solve has balanced internal reactions.

However, a finite number of independent atom-local Jacobi or VBD updates does **not generally preserve pairwise momentum exactly**:

[
\Delta z_i=-H_i^{-1}g_i
]

and

[
\Delta z_j=-H_j^{-1}g_j
]

need not be equal-and-opposite because (H_i\ne H_j).

Thus:

* SI and XPBD with explicit pair reactions give the cleanest exact momentum bookkeeping.
* PD/VBD may give better local convergence, but incomplete solves can introduce small rigid-body drift.

AVBD improves convergence for high stiffness ratios and hard constraints by adding augmented multipliers, but it does not remove this finite-iteration momentum issue automatically. ([graphics.cs.utah.edu][5])

---

# 10. Collision evaluation: pair once versus symmetric double gather

There are two valid GPU patterns.

## A. Compute every collision pair once

Compute (P_{ij}) once and store:

[
+P_{ij}\quad\text{for }i,
]

[
-P_{ij}\quad\text{for }j.
]

Then both atoms gather their slots.

Advantages:

* exactly the same floating-point (P_{ij}) is used on both atoms;
* exact Newton-third-law recoil;
* asymmetric neighbor maps cannot break conservation.

Cost:

* temporary pair/recoil storage;
* reverse mapping or atomics.

## B. Symmetric atom-centered double evaluation

Atom (i) scans (j) and computes its own contribution. Atom (j) independently scans (i).

This remains analytically momentum-conserving if:

1. both neighbor lists contain the pair;
2. both use identical positions, radii, masses and parameters;
3. both evaluate the same scalar (p_{ij});
4. the directions are exactly opposite.

This needs no collision recoil buffer.

### About the current factor `0.5`

Your collision kernel currently computes

```cpp
dl = penetration / (wi + wj);
dl *= 0.5f;
sum += n * (dl * wi);
```

Because each atom only writes its **own** correction, the duplicated pair evaluation does not itself double the correction on either atom.

Therefore the `0.5` is not required merely because both atoms evaluate the pair.

It should be interpreted as an explicit collision under-relaxation parameter:

```cpp
dl *= collision_relaxation;  // e.g. 0.5
```

That is perfectly reasonable for stability, but it should not be described as double-counting compensation.

---

# 11. Common atom-centered microstep pseudocode

```cpp
struct Reaction {
    Vec3 linear;    // impulse or pseudo-impulse
    Vec3 angular;
};

void computeHardReactions(float h, Method method)
{
    clear(ownReaction);
    clear(portRecoil);
    clear(collisionReaction);

    // One thread per central oriented atom
    parallel_for(int i=0; i<nAtoms; i++){

        Reaction acc = {};

        // Outgoing ports
        for(int k=0; k<nPorts[i]; k++){
            int j = neigh[i][k];
            if(j < 0) continue;

            Vec3 r = rotate(q[i], portLocal[i][k]);
            Vec3 e = x[j] - (x[i] + r);

            float l = length(e);
            if(l < eps) continue;

            Vec3 n = e/l;

            float p = computePortReaction(
                method, i, j, r, n, l, h
            );

            Vec3 P = p*n;

            acc.linear  += P;
            acc.angular += cross(r,P);

            portRecoil[i][k] = -P;
        }

        ownReaction[i] = acc;

        // Sphere collisions
        Vec3 Jcoll = {};

        for(int j : nearbyAtoms(i)){
            if(j == i || excluded(i,j)) continue;

            Vec3 d = x[i]-x[j];
            float r = length(d);
            float penetration = radius[i]+radius[j]-r;

            if(penetration <= 0.0f) continue;

            Vec3 n = d/r;

            float p = computeCollisionReaction(
                method, i, j, n, penetration, h
            );

            Jcoll += p*n;
        }

        collisionReaction[i] = Jcoll;
    }
}
```

Apply everything once:

```cpp
void applyHardReactionsSI(float h)
{
    parallel_for(int i=0; i<nAtoms; i++){

        Vec3 J =
              ownReaction[i].linear
            + collisionReaction[i];

        for(int slot : incomingPortSlots(i)){
            J += portRecoil[slot];
        }

        Vec3 L = ownReaction[i].angular;

        v[i]     += invMass[i] * J;
        omega[i] += invInertia[i] * L;

        x[i] += h*v[i];
        q[i]  = integrateQuaternion(q[i],omega[i],h);
    }
}
```

For XPBD, the same gather is used, but the final operation is:

```cpp
dx     = invMass[i]    * J;
dtheta = invInertia[i] * L;

x[i] += dx;
q[i]  = applyRotation(q[i],dtheta);

v[i]     = vFree[i]     + dx/h;
omega[i] = omegaFree[i] + dtheta/h;
```

---

# 12. Operations that can silently break conservation

Even with correct pair reactions, these later operations can destroy momentum conservation:

### Per-atom clipping after gathering

Bad:

```cpp
dx[i] = clampLength(dx[i], maxMove);
```

The opposite atom may be clipped differently.

Better: clamp or scale the scalar pair reaction (p_{ij}) before writing both recoils.

### Different relaxation on the two pair endpoints

The same scalar relaxation factor must multiply both (+P) and (-P).

### Full vector port damping

Use radial damping

[
D,n n^T v_{\rm relative},
]

not arbitrary tangential damping, unless you explicitly add compensating torques.

### Independent position heavy-ball terms

Your current

```cpp
move = dx*relaxation + dpos_mom*beta;
```

is not manifestly pairwise.

It may preserve zero total linear momentum if the stored displacement momentum already sums to zero, but angular conservation is much less transparent.

For the momentum-conserving version, maintain actual (v,\omega) and apply damping to:

* pair-relative radial velocities, or
* internal velocity modes after removing center-of-mass translation and rigid rotation.

### Fixed atoms

If one atom is fixed, the reaction is physically transferred to the external environment. Momentum of the mobile subsystem is then not expected to be conserved.

---

# Recommended architecture

For your immediate implementation I would use:

[
\boxed{
\text{cached soft half-kick}
\rightarrow
s\times
\begin{cases}
\text{gather port impulses}\
+\text{gather collision impulses}\
+\text{gather recoils}\
+\text{one atom velocity/rotation update}\
+\text{drift by }h
\end{cases}
\rightarrow
\text{new soft half-kick}
}
]

with:

* soft-step local implicit SI for ports and collisions;
* one scalar central reaction per interaction;
* no warm-start storage initially;
* no per-atom clipping;
* radial damping only;
* symmetric collision neighbor lists or one-evaluation pair buffers;
* one common (J_i,L_i) application per atom.

That version gives the most transparent combination of:

* very long soft-force macrosteps;
* stable cheap hard microsteps;
* exact Newton-third-law recoils;
* exact linear momentum conservation;
* exact angular momentum conservation for every individual port and sphere collision.

The relevant kernel is [RRsp3(2).cl](sandbox:/mnt/data/RRsp3%282%29.cl).

[1]: https://www.columbia.edu/cu/chemistry/groups/berne/papers/jcp_97_1990_1992.pdf?utm_source=chatgpt.com "Reversible multiple time scale molecular dynamics"
[2]: https://dl.acm.org/doi/10.1145/2601097.2601116?utm_source=chatgpt.com "Projective dynamics: fusing constraint projections for fast ..."
[3]: https://matthias-research.github.io/pages/publications/XPBD.pdf?utm_source=chatgpt.com "XPBD: Position-Based Simulation of Compliant ..."
[4]: https://box2d.org/posts/2024/02/solver2d/?utm_source=chatgpt.com "Solver2D"
[5]: https://graphics.cs.utah.edu/research/projects/avbd/?utm_source=chatgpt.com "Augmented Vertex Block Descent (AVBD)"

---

# USER

this is what we produced based on your instructions, but it crashes with NaNs and huge velocities can you find errors, generally correct it and improve it?

---

# Chat GPT 5.5

The primary failure is clear and deterministic: **the port impulse has the wrong sign**.

With your definitions

[
e=x_j-(x_i+r_i),\qquad n=e/|e|,
]

and

[
u=n\cdot\left[v_j-(v_i+\omega_i\times r_i)\right]=\dot C,
]

you apply (+P=+pn) to the port atom and (-P) to the neighbor. This changes the separation rate as

[
u' = u-w_{\rm eff}p.
]

Therefore, for (C>0), the restoring impulse must have (p>0). The current minus sign makes every port push its atom **away** from its neighbor, exponentially amplifying the error. 

I prepared corrected source, regression tests, and a unified patch:

* [Corrected `rigid_atom_ff.c`](sandbox:/mnt/data/rigid_atom_ff_fixed.c)
* [Corrected tests](sandbox:/mnt/data/test_rigid_atom_ff_fixed.c)
* [Unified patch](sandbox:/mnt/data/rigid_atom_ff_complete_fix.patch)

## 1. Catastrophic port-sign error

Current code:

```c
float p = -( a * C / H + ( a + d ) * u )
        / ( wport * ( 1.0f + a + d ) );
```

Correct local implicit impulse:

```c
float p = ( a * C / h + ( a + damp ) * u )
        / ( wport * ( 1.0f + a + damp ) );
```

The derivation is:

[
u'=u-wp,
]

[
C'=C+hu',
]

and for an implicit spring-damper impulse,

[
p=h(kC'+cu').
]

Solving gives

[
p=
\frac{hkC+h(c+hk)u}
{1+(hc+h^2k)w}.
]

Introducing

[
a=h^2kw,\qquad d=hcw
]

gives

[
\boxed{
p=
\frac{aC/h+(a+d)u}
{w(1+a+d)}
}.
]

Your previous expression had exactly the opposite restoring direction.

This explains the output:

```text
vel=(-5111, ...)
vel=(+8222, ...)
```

after the first macrostep. Each microstep enlarges the port error, which enlarges the next impulse.

## 2. Use the microstep (h), not macrostep (H)

Both hard solvers receive `h`, but ignore it and divide positional errors by `w->H`. 

For a genuine hard microstep,

[
h=\frac{H}{n_{\rm sub}},
]

the inertial stiffness is

[
\frac{M}{h^2},
]

and the velocity bias is proportional to

[
\frac{C}{h}.
]

Thus:

```c
a * C / h
```

not:

```c
a * C / H
```

The collision formula needs the same correction:

```c
float p =
    ( a * penetration / h - ( a + damp ) * separatingVelocity )
    / ( effectiveInvMass * ( 1.0f + a + damp ) );
```

A subtle point: with fixed dimensionless `a`, velocities scale approximately as (C/h), so more substeps can produce larger printed velocities while producing bounded displacements (h v). That is not itself instability.

For physical stiffness independent of `nsub`, store (k,c), and calculate:

```c
float a    = h * h * k * wEff;
float damp = h * c * wEff;
```

Then changing `nsub` changes integration accuracy but not the underlying spring parameters.

## 3. Soft force is applied twice

The current macrostep performs:

1. old-force half kick: (H/2),
2. old-force kick in every microstep: total (H),
3. new-force half kick: (H/2).

So the soft force contributes a total impulse corresponding to approximately

[
2H
]

per macrostep rather than (H). The full soft kick inside the microstep loop must be removed. 

Correct splitting:

```c
applySoftKick(w, 0.5f * H, oldForce, oldTorque);

for (int isub = 0; isub < nsub; isub++)
{
    updateCollisions(w);
    hardMicrostep(w, h);
}

evaluateSoftForces(w);

applySoftKick(w, 0.5f * H, newForce, newTorque);
```

After initialization, the ending force becomes the cached starting force for the next step, so this still costs one new (O(N^2)) evaluation per macrostep.

Call this once after setting the initial geometry:

```c
RAffWorld_EvaluateSoftForces(w);
```

Otherwise the first macrostep starts with a zero cached force and receives only the final half kick.

## 4. The Lennard–Jones force is reversed

The current code uses:

```c
float fmag = 12.0f * eps / r * (sr11 - sr5);
```

with `n` pointing from atom (i) to (j). This produces:

* attraction at short range,
* repulsion at long range,

which is the opposite of the intended interaction. The powers are also inconsistent with the stated (12-6) potential. 

For

[
U(r)=\epsilon
\left[
\left(\frac{r_{\min}}r\right)^{12}
----------------------------------

2\left(\frac{r_{\min}}r\right)^6
\right],
]

the force on atom (i), along (n=(x_j-x_i)/r), is

[
F_i=
\frac{12\epsilon}{r}
\left[
\left(\frac{r_{\min}}r\right)^6
-------------------------------

\left(\frac{r_{\min}}r\right)^{12}
\right]n.
]

Since hard collisions handle repulsion, the fixed version retains only the attractive branch:

```c
if (r <= softRmin)
    continue;

float sr   = softRmin / r;
float sr2  = sr * sr;
float sr6  = sr2 * sr2 * sr2;
float sr12 = sr6 * sr6;

float fmag = 12.0f * softEps / r * (sr6 - sr12);
```

This is positive, hence attractive, for (r>r_{\min}).

## 5. Update contacts during the macrostep

The current collision list is updated only at `s == 0`. If two atoms enter contact during later hard microsteps, the solver will not know until the next macrostep. 

For the test implementation:

```c
for (int s = 0; s < nsub; s++)
{
    RAffWorld_UpdateCollisions(w);
    hardMicrostep(w, h);
}
```

In the optimized version, the AABB list can be persistent with a Verlet/swept margin, but it must conservatively cover motion during the whole macrostep.

## 6. Exact coincident collisions were ignored

Current code:

```c
if (r < 1e-8f)
    continue;
```

leaves exactly coincident atoms permanently unresolved.

The fixed version assigns a deterministic fallback normal:

```c
if (r2 > 1e-16f)
{
    r = sqrtf(r2);
    n = d / r;
}
else
{
    r = 0.0f;
    n = ((i ^ j) & 1) ? axisX : axisY;
}
```

The applied impulses remain exactly (+P,-P), so this preserves linear momentum. Since the centers coincide, the orbital angular impulse is also zero.

## 7. Prefer one solver sweep per microstep initially

Your chain test already uses:

```c
w->nsub = 16;
w->nIter = 1;
```

so its failure is caused by the impulse sign, not by repeated iterations. 

Nevertheless, the world default was `nIter = 4`. A repeated iteration:

* keeps the same positions and port errors,
* applies another full (C/h) bias,
* changes velocities only.

That alters effective damping/stiffness and can become aggressive in a high-valence Jacobi gather. I changed the default to:

```c
w->nIter = 1;
```

Use more actual microsteps first:

```c
nsub ↑
nIter = 1
```

because those microsteps re-integrate and relinearize the port directions and contacts.

More than one velocity sweep can be added later if needed, but then it is worth using a proper accumulated-impulse Soft-Step formula or reducing the bias in subsequent sweeps.

## 8. Momentum conservation remains correct

For each directed port:

```c
Jown[i]       += P;
Lown[i]       += cross(r, P);
Jrecoil[port]  = -P;
```

The total linear impulse is

[
P-P=0.
]

The total angular impulse is

[
x_i\times P+r\times P-x_j\times P
=================================

(x_i+r-x_j)\times P.
]

Since (P) is parallel to the port error,

[
(x_i+r-x_j)\times P=0.
]

Thus both total linear and angular momentum are conserved per directed port, apart from floating-point roundoff.

Collisions likewise use exactly (+P,-P) along the center-to-center direction, so they conserve both momenta.

The fixed implementation preserves this structure. It does not clamp velocities independently per atom, because such clipping would destroy the equal-and-opposite recoil.

## 9. Reciprocal ports double the physical bond energy

Your tests create both:

[
i\rightarrow j
]

and

[
j\rightarrow i.
]

This is legitimate because each directed port constrains the orientation of its own atom. But it means a physical bond contributes two port energies.

For a physical stiffness (K_{\rm bond}), a natural assignment is approximately

[
K_{i\rightarrow j}
==================

# K_{j\rightarrow i}

\frac12K_{\rm bond}.
]

With algorithmic dimensionless parameters, start more gently for larger-valence structures:

```c
aBond = 0.25f to 0.5f;
dBond = 0.5f to 1.0f;
nIter = 1;
```

Then increase `aBond` after confirming stability.

## 10. Regression test improvement

The old bond test only inspected the result after 50 steps. I added a one-step restoring-direction test:

```c
float err0 = fabsf(length(x1 - x0) - bondLength);

RAffWorld_Step(w);

float err1 = fabsf(length(x1 - x0) - bondLength);

ENSURE(err1 < err0);
```

This catches the sign error immediately rather than waiting for overflow.

The chain test now checks `isfinite()` after every step, so it stops at the first invalid state instead of printing NaNs at steps 50, 100, and 150.

## Suggested first run

Apply the patch and force a clean rebuild. Your log shows only `test_rigid_atom_ff.c.o` being recompiled while `libbox3d.a` says “Nothing to be done,” so make sure the corrected solver source is actually part of the target:

```bash
cd /home/prokop/git/box3d

cp /path/to/rigid_atom_ff_fixed.c src/rigid_atom_ff.c
cp /path/to/test_rigid_atom_ff_fixed.c test/test_rigid_atom_ff.c

cmake --build build --clean-first --target test -j
./build/bin/test RigidAtomFFTest
```

For the first stable test:

```c
w->H      = 0.02f;
w->nsub   = 16;
w->nIter  = 1;
w->aBond  = 0.5f;
w->dBond  = 0.5f;
w->aColl  = 0.5f;
w->dColl  = 0.5f;
```

Then try `aBond = 1.0f` once the restoring sign is verified.

I could not compile against your local Box3D tree in this environment, but I validated the corrected impulse algebra and reproduced the six-atom chain update independently: reversing the sign causes exponential growth; the corrected sign remains finite and converges.
