// RigidAtomFF: atom-centered rigid-body dynamics with local implicit SI.
// See docs/prokop/RigidAtomFF_SI_short.plan.md for the full design.

#pragma once

#include "box3d/math_functions.h"
#include "mol_loader.h"

#ifdef __cplusplus
extern "C" {
#endif

#define RAFF_MAX_PORTS 4

typedef struct RAffAtom
{
	b3Vec3 pos;
	b3Quat quat;
	b3Vec3 vel;
	b3Vec3 omega;
	float mass;
	float invMass;
	b3Matrix3 invInertia;
	float radius;
	int nPorts;
} RAffAtom;

typedef struct RAffPort
{
	int neighIdx;
	b3Vec3 localDir;
	float kPort;     // per-port stiffness override (0 = use world default)
	float dPort;     // per-port damping override (0 = use world default)
} RAffPort;

typedef struct RAffCollisionPair
{
	int i;
	int j;
	float Rij;
} RAffCollisionPair;

typedef struct RAffWorld
{
	RAffAtom* atoms;
	int nAtoms;

	RAffPort* ports; // flat: nAtoms * RAFF_MAX_PORTS

	// Reverse adjacency: for each atom, which (sourceAtom, sourcePort) target it
	int* incomingCount;
	int* incomingOffset;
	int* incomingAtoms;
	int* incomingPorts;
	int totalIncoming;

	// Collision pairs (dynamic)
	RAffCollisionPair* collisionPairs;
	int nCollisionPairs;
	int collisionPairsCap;

	// Soft forces (cached, evaluated once per macrostep)
	b3Vec3* fsoft;
	b3Vec3* tsoft;

	// Temp accumulators (per microstep)
	b3Vec3* Jown;
	b3Vec3* Lown;
	b3Vec3* Jrecoil; // flat: nAtoms * RAFF_MAX_PORTS
	b3Vec3* Jcoll;

	// Parameters
	float H;          // soft macrostep
	int nsub;         // number of hard microsteps
	float kBond;      // port spring stiffness (physical units)
	float dBond;      // port radial damping (physical units)
	float kColl;      // collision spring stiffness
	float dColl;      // collision damping
	float relaxation; // Jacobi relaxation factor (0-1)

	// Soft force parameters
	float softCutoff;
	float softEps;   // LJ epsilon
	float softRmin;  // LJ rmin
} RAffWorld;

RAffWorld* RAffWorld_Create( int nAtoms );
void RAffWorld_Destroy( RAffWorld* w );

// Build reverse adjacency from current port topology.
void RAffWorld_BuildAdjacency( RAffWorld* w );

// Update collision pairs using simple O(N^2) broad phase.
void RAffWorld_UpdateCollisions( RAffWorld* w );

// Evaluate soft forces (Lennard-Jones attractive + Coulomb stub).
void RAffWorld_EvaluateSoftForces( RAffWorld* w );

// One full macrostep: soft half-kick -> nsub hard microsteps -> soft half-kick.
void RAffWorld_Step( RAffWorld* w );

// Diagnostic: compute total linear and angular momentum.
b3Vec3 RAffWorld_TotalLinearMomentum( RAffWorld* w );
b3Vec3 RAffWorld_TotalAngularMomentum( RAffWorld* w );

// Build RAffWorld atoms and ports from a parsed MOL file.
// offset shifts all atom positions. scale multiplies coordinates (e.g. 0.1 for Angstrom→nm).
// Uses covalent radii for atom sizes and atomic masses for inertia.
void RAffWorld_LoadMol( RAffWorld* w, const b3MolFile* mol, b3Vec3 offset, float scale );

#ifdef __cplusplus
}
#endif
