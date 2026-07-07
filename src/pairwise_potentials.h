// PairwisePotentials: non-covalent pairwise force models for RigidAtomFF.
// See docs/prokop/RigidAtomFF_SI_short.plan.md — "Soft/Hard Potential Splitting".
//
// Each potential is split into:
//   - Hard parabolic part (stiff, evaluated every microstep via collision solver)
//   - Soft residual part (smooth, evaluated once per macrostep via fsoft)
//
// The split is C²-continuous at the equilibrium distance r0 where dE/dr = 0.

#pragma once

#include "box3d/math_functions.h"

#ifdef __cplusplus
extern "C" {
#endif

/// Per-pair parameters for the Coulomb + power-law repulsion split potential.
/// E(r) = qi_qj / r + C / r^n
///
/// Precompute once per pair (or when charges change). Store r0 and k_hard
/// in the collision pair struct for the hard solver. Use SoftForce() once
/// per macrostep to accumulate into fsoft.
typedef struct RAffPairPot
{
	float qi_qj;   ///< charge product q_i * q_j
	float C;       ///< repulsion coefficient (C > 0)
	int   n;       ///< repulsion power (4, 6, or 8 recommended)

	// Derived (computed by RAffPairPot_Setup):
	float r0;      ///< equilibrium distance where dE/dr = 0
	float k_hard;  ///< parabolic stiffness = d²E/dr² at r0
	float E_r0;    ///< E(r0) = qi_qj/r0 + C/r0^n
} RAffPairPot;

/// Compute r0, k_hard, E_r0 from qi_qj, C, n.
/// Requires qi_qj < 0 (attractive electrostatics) and C > 0.
/// If qi_qj >= 0 (no equilibrium), sets r0 = 0 and k_hard = 0 (no hard part).
void RAffPairPot_Setup( RAffPairPot* pot );

/// Hard force magnitude (positive = repulsive) for the parabolic part.
/// Returns 0 for r >= r0 (unilateral). For r < r0: F = k_hard * (r0 - r).
/// This is meant to be used as per-pair k and Rij in the collision impulse solver.
static inline float RAffPairPot_HardForce( const RAffPairPot* pot, float r )
{
	if ( pot->r0 <= 0.0f || r >= pot->r0 ) return 0.0f;
	return pot->k_hard * ( pot->r0 - r );
}

/// Soft residual force magnitude (positive = repulsive) along the pair axis.
/// F_soft = -dE_soft/dr. For r >= r0 this equals the full force.
/// For r < r0 the parabolic part is subtracted.
/// Call this once per macrostep and add to fsoft[i] -= F*n, fsoft[j] += F*n
/// (where n points from i to j).
float RAffPairPot_SoftForce( const RAffPairPot* pot, float r );

/// Soft residual potential energy. Useful for diagnostics.
float RAffPairPot_SoftEnergy( const RAffPairPot* pot, float r );

// ---- Future pairwise potentials can be added here ----

#ifdef __cplusplus
}
#endif
