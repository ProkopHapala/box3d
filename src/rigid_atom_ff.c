#include "rigid_atom_ff.h"
#include "math_internal.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

// --- Create / Destroy ---

RAffWorld* RAffWorld_Create( int nAtoms )
{
	RAffWorld* w = (RAffWorld*)calloc( 1, sizeof( RAffWorld ) );
	w->nAtoms = nAtoms;
	w->atoms = (RAffAtom*)calloc( nAtoms, sizeof( RAffAtom ) );
	w->ports = (RAffPort*)calloc( (size_t)nAtoms * RAFF_MAX_PORTS, sizeof( RAffPort ) );
	w->incomingCount = (int*)calloc( nAtoms, sizeof( int ) );
	w->incomingOffset = (int*)calloc( nAtoms, sizeof( int ) );
	// incomingAtoms/Ports allocated after adjacency build
	w->fsoft = (b3Vec3*)calloc( nAtoms, sizeof( b3Vec3 ) );
	w->tsoft = (b3Vec3*)calloc( nAtoms, sizeof( b3Vec3 ) );
	w->Jown = (b3Vec3*)calloc( nAtoms, sizeof( b3Vec3 ) );
	w->Lown = (b3Vec3*)calloc( nAtoms, sizeof( b3Vec3 ) );
	w->Jrecoil = (b3Vec3*)calloc( (size_t)nAtoms * RAFF_MAX_PORTS, sizeof( b3Vec3 ) );
	w->Jcoll = (b3Vec3*)calloc( nAtoms, sizeof( b3Vec3 ) );
	w->lambdaPort = (float*)calloc( (size_t)nAtoms * RAFF_MAX_PORTS, sizeof( float ) );
	w->Hblock = (float*)calloc( (size_t)nAtoms * 36, sizeof( float ) );
	w->gblock = (float*)calloc( (size_t)nAtoms * 6, sizeof( float ) );
	w->solverMethod = RAFF_SOLVER_SI;

	// Default parameters — physical stiffness and damping
	// K must be large enough that h^2*K*w_eff ~ O(1) for effective correction
	w->H = 0.02f;
	w->nsub = 8;
	w->kBond = 1e6f;    // port spring stiffness
	w->dBond = 1e4f;    // port radial damping
	w->kColl = 5e6f;    // collision spring stiffness
	w->dColl = 5e4f;    // collision damping
	w->relaxation = 0.5f; // damp Jacobi overshoot from multiple constraints per atom
	w->softCutoff = 8.0f;
	w->softEps = 0.1f;
	w->softRmin = 3.0f;

	// Collision pairs — start with capacity for nAtoms^2/2
	w->collisionPairsCap = nAtoms * nAtoms / 2 + 16;
	w->collisionPairs = (RAffCollisionPair*)calloc( w->collisionPairsCap, sizeof( RAffCollisionPair ) );
	w->nCollisionPairs = 0;

	// Initialize atoms with identity orientation
	for ( int i = 0; i < nAtoms; i++ )
	{
		w->atoms[i].quat = b3Quat_identity;
		w->atoms[i].invInertia = b3Mat3_identity;
	}

	return w;
}

void RAffWorld_Destroy( RAffWorld* w )
{
	if ( w == NULL ) return;
	free( w->atoms );
	free( w->ports );
	free( w->incomingCount );
	free( w->incomingOffset );
	free( w->incomingAtoms );
	free( w->incomingPorts );
	free( w->collisionPairs );
	free( w->fsoft );
	free( w->tsoft );
	free( w->Jown );
	free( w->Lown );
	free( w->Jrecoil );
	free( w->Jcoll );
	free( w->lambdaPort );
	free( w->Hblock );
	free( w->gblock );
	free( w );
}

// --- Reverse adjacency ---

void RAffWorld_BuildAdjacency( RAffWorld* w )
{
	// Count incoming ports per atom
	memset( w->incomingCount, 0, w->nAtoms * sizeof( int ) );
	for ( int i = 0; i < w->nAtoms; i++ )
	{
		for ( int k = 0; k < w->atoms[i].nPorts; k++ )
		{
			int j = w->ports[i * RAFF_MAX_PORTS + k].neighIdx;
			if ( j >= 0 && j < w->nAtoms )
				w->incomingCount[j]++;
		}
	}

	// Compute offsets (prefix sum)
	w->totalIncoming = 0;
	for ( int i = 0; i < w->nAtoms; i++ )
	{
		w->incomingOffset[i] = w->totalIncoming;
		w->totalIncoming += w->incomingCount[i];
	}

	// Allocate flat arrays
	w->incomingAtoms = (int*)realloc( w->incomingAtoms, w->totalIncoming * sizeof( int ) );
	w->incomingPorts = (int*)realloc( w->incomingPorts, w->totalIncoming * sizeof( int ) );

	// Fill using a temporary cursor
	int* cursor = (int*)calloc( w->nAtoms, sizeof( int ) );
	for ( int i = 0; i < w->nAtoms; i++ )
	{
		for ( int k = 0; k < w->atoms[i].nPorts; k++ )
		{
			int j = w->ports[i * RAFF_MAX_PORTS + k].neighIdx;
			if ( j >= 0 && j < w->nAtoms )
			{
				int slot = w->incomingOffset[j] + cursor[j];
				w->incomingAtoms[slot] = i;
				w->incomingPorts[slot] = k;
				cursor[j]++;
			}
		}
	}
	free( cursor );
}

// --- Broad phase: simple O(N^2) pair search ---

void RAffWorld_UpdateCollisions( RAffWorld* w )
{
	w->nCollisionPairs = 0;
	for ( int i = 0; i < w->nAtoms; i++ )
	{
		for ( int j = i + 1; j < w->nAtoms; j++ )
		{
			// Skip bonded pairs (they have port constraints, not collision)
			bool bonded = false;
			for ( int k = 0; k < w->atoms[i].nPorts; k++ )
			{
				if ( w->ports[i * RAFF_MAX_PORTS + k].neighIdx == j )
				{
					bonded = true;
					break;
				}
			}
			if ( bonded ) continue;

			float Rij = w->atoms[i].radius + w->atoms[j].radius;
			b3Vec3 d = b3Sub( w->atoms[i].pos, w->atoms[j].pos );
			float r = b3Length( d );
			// Use a small margin to catch near-contact pairs early
			if ( r < Rij * 1.1f + w->softCutoff * 0.5f )
			{
				if ( w->nCollisionPairs < w->collisionPairsCap )
				{
					w->collisionPairs[w->nCollisionPairs].i = i;
					w->collisionPairs[w->nCollisionPairs].j = j;
					w->collisionPairs[w->nCollisionPairs].Rij = Rij;
					w->nCollisionPairs++;
				}
			}
		}
	}
}

// --- Soft forces: Lennard-Jones attractive + simple Coulomb stub ---

void RAffWorld_EvaluateSoftForces( RAffWorld* w )
{
	for ( int i = 0; i < w->nAtoms; i++ )
	{
		w->fsoft[i] = b3Vec3_zero;
		w->tsoft[i] = b3Vec3_zero;
	}

	float rcut2 = w->softCutoff * w->softCutoff;

	for ( int i = 0; i < w->nAtoms; i++ )
	{
		for ( int j = i + 1; j < w->nAtoms; j++ )
		{
			// Skip bonded pairs
			bool bonded = false;
			for ( int k = 0; k < w->atoms[i].nPorts; k++ )
			{
				if ( w->ports[i * RAFF_MAX_PORTS + k].neighIdx == j )
				{
					bonded = true;
					break;
				}
			}
			if ( bonded ) continue;

			b3Vec3 d = b3Sub( w->atoms[j].pos, w->atoms[i].pos );
			float r2 = b3Dot( d, d );
			if ( r2 > rcut2 || r2 < 1e-12f ) continue;

			float r = sqrtf( r2 );
			b3Vec3 n = b3MulSV( 1.0f / r, d );

			// Lennard-Jones: F = 12*eps/r * [ (rmin/r)^11 - (rmin/r)^5 ] * n
			// (repulsive + attractive, but we only use it for soft attraction here)
			float sr = w->softRmin / r;
			float sr5 = sr * sr * sr * sr * sr;
			float sr11 = sr5 * sr5 * sr;
			float fmag = 12.0f * w->softEps / r * ( sr11 - sr5 );

			// Clamp to avoid singularities at very short range
			if ( fmag > 100.0f ) fmag = 100.0f;
			if ( fmag < -100.0f ) fmag = -100.0f;

			b3Vec3 f = b3MulSV( fmag, n );
			w->fsoft[i] = b3Add( w->fsoft[i], f );
			w->fsoft[j] = b3Sub( w->fsoft[j], f );
		}
	}
}

// --- Hard microstep internals ---

static void applySoftKick( RAffWorld* w, float dt, b3Vec3* f, b3Vec3* t )
{
	for ( int i = 0; i < w->nAtoms; i++ )
	{
		w->atoms[i].vel = b3Add( w->atoms[i].vel, b3MulSV( dt * w->atoms[i].invMass, f[i] ) );
		w->atoms[i].omega = b3Add( w->atoms[i].omega, b3MulMV( w->atoms[i].invInertia, b3MulSV( dt, t[i] ) ) );
	}
}

// --- Gather-apply helpers ---

static void clearAccumulators( RAffWorld* w )
{
	for ( int i = 0; i < w->nAtoms; i++ )
	{
		w->Jown[i] = b3Vec3_zero;
		w->Lown[i] = b3Vec3_zero;
		w->Jcoll[i] = b3Vec3_zero;
	}
	for ( int i = 0; i < w->nAtoms * RAFF_MAX_PORTS; i++ )
		w->Jrecoil[i] = b3Vec3_zero;
}

// --- SI: Sequential Impulses (gather phase) ---

static void gatherPortImpulsesSI( RAffWorld* w, float h )
{
	for ( int i = 0; i < w->nAtoms; i++ )
	{
		RAffAtom* ai = &w->atoms[i];
		b3Vec3 Jacc = b3Vec3_zero;
		b3Vec3 Lacc = b3Vec3_zero;

		for ( int k = 0; k < ai->nPorts; k++ )
		{
			RAffPort* port = &w->ports[i * RAFF_MAX_PORTS + k];
			int j = port->neighIdx;
			if ( j < 0 || j >= w->nAtoms ) continue;

			RAffAtom* aj = &w->atoms[j];

			b3Vec3 r = b3RotateVector( ai->quat, port->localDir );
			b3Vec3 tip = b3Add( ai->pos, r );
			b3Vec3 e = b3Sub( aj->pos, tip );
			float C = b3Length( e );
			if ( C < 1e-8f ) continue;
			b3Vec3 n = b3MulSV( 1.0f / C, e );

			b3Vec3 tipVel = b3Add( ai->vel, b3Cross( ai->omega, r ) );
			float u = b3Dot( n, b3Sub( aj->vel, tipVel ) );

			b3Vec3 s = b3Cross( r, n );
			float wport = ai->invMass + aj->invMass + b3Dot( s, b3MulMV( ai->invInertia, s ) );

			// K/2 per directed port (reciprocal ports double bond energy)
			float K = port->kPort > 0.0f ? port->kPort : w->kBond * 0.5f;
			float D = port->dPort > 0.0f ? port->dPort : w->dBond * 0.5f;

			float hK = h * K;
			float hD = h * D;
			float h2K = h * hK;
			// Local implicit (backward Euler) impulse — positive sign restores bond
			float p = ( hK * C + h * ( D + hK ) * u ) / ( 1.0f + ( hD + h2K ) * wport );

			b3Vec3 P = b3MulSV( p, n );

			Jacc = b3Add( Jacc, P );
			Lacc = b3Add( Lacc, b3Cross( r, P ) );
			w->Jrecoil[i * RAFF_MAX_PORTS + k] = b3MulSV( -1.0f, P );
		}

		w->Jown[i] = Jacc;
		w->Lown[i] = Lacc;
	}
}

// --- XPBD: warm-started Lagrange multipliers (gather phase) ---

static void gatherPortImpulsesXPBD( RAffWorld* w, float h )
{
	for ( int i = 0; i < w->nAtoms; i++ )
	{
		RAffAtom* ai = &w->atoms[i];
		b3Vec3 Jacc = b3Vec3_zero;
		b3Vec3 Lacc = b3Vec3_zero;

		for ( int k = 0; k < ai->nPorts; k++ )
		{
			RAffPort* port = &w->ports[i * RAFF_MAX_PORTS + k];
			int j = port->neighIdx;
			if ( j < 0 || j >= w->nAtoms ) continue;

			RAffAtom* aj = &w->atoms[j];

			b3Vec3 r = b3RotateVector( ai->quat, port->localDir );
			b3Vec3 tip = b3Add( ai->pos, r );
			b3Vec3 e = b3Sub( aj->pos, tip );
			float C = b3Length( e );
			if ( C < 1e-8f ) continue;
			b3Vec3 n = b3MulSV( 1.0f / C, e );

			b3Vec3 s = b3Cross( r, n );
			float wport = ai->invMass + aj->invMass + b3Dot( s, b3MulMV( ai->invInertia, s ) );

			// XPBD compliance: alphaTilde = 1 / (K * h^2)
			// Using K/2 for reciprocal ports
			float K = port->kPort > 0.0f ? port->kPort : w->kBond * 0.5f;
			float alphaTilde = 1.0f / ( K * h * h );

			float* lambda = &w->lambdaPort[i * RAFF_MAX_PORTS + k];

			// XPBD: dlambda = -(C + alphaTilde * lambda) / (wport + alphaTilde)
			// C > 0 (stretched) → dlambda < 0. Impulse must pull atoms together → p > 0 along n.
			// Position correction: dx_i = -dlambda * grad_C_i = -dlambda * (-n) = dlambda * n
			// Since dlambda < 0, dx_i is along -n (toward i) — WRONG for our convention where n points tip→j.
			// Our C = |error| (always positive), n points tip→j. We want atom i to move toward j (+n).
			// So impulse on i: P = -dlambda/h * n (positive since dlambda < 0)
			float dlambda = -( C + alphaTilde * ( *lambda ) ) / ( wport + alphaTilde );
			*lambda += dlambda;

			float p = -dlambda / h;
			b3Vec3 P = b3MulSV( p, n );

			Jacc = b3Add( Jacc, P );
			Lacc = b3Add( Lacc, b3Cross( r, P ) );
			w->Jrecoil[i * RAFF_MAX_PORTS + k] = b3MulSV( -1.0f, P );
		}

		w->Jown[i] = Jacc;
		w->Lown[i] = Lacc;
	}
}

// --- PD: Projective Dynamics (gather phase) ---

static void gatherPortImpulsesPD( RAffWorld* w, float h )
{
	for ( int i = 0; i < w->nAtoms; i++ )
	{
		RAffAtom* ai = &w->atoms[i];
		b3Vec3 Jacc = b3Vec3_zero;
		b3Vec3 Lacc = b3Vec3_zero;

		for ( int k = 0; k < ai->nPorts; k++ )
		{
			RAffPort* port = &w->ports[i * RAFF_MAX_PORTS + k];
			int j = port->neighIdx;
			if ( j < 0 || j >= w->nAtoms ) continue;

			RAffAtom* aj = &w->atoms[j];

			b3Vec3 r = b3RotateVector( ai->quat, port->localDir );
			b3Vec3 tip = b3Add( ai->pos, r );
			b3Vec3 e = b3Sub( aj->pos, tip );
			float C = b3Length( e );
			if ( C < 1e-8f ) continue;
			b3Vec3 n = b3MulSV( 1.0f / C, e );

			b3Vec3 tipVel = b3Add( ai->vel, b3Cross( ai->omega, r ) );
			float u = b3Dot( n, b3Sub( aj->vel, tipVel ) );

			b3Vec3 s = b3Cross( r, n );
			float wport = ai->invMass + aj->invMass + b3Dot( s, b3MulMV( ai->invInertia, s ) );

			// PD: position-level correction with weight = K * h^2, plus velocity damping
			// Using K/2 for reciprocal ports
			float K = port->kPort > 0.0f ? port->kPort : w->kBond * 0.5f;
			float D = port->dPort > 0.0f ? port->dPort : w->dBond * 0.5f;
			float a = K * h * h * wport; // dimensionless stiffness
			float hD_w = h * D * wport;

			// Position correction + velocity damping
			float p = ( a * C / h + ( a + hD_w ) * u ) / ( ( 1.0f + a + hD_w ) * wport );

			b3Vec3 P = b3MulSV( p, n );

			Jacc = b3Add( Jacc, P );
			Lacc = b3Add( Lacc, b3Cross( r, P ) );
			w->Jrecoil[i * RAFF_MAX_PORTS + k] = b3MulSV( -1.0f, P );
		}

		w->Jown[i] = Jacc;
		w->Lown[i] = Lacc;
	}
}

// --- VBD: 6x6 block solve (gather phase) ---

// Solve 6x6 linear system H * x = g using Schur complement:
// H = [Htt Htr; Htr^T Hrr], split into 3x3 blocks.
// Solve: S * dtheta = (gr - Htr^T * Htt^{-1} * gt), where S = Hrr - Htr^T * Htt^{-1} * Htr
// Then: dx = Htt^{-1} * (gt - Htr * dtheta)
static void solve6x6( const float H[6][6], const float g[6], float dq[6] )
{
	// Extract 3x3 blocks from the 6x6 Hessian (row-major H → column-major b3Matrix3)
	// b3Matrix3 stores columns as cx, cy, cz. H[row][col] is row-major.
	// So Htt.cx (column 0 of Htt) = {H[0][0], H[1][0], H[2][0]}
	b3Matrix3 Htt = { { H[0][0], H[1][0], H[2][0] }, { H[0][1], H[1][1], H[2][1] }, { H[0][2], H[1][2], H[2][2] } };
	b3Matrix3 Htr = { { H[0][3], H[1][3], H[2][3] }, { H[0][4], H[1][4], H[2][4] }, { H[0][5], H[1][5], H[2][5] } };
	b3Matrix3 Hrt = { { H[3][0], H[4][0], H[5][0] }, { H[3][1], H[4][1], H[5][1] }, { H[3][2], H[4][2], H[5][2] } };
	b3Matrix3 Hrr = { { H[3][3], H[4][3], H[5][3] }, { H[3][4], H[4][4], H[5][4] }, { H[3][5], H[4][5], H[5][5] } };

	b3Vec3 gt = { g[0], g[1], g[2] };
	b3Vec3 gr = { g[3], g[4], g[5] };

	// Htt_inv
	b3Matrix3 HttInv = b3InvertMatrix( Htt );

	// Schur complement: S = Hrr - Hrt * HttInv * Htr
	b3Matrix3 HttInv_Htr = { b3MulMV( HttInv, Htr.cx ), b3MulMV( HttInv, Htr.cy ), b3MulMV( HttInv, Htr.cz ) };
	b3Matrix3 S = { b3Sub( Hrr.cx, b3MulMV( Hrt, HttInv_Htr.cx ) ),
					b3Sub( Hrr.cy, b3MulMV( Hrt, HttInv_Htr.cy ) ),
					b3Sub( Hrr.cz, b3MulMV( Hrt, HttInv_Htr.cz ) ) };

	// rhs = gr - Hrt * HttInv * gt
	b3Vec3 HttInv_gt = b3MulMV( HttInv, gt );
	b3Vec3 rhs = b3Sub( gr, b3MulMV( Hrt, HttInv_gt ) );

	// Solve S * dtheta = rhs
	b3Vec3 dtheta = b3Solve3( S, rhs );

	// dx = HttInv * (gt - Htr * dtheta)
	b3Vec3 dx = b3MulMV( HttInv, b3Sub( gt, b3MulMV( Htr, dtheta ) ) );

	dq[0] = dx.x; dq[1] = dx.y; dq[2] = dx.z;
	dq[3] = dtheta.x; dq[4] = dtheta.y; dq[5] = dtheta.z;
}

static void gatherPortImpulsesVBD( RAffWorld* w, float h )
{
	float h2 = h * h;

	for ( int i = 0; i < w->nAtoms; i++ )
	{
		RAffAtom* ai = &w->atoms[i];

		float H[6][6];
		float g[6];
		for ( int a = 0; a < 6; a++ )
		{
			g[a] = 0.0f;
			for ( int b = 0; b < 6; b++ )
				H[a][b] = 0.0f;
		}

		// Inertial term: H += M (mass matrix, NOT inverse mass)
		// Translation block: mass on diagonal
		H[0][0] += ai->mass;
		H[1][1] += ai->mass;
		H[2][2] += ai->mass;

		// Rotation block: I (inertia tensor)
		// b3Matrix3 is column-major, H is row-major, so transpose: H[i][j] = Imat.col[j].row[i]
		b3Matrix3 Imat = b3InvertMatrix( ai->invInertia );
		H[3][3] += Imat.cx.x; H[3][4] += Imat.cy.x; H[3][5] += Imat.cz.x;
		H[4][3] += Imat.cx.y; H[4][4] += Imat.cy.y; H[4][5] += Imat.cz.y;
		H[5][3] += Imat.cx.z; H[5][4] += Imat.cy.z; H[5][5] += Imat.cz.z;

		// Inertial gradient: g += M * v * h
		g[0] += ai->mass * ai->vel.x * h;
		g[1] += ai->mass * ai->vel.y * h;
		g[2] += ai->mass * ai->vel.z * h;
		b3Vec3 grot = b3MulMV( Imat, b3MulSV( h, ai->omega ) );
		g[3] += grot.x;
		g[4] += grot.y;
		g[5] += grot.z;

		// Port contributions
		for ( int k = 0; k < ai->nPorts; k++ )
		{
			RAffPort* port = &w->ports[i * RAFF_MAX_PORTS + k];
			int j = port->neighIdx;
			if ( j < 0 || j >= w->nAtoms ) continue;

			RAffAtom* aj = &w->atoms[j];

			b3Vec3 r = b3RotateVector( ai->quat, port->localDir );
			b3Vec3 tip = b3Add( ai->pos, r );
			b3Vec3 e = b3Sub( aj->pos, tip );
			float C = b3Length( e );
			if ( C < 1e-8f ) continue;
			b3Vec3 n = b3MulSV( 1.0f / C, e );

			// K/2 for reciprocal ports
			float K = port->kPort > 0.0f ? port->kPort : w->kBond * 0.5f;
			float D = port->dPort > 0.0f ? port->dPort : w->dBond * 0.5f;
			float Kh2 = K * h2;
			float Dh2 = D * h2;

			b3Vec3 rxn = b3Cross( r, n );
			float Jv[6] = { n.x, n.y, n.z, rxn.x, rxn.y, rxn.z };

			// Relative velocity along constraint
			b3Vec3 tipVel = b3Add( ai->vel, b3Cross( ai->omega, r ) );
			float u = b3Dot( n, b3Sub( aj->vel, tipVel ) );

			// Hessian: H += h^2 * K * J * J^T  (rank-1 update)
			for ( int a = 0; a < 6; a++ )
				for ( int b = 0; b < 6; b++ )
					H[a][b] += Kh2 * Jv[a] * Jv[b];

			// RHS: g -= h^2 * (K * C + D * u) * J (gradient sign; negated before solve)
			float stiffDamp = Kh2 * C + Dh2 * u;
			for ( int a = 0; a < 6; a++ )
				g[a] -= stiffDamp * Jv[a];

			// Recoil: use scalar approximation for neighbor
			b3Vec3 s = b3Cross( r, n );
			float wport = ai->invMass + aj->invMass + b3Dot( s, b3MulMV( ai->invInertia, s ) );
			float a_dim = K * h * h * wport;
			float p = a_dim * C / ( 1.0f + a_dim ) / ( h * wport );
			b3Vec3 P = b3MulSV( p, n );
			w->Jrecoil[i * RAFF_MAX_PORTS + k] = b3MulSV( -1.0f, P );
		}

		// Solve H * dq = -g (negate gradient to get RHS)
		for ( int a = 0; a < 6; a++ ) g[a] = -g[a];
		float dq[6];
		solve6x6( H, g, dq );

		// v_new = v + dq/h (ADD correction velocity)
		b3Vec3 dvel = { dq[0] / h, dq[1] / h, dq[2] / h };
		b3Vec3 domega = { dq[3] / h, dq[4] / h, dq[5] / h };

		b3Vec3 J = b3MulSV( ai->mass, dvel );
		b3Vec3 L = b3MulMV( Imat, domega );

		w->Jown[i] = J;
		w->Lown[i] = L;
	}
}

// --- AVBD: Augmented VBD (gather phase) ---

static void gatherPortImpulsesAVBD( RAffWorld* w, float h )
{
	float h2 = h * h;

	for ( int i = 0; i < w->nAtoms; i++ )
	{
		RAffAtom* ai = &w->atoms[i];

		float H[6][6];
		float g[6];
		for ( int a = 0; a < 6; a++ )
		{
			g[a] = 0.0f;
			for ( int b = 0; b < 6; b++ )
				H[a][b] = 0.0f;
		}

		// Inertial term (same as VBD — mass matrix, NOT inverse mass)
		H[0][0] += ai->mass;
		H[1][1] += ai->mass;
		H[2][2] += ai->mass;

		// Transpose: b3Matrix3 column-major → H row-major
		b3Matrix3 Imat = b3InvertMatrix( ai->invInertia );
		H[3][3] += Imat.cx.x; H[3][4] += Imat.cy.x; H[3][5] += Imat.cz.x;
		H[4][3] += Imat.cx.y; H[4][4] += Imat.cy.y; H[4][5] += Imat.cz.y;
		H[5][3] += Imat.cx.z; H[5][4] += Imat.cy.z; H[5][5] += Imat.cz.z;

		// Inertial gradient: g += M * v * h
		g[0] += ai->mass * ai->vel.x * h;
		g[1] += ai->mass * ai->vel.y * h;
		g[2] += ai->mass * ai->vel.z * h;
		b3Vec3 grot = b3MulMV( Imat, b3MulSV( h, ai->omega ) );
		g[3] += grot.x;
		g[4] += grot.y;
		g[5] += grot.z;

		for ( int k = 0; k < ai->nPorts; k++ )
		{
			RAffPort* port = &w->ports[i * RAFF_MAX_PORTS + k];
			int j = port->neighIdx;
			if ( j < 0 || j >= w->nAtoms ) continue;

			RAffAtom* aj = &w->atoms[j];

			b3Vec3 r = b3RotateVector( ai->quat, port->localDir );
			b3Vec3 tip = b3Add( ai->pos, r );
			b3Vec3 e = b3Sub( aj->pos, tip );
			float C = b3Length( e );
			if ( C < 1e-8f ) continue;
			b3Vec3 n = b3MulSV( 1.0f / C, e );

			float K = port->kPort > 0.0f ? port->kPort : w->kBond * 0.5f;
			float D = port->dPort > 0.0f ? port->dPort : w->dBond * 0.5f;
			float Kh2 = K * h2;
			float Dh2 = D * h2;

			b3Vec3 rxn = b3Cross( r, n );
			float Jv[6] = { n.x, n.y, n.z, rxn.x, rxn.y, rxn.z };

			b3Vec3 tipVel = b3Add( ai->vel, b3Cross( ai->omega, r ) );
			float u = b3Dot( n, b3Sub( aj->vel, tipVel ) );

			// Standard VBD terms: H += h^2 * K * J * J^T, g -= h^2 * (K*C + D*u) * J
			for ( int a = 0; a < 6; a++ )
				for ( int b = 0; b < 6; b++ )
					H[a][b] += Kh2 * Jv[a] * Jv[b];
			float stiffDamp = Kh2 * C + Dh2 * u;
			for ( int a = 0; a < 6; a++ )
				g[a] -= stiffDamp * Jv[a];

			// Augmented Lagrangian: H += lambda * J * J^T, g -= lambda * C * J
			// Energy term: lambda * C^2 / 2. Gradient: lambda * C * dC/dq = lambda * C * J
			float lambda = w->lambdaPort[i * RAFF_MAX_PORTS + k];
			for ( int a = 0; a < 6; a++ )
				for ( int b = 0; b < 6; b++ )
					H[a][b] += lambda * Jv[a] * Jv[b];
			for ( int a = 0; a < 6; a++ )
				g[a] -= lambda * C * Jv[a];

			// Recoil (same scalar approximation as VBD)
			b3Vec3 s = b3Cross( r, n );
			float wport = ai->invMass + aj->invMass + b3Dot( s, b3MulMV( ai->invInertia, s ) );
			float a_dim = K * h * h * wport;
			float p = a_dim * C / ( 1.0f + a_dim ) / ( h * wport );
			b3Vec3 P = b3MulSV( p, n );
			w->Jrecoil[i * RAFF_MAX_PORTS + k] = b3MulSV( -1.0f, P );
		}

		// Solve H * dq = -g (negate gradient to get RHS)
		for ( int a = 0; a < 6; a++ ) g[a] = -g[a];
		float dq[6];
		solve6x6( H, g, dq );

		// Update multipliers: lambda += K * C_post
		// C_post is the constraint violation AFTER applying the correction dq
		for ( int k = 0; k < ai->nPorts; k++ )
		{
			RAffPort* port = &w->ports[i * RAFF_MAX_PORTS + k];
			int j = port->neighIdx;
			if ( j < 0 || j >= w->nAtoms ) continue;

			// Compute post-solve tip position: pos + dx, r + cross(dtheta, r)
			b3Vec3 dx = { dq[0], dq[1], dq[2] };
			b3Vec3 dtheta = { dq[3], dq[4], dq[5] };
			b3Vec3 r = b3RotateVector( ai->quat, port->localDir );
			b3Vec3 r_new = b3Add( r, b3Cross( dtheta, r ) );
			b3Vec3 tip_new = b3Add( b3Add( ai->pos, dx ), r_new );
			b3Vec3 e_new = b3Sub( w->atoms[j].pos, tip_new );
			float C_post = b3Length( e_new );

			float K = port->kPort > 0.0f ? port->kPort : w->kBond * 0.5f;
			w->lambdaPort[i * RAFF_MAX_PORTS + k] += K * h2 * C_post;
		}

		// v_new = v + dq/h (ADD correction velocity, same as VBD)
		b3Vec3 dvel = { dq[0] / h, dq[1] / h, dq[2] / h };
		b3Vec3 domega = { dq[3] / h, dq[4] / h, dq[5] / h };
		b3Vec3 J = b3MulSV( ai->mass, dvel );
		b3Vec3 L = b3MulMV( Imat, domega );

		w->Jown[i] = J;
		w->Lown[i] = L;
	}
}

// --- Collision impulses (shared across all methods) ---

static void gatherCollisionImpulses( RAffWorld* w, float h )
{
	float K = w->kColl;
	float D = w->dColl;

	for ( int p = 0; p < w->nCollisionPairs; p++ )
	{
		int i = w->collisionPairs[p].i;
		int j = w->collisionPairs[p].j;
		float Rij = w->collisionPairs[p].Rij;

		RAffAtom* ai = &w->atoms[i];
		RAffAtom* aj = &w->atoms[j];

		b3Vec3 d = b3Sub( ai->pos, aj->pos );
		float r = b3Length( d );
		if ( r < 1e-8f ) continue;

		float delta = Rij - r;
		if ( delta <= 0.0f ) continue;

		b3Vec3 n = b3MulSV( 1.0f / r, d );
		float uc = b3Dot( n, b3Sub( ai->vel, aj->vel ) );
		float wc = ai->invMass + aj->invMass;

		float omega = sqrtf( K * wc );
		float zeta = ( omega > 1e-8f ) ? D * wc / ( 2.0f * omega ) : 0.0f;
		float a1 = 2.0f * zeta + h * omega;
		float a2 = h * omega * a1;
		float a3 = 1.0f / ( 1.0f + a2 );
		float biasRate = omega / a1;
		float massScale = a2 * a3;

		float bias = biasRate * delta;
		float effectiveMass = 1.0f / wc;
		float pImp = -effectiveMass * ( massScale * uc + bias );
		if ( pImp < 0.0f ) pImp = 0.0f;

		b3Vec3 P = b3MulSV( pImp, n );
		w->Jcoll[i] = b3Add( w->Jcoll[i], P );
		w->Jcoll[j] = b3Sub( w->Jcoll[j], P );
	}
}

// --- Apply phase: aggregate impulses and update atoms ---

static void applyAtomUpdate( RAffWorld* w, float h )
{
	for ( int i = 0; i < w->nAtoms; i++ )
	{
		RAffAtom* ai = &w->atoms[i];

		// Gather all linear impulses: own + incoming recoils + collisions
		b3Vec3 J = w->Jown[i];

		for ( int s = 0; s < w->incomingCount[i]; s++ )
		{
			int slot = w->incomingOffset[i] + s;
			int srcAtom = w->incomingAtoms[slot];
			int srcPort = w->incomingPorts[slot];
			J = b3Add( J, w->Jrecoil[srcAtom * RAFF_MAX_PORTS + srcPort] );
		}

		J = b3Add( J, w->Jcoll[i] );

		// Angular impulse (only from own ports; collisions have no torque)
		b3Vec3 L = w->Lown[i];

		// Update velocities
		ai->vel = b3Add( ai->vel, b3MulSV( ai->invMass, J ) );
		ai->omega = b3Add( ai->omega, b3MulMV( ai->invInertia, L ) );

		// Integrate position and rotation
		ai->pos = b3Add( ai->pos, b3MulSV( h, ai->vel ) );
		ai->quat = b3NormalizeQuat( b3IntegrateRotation( ai->quat, b3MulSV( h, ai->omega ) ) );
	}
}

// --- Hard microstep: dispatch by solver method ---

static void hardMicrostep( RAffWorld* w, float h )
{
	clearAccumulators( w );

	// Update collision pairs every microstep
	RAffWorld_UpdateCollisions( w );

	// Gather port impulses — method-specific
	switch ( w->solverMethod )
	{
		case RAFF_SOLVER_SI:   gatherPortImpulsesSI( w, h );   break;
		case RAFF_SOLVER_XPBD: gatherPortImpulsesXPBD( w, h ); break;
		case RAFF_SOLVER_PD:   gatherPortImpulsesPD( w, h );   break;
		case RAFF_SOLVER_VBD:  gatherPortImpulsesVBD( w, h );  break;
		case RAFF_SOLVER_AVBD: gatherPortImpulsesAVBD( w, h ); break;
		default:               gatherPortImpulsesSI( w, h );   break;
	}

	// Gather collision impulses (shared)
	gatherCollisionImpulses( w, h );

	// Apply: aggregate and update
	applyAtomUpdate( w, h );
}

// --- Macrostep ---

void RAffWorld_Step( RAffWorld* w )
{
	// Zero XPBD/AVBD Lagrange multipliers at macrostep start (warm start within microsteps only)
	if ( w->solverMethod == RAFF_SOLVER_XPBD || w->solverMethod == RAFF_SOLVER_AVBD )
		memset( w->lambdaPort, 0, (size_t)w->nAtoms * RAFF_MAX_PORTS * sizeof( float ) );

	// First half soft kick (using cached force from previous step)
	applySoftKick( w, 0.5f * w->H, w->fsoft, w->tsoft );

	float h = w->H / w->nsub;

	for ( int s = 0; s < w->nsub; s++ )
	{
		// No soft kick inside microsteps — the half kicks at top and bottom
		// already account for the full soft force impulse H.
		hardMicrostep( w, h );
	}

	// Evaluate new soft forces at final configuration
	RAffWorld_EvaluateSoftForces( w );

	// Second half soft kick with new force
	applySoftKick( w, 0.5f * w->H, w->fsoft, w->tsoft );
}

// --- Diagnostics ---

b3Vec3 RAffWorld_TotalLinearMomentum( RAffWorld* w )
{
	b3Vec3 P = b3Vec3_zero;
	for ( int i = 0; i < w->nAtoms; i++ )
		P = b3Add( P, b3MulSV( w->atoms[i].mass, w->atoms[i].vel ) );
	return P;
}

b3Vec3 RAffWorld_TotalAngularMomentum( RAffWorld* w )
{
	b3Vec3 L = b3Vec3_zero;
	for ( int i = 0; i < w->nAtoms; i++ )
	{
		RAffAtom* a = &w->atoms[i];
		L = b3Add( L, b3Cross( a->pos, b3MulSV( a->mass, a->vel ) ) );
		// Angular momentum from rotation: I * omega
		b3Vec3 Iw = b3MulMV( a->invInertia, a->omega );
		// Iw is invI*omega, so omega/invI = I*omega; but we store invI, not I.
		// For diagonal inertia, I*omega = omega / invI_diag.
		// For simplicity, use mass*radius^2 * omega as approximation.
		// For exact: L_rot = I * omega. Since we store invI, L_rot = omega * (1/invI).
		// We compute it as: L_rot_i = omega_i / invI_ii (component-wise for diagonal).
		// For general invI, we need the inverse: L_rot = inv(invI) * omega.
		// Use b3InvertMatrix to get I from invI.
		b3Matrix3 I = b3InvertMatrix( a->invInertia );
		L = b3Add( L, b3MulMV( I, a->omega ) );
	}
	return L;
}

void RAffWorld_LoadMol( RAffWorld* w, const b3MolFile* mol, b3Vec3 offset, float scale )
{
	int n = mol->atomCount;
	if ( n > w->nAtoms ) n = w->nAtoms;

	for ( int i = 0; i < n; i++ )
	{
		const b3MolAtom* ma = &mol->atoms[i];
		float radius = b3GetCovalentRadius( ma->element ) * scale;
		float mass = b3GetAtomicMass( ma->element );

		w->atoms[i].pos.x = offset.x + ma->x * scale;
		w->atoms[i].pos.y = offset.y + ma->y * scale;
		w->atoms[i].pos.z = offset.z + ma->z * scale;
		w->atoms[i].mass = mass;
		w->atoms[i].invMass = 1.0f / mass;
		w->atoms[i].radius = radius;
		w->atoms[i].nPorts = 0;
		w->atoms[i].quat = b3Quat_identity;
		w->atoms[i].vel = b3Vec3_zero;
		w->atoms[i].omega = b3Vec3_zero;

		float I = 0.4f * mass * radius * radius;
		b3Matrix3 Imat = { { I, 0, 0 }, { 0, I, 0 }, { 0, 0, I } };
		w->atoms[i].invInertia = b3InvertMatrix( Imat );
	}

	// Create reciprocal ports for each bond
	for ( int b = 0; b < mol->bondCount; b++ )
	{
		int a = mol->bonds[b].a;
		int bb = mol->bonds[b].b;
		if ( a < 0 || a >= n || bb < 0 || bb >= n ) continue;

		b3Vec3 ab = { mol->atoms[bb].x - mol->atoms[a].x, mol->atoms[bb].y - mol->atoms[a].y, mol->atoms[bb].z - mol->atoms[a].z };
		ab = b3MulSV( scale, ab );
		b3Vec3 ba = b3MulSV( -1.0f, ab );

		int k = w->atoms[a].nPorts;
		if ( k < RAFF_MAX_PORTS )
		{
			w->ports[a * RAFF_MAX_PORTS + k].neighIdx = bb;
			w->ports[a * RAFF_MAX_PORTS + k].localDir = ab;
			w->ports[a * RAFF_MAX_PORTS + k].kPort = 0.0f;
			w->ports[a * RAFF_MAX_PORTS + k].dPort = 0.0f;
			w->atoms[a].nPorts++;
		}

		k = w->atoms[bb].nPorts;
		if ( k < RAFF_MAX_PORTS )
		{
			w->ports[bb * RAFF_MAX_PORTS + k].neighIdx = a;
			w->ports[bb * RAFF_MAX_PORTS + k].localDir = ba;
			w->ports[bb * RAFF_MAX_PORTS + k].kPort = 0.0f;
			w->ports[bb * RAFF_MAX_PORTS + k].dPort = 0.0f;
			w->atoms[bb].nPorts++;
		}
	}

	RAffWorld_BuildAdjacency( w );

	for ( int i = 0; i < w->nAtoms; i++ )
	{
		w->fsoft[i] = b3Vec3_zero;
		w->tsoft[i] = b3Vec3_zero;
	}
}
