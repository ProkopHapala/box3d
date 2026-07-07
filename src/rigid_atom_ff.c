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

static void solvePortImpulsesSI( RAffWorld* w, float h )
{
	// Sequential Impulses: process each directed port, update velocities immediately.
	// Each port uses K/2 to avoid doubling bond energy (reciprocal ports).
	for ( int i = 0; i < w->nAtoms; i++ )
	{
		RAffAtom* ai = &w->atoms[i];
		for ( int k = 0; k < ai->nPorts; k++ )
		{
			RAffPort* port = &w->ports[i * RAFF_MAX_PORTS + k];
			int j = port->neighIdx;
			if ( j < 0 || j >= w->nAtoms ) continue;

			RAffAtom* aj = &w->atoms[j];

			b3Vec3 r = b3RotateVector( ai->quat, port->localDir );
			b3Vec3 tip = b3Add( ai->pos, r );
			// Design doc convention: e = x_j - tip, n points from tip toward neighbor
			b3Vec3 e = b3Sub( aj->pos, tip );
			float C = b3Length( e );
			if ( C < 1e-8f ) continue;
			b3Vec3 n = b3MulSV( 1.0f / C, e );

			// u = n · (v_j - tipVel) = dC/dt (positive when separating)
			b3Vec3 tipVel = b3Add( ai->vel, b3Cross( ai->omega, r ) );
			float u = b3Dot( n, b3Sub( aj->vel, tipVel ) );

			b3Vec3 s = b3Cross( r, n );
			float wport = ai->invMass + aj->invMass + b3Dot( s, b3MulMV( ai->invInertia, s ) );

			// Use K/2 per directed port (reciprocal ports double the energy)
			float K = port->kPort > 0.0f ? port->kPort : w->kBond * 0.5f;
			float D = port->dPort > 0.0f ? port->dPort : w->dBond * 0.5f;

			// Local implicit (backward Euler) impulse — NO minus sign:
			// p = (hK*C + h*(D+hK)*u) / (1 + (hD+h²K)*w)
			float hK = h * K;
			float hD = h * D;
			float h2K = h * hK;
			float p = ( hK * C + h * ( D + hK ) * u ) / ( 1.0f + ( hD + h2K ) * wport );

			b3Vec3 P = b3MulSV( p, n );

			// Apply +P to atom i (at port tip), -P to atom j (at center)
			ai->vel = b3Add( ai->vel, b3MulSV( ai->invMass, P ) );
			ai->omega = b3Add( ai->omega, b3MulMV( ai->invInertia, b3Cross( r, P ) ) );
			aj->vel = b3Sub( aj->vel, b3MulSV( aj->invMass, P ) );
		}
	}
}

static void solveCollisionImpulsesSI( RAffWorld* w, float h )
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
		ai->vel = b3Add( ai->vel, b3MulSV( ai->invMass, P ) );
		aj->vel = b3Sub( aj->vel, b3MulSV( aj->invMass, P ) );
	}
}

static void integratePositions( RAffWorld* w, float h )
{
	for ( int i = 0; i < w->nAtoms; i++ )
	{
		RAffAtom* ai = &w->atoms[i];
		ai->pos = b3Add( ai->pos, b3MulSV( h, ai->vel ) );
		ai->quat = b3NormalizeQuat( b3IntegrateRotation( ai->quat, b3MulSV( h, ai->omega ) ) );
	}
}

static void hardMicrostep( RAffWorld* w, float h )
{
	// Single SI sweep per microstep (design doc: one gather-apply per microstep)
	solvePortImpulsesSI( w, h );
	// Collisions disabled until bond interaction is stable
	// solveCollisionImpulsesSI( w, h );
	integratePositions( w, h );
}

// --- Macrostep ---

void RAffWorld_Step( RAffWorld* w )
{
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
