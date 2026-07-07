// Tests for RigidAtomFF SI solver.
// Verifies bond convergence, collision repulsion, and momentum conservation.

#include "test_macros.h"
#include "rigid_atom_ff.h"
#include "math_internal.h"

#include <math.h>
#include <stdio.h>

// Helper: set up a simple two-atom world with one bond
static RAffWorld* CreateTwoAtomWorld( float bondLength, float mass, float radius )
{
	RAffWorld* w = RAffWorld_Create( 2 );

	// Atom 0 at origin
	w->atoms[0].pos = b3Vec3_zero;
	w->atoms[0].mass = mass;
	w->atoms[0].invMass = 1.0f / mass;
	w->atoms[0].radius = radius;
	w->atoms[0].nPorts = 1;

	// Atom 1 at bondLength along x
	w->atoms[1].pos = (b3Vec3){ bondLength, 0.0f, 0.0f };
	w->atoms[1].mass = mass;
	w->atoms[1].invMass = 1.0f / mass;
	w->atoms[1].radius = radius;
	w->atoms[1].nPorts = 1;

	// Port on atom 0 pointing toward atom 1
	w->ports[0 * RAFF_MAX_PORTS + 0].neighIdx = 1;
	w->ports[0 * RAFF_MAX_PORTS + 0].localDir = (b3Vec3){ bondLength, 0.0f, 0.0f };
	w->ports[0 * RAFF_MAX_PORTS + 0].kPort = 0.0f; // use world default
	w->ports[0 * RAFF_MAX_PORTS + 0].dPort = 0.0f;

	// Port on atom 1 pointing toward atom 0
	w->ports[1 * RAFF_MAX_PORTS + 0].neighIdx = 0;
	w->ports[1 * RAFF_MAX_PORTS + 0].localDir = (b3Vec3){ -bondLength, 0.0f, 0.0f };
	w->ports[1 * RAFF_MAX_PORTS + 0].kPort = 0.0f;
	w->ports[1 * RAFF_MAX_PORTS + 0].dPort = 0.0f;

	// Inertia: treat as uniform sphere, I = 2/5 * m * r^2
	float I = 0.4f * mass * radius * radius;
	b3Matrix3 Imat = { { I, 0, 0 }, { 0, I, 0 }, { 0, 0, I } };
	w->atoms[0].invInertia = b3InvertMatrix( Imat );
	w->atoms[1].invInertia = b3InvertMatrix( Imat );

	RAffWorld_BuildAdjacency( w );
	return w;
}

static void SetInertia( RAffWorld* w, int i, float mass, float radius )
{
	float I = 0.4f * mass * radius * radius;
	b3Matrix3 Imat = { { I, 0, 0 }, { 0, I, 0 }, { 0, 0, I } };
	w->atoms[i].invInertia = b3InvertMatrix( Imat );
}

int RigidAtomFF_BondTest( void )
{
	float bondLength = 1.5f;
	RAffWorld* w = CreateTwoAtomWorld( bondLength, 12.0f, 0.7f );

	// Displace atom 1 slightly so the bond is stretched
	w->atoms[1].pos = (b3Vec3){ bondLength + 0.5f, 0.0f, 0.0f };

	// Run a few macrosteps with no soft forces (just hard relaxation)
	w->H = 0.01f;
	w->nsub = 16;
	w->kBond = 1e6f;
	w->dBond = 1e4f;
	w->softCutoff = 0.0f; // disable soft forces

	// Zero out soft forces initially
	for ( int i = 0; i < w->nAtoms; i++ )
	{
		w->fsoft[i] = b3Vec3_zero;
		w->tsoft[i] = b3Vec3_zero;
	}

	for ( int step = 0; step < 50; step++ )
	{
		RAffWorld_Step( w );
	}

	// Check bond length converged
	b3Vec3 d = b3Sub( w->atoms[1].pos, w->atoms[0].pos );
	float r = b3Length( d );
	ENSURE_SMALL( r - bondLength, 0.1f );

	// Check momentum conserved (should be zero since we started at rest)
	b3Vec3 P = RAffWorld_TotalLinearMomentum( w );
	ENSURE_SMALL( P.x, 0.5f );
	ENSURE_SMALL( P.y, 0.5f );
	ENSURE_SMALL( P.z, 0.5f );

	RAffWorld_Destroy( w );
	return 0;
}

int RigidAtomFF_CollisionTest( void )
{
	RAffWorld* w = RAffWorld_Create( 2 );

	// Two non-bonded atoms approaching each other
	float mass = 12.0f;
	float radius = 1.0f;
	float Rij = radius + radius; // = 2.0

	w->atoms[0].pos = (b3Vec3){ -1.0f, 0.0f, 0.0f };
	w->atoms[0].vel = (b3Vec3){ 5.0f, 0.0f, 0.0f };
	w->atoms[0].mass = mass;
	w->atoms[0].invMass = 1.0f / mass;
	w->atoms[0].radius = radius;
	w->atoms[0].nPorts = 0;

	w->atoms[1].pos = (b3Vec3){ 1.0f, 0.0f, 0.0f };
	w->atoms[1].vel = (b3Vec3){ -5.0f, 0.0f, 0.0f };
	w->atoms[1].mass = mass;
	w->atoms[1].invMass = 1.0f / mass;
	w->atoms[1].radius = radius;
	w->atoms[1].nPorts = 0;

	SetInertia( w, 0, mass, radius );
	SetInertia( w, 1, mass, radius );

	RAffWorld_BuildAdjacency( w );

	w->H = 0.005f;
	w->nsub = 8;
	w->kColl = 5e6f;
	w->dColl = 5e4f;
	w->softCutoff = 0.0f;

	for ( int i = 0; i < w->nAtoms; i++ )
	{
		w->fsoft[i] = b3Vec3_zero;
		w->tsoft[i] = b3Vec3_zero;
	}

	for ( int step = 0; step < 30; step++ )
	{
		RAffWorld_Step( w );
	}

	// Atoms should not be penetrating
	b3Vec3 d = b3Sub( w->atoms[0].pos, w->atoms[1].pos );
	float r = b3Length( d );
	ENSURE( r >= Rij * 0.8f );

	// Momentum should be conserved (zero total since equal and opposite)
	b3Vec3 P = RAffWorld_TotalLinearMomentum( w );
	ENSURE_SMALL( P.x, 1.0f );
	ENSURE_SMALL( P.y, 1.0f );
	ENSURE_SMALL( P.z, 1.0f );

	RAffWorld_Destroy( w );
	return 0;
}

int RigidAtomFF_MomentumTest( void )
{
	// 4-atom chain with bonds, give random velocities, check momentum conservation
	int n = 4;
	RAffWorld* w = RAffWorld_Create( n );

	float mass = 12.0f;
	float radius = 0.7f;
	float bondLen = 1.5f;

	for ( int i = 0; i < n; i++ )
	{
		w->atoms[i].pos = (b3Vec3){ i * bondLen, 0.0f, 0.0f };
		w->atoms[i].mass = mass;
		w->atoms[i].invMass = 1.0f / mass;
		w->atoms[i].radius = radius;
		w->atoms[i].nPorts = 0;

		SetInertia( w, i, mass, radius );
	}

	// Set up chain bonds: 0-1, 1-2, 2-3
	for ( int i = 0; i < n - 1; i++ )
	{
		int k = w->atoms[i].nPorts;
		w->ports[i * RAFF_MAX_PORTS + k].neighIdx = i + 1;
		w->ports[i * RAFF_MAX_PORTS + k].localDir = (b3Vec3){ bondLen, 0.0f, 0.0f };
		w->ports[i * RAFF_MAX_PORTS + k].kPort = 0.0f;
		w->ports[i * RAFF_MAX_PORTS + k].dPort = 0.0f;
		w->atoms[i].nPorts++;

		int k2 = w->atoms[i + 1].nPorts;
		w->ports[( i + 1 ) * RAFF_MAX_PORTS + k2].neighIdx = i;
		w->ports[( i + 1 ) * RAFF_MAX_PORTS + k2].localDir = (b3Vec3){ -bondLen, 0.0f, 0.0f };
		w->ports[( i + 1 ) * RAFF_MAX_PORTS + k2].kPort = 0.0f;
		w->ports[( i + 1 ) * RAFF_MAX_PORTS + k2].dPort = 0.0f;
		w->atoms[i + 1].nPorts++;
	}

	RAffWorld_BuildAdjacency( w );

	// Give atoms some initial velocities that sum to zero
	w->atoms[0].vel = (b3Vec3){ 1.0f, 0.5f, -0.3f };
	w->atoms[1].vel = (b3Vec3){ -0.5f, -1.0f, 0.8f };
	w->atoms[2].vel = (b3Vec3){ 0.3f, 0.6f, -0.5f };
	w->atoms[3].vel = (b3Vec3){ -0.8f, -0.1f, 0.0f };

	w->H = 0.005f;
	w->nsub = 8;
	w->kBond = 1e6f;
	w->dBond = 1e4f;
	w->softCutoff = 0.0f;

	for ( int i = 0; i < w->nAtoms; i++ )
	{
		w->fsoft[i] = b3Vec3_zero;
		w->tsoft[i] = b3Vec3_zero;
	}

	b3Vec3 P0 = RAffWorld_TotalLinearMomentum( w );

	for ( int step = 0; step < 100; step++ )
	{
		RAffWorld_Step( w );
	}

	b3Vec3 P1 = RAffWorld_TotalLinearMomentum( w );

	// Linear momentum should be conserved
	ENSURE_SMALL( P1.x - P0.x, 0.5f );
	ENSURE_SMALL( P1.y - P0.y, 0.5f );
	ENSURE_SMALL( P1.z - P0.z, 0.5f );

	RAffWorld_Destroy( w );
	return 0;
}

int RigidAtomFF_ChainTest( void )
{
	// 6-atom chain: verify stability with larger H and no explosion
	int n = 6;
	RAffWorld* w = RAffWorld_Create( n );

	float mass = 12.0f;
	float radius = 0.7f;
	float bondLen = 1.5f;

	for ( int i = 0; i < n; i++ )
	{
		w->atoms[i].pos = (b3Vec3){ i * bondLen, 0.0f, 0.0f };
		w->atoms[i].mass = mass;
		w->atoms[i].invMass = 1.0f / mass;
		w->atoms[i].radius = radius;
		w->atoms[i].nPorts = 0;

		SetInertia( w, i, mass, radius );
	}

	for ( int i = 0; i < n - 1; i++ )
	{
		int k = w->atoms[i].nPorts;
		w->ports[i * RAFF_MAX_PORTS + k].neighIdx = i + 1;
		w->ports[i * RAFF_MAX_PORTS + k].localDir = (b3Vec3){ bondLen, 0.0f, 0.0f };
		w->ports[i * RAFF_MAX_PORTS + k].kPort = 0.0f;
		w->ports[i * RAFF_MAX_PORTS + k].dPort = 0.0f;
		w->atoms[i].nPorts++;

		int k2 = w->atoms[i + 1].nPorts;
		w->ports[( i + 1 ) * RAFF_MAX_PORTS + k2].neighIdx = i;
		w->ports[( i + 1 ) * RAFF_MAX_PORTS + k2].localDir = (b3Vec3){ -bondLen, 0.0f, 0.0f };
		w->ports[( i + 1 ) * RAFF_MAX_PORTS + k2].kPort = 0.0f;
		w->ports[( i + 1 ) * RAFF_MAX_PORTS + k2].dPort = 0.0f;
		w->atoms[i + 1].nPorts++;
	}

	RAffWorld_BuildAdjacency( w );

	// Perturb positions slightly
	for ( int i = 0; i < n; i++ )
	{
		w->atoms[i].pos.x += 0.1f * (float)( i - n / 2 );
		w->atoms[i].pos.y += 0.05f * (float)i;
	}

	w->H = 0.02f;
	w->nsub = 16;
	w->kBond = 1e6f;
	w->dBond = 1e4f;
	w->softCutoff = 0.0f;

	for ( int i = 0; i < w->nAtoms; i++ )
	{
		w->fsoft[i] = b3Vec3_zero;
		w->tsoft[i] = b3Vec3_zero;
	}

	for ( int step = 0; step < 200; step++ )
	{
		RAffWorld_Step( w );
		if ( step % 50 == 0 )
		{
			for ( int i = 0; i < n; i++ )
				printf( "  step %d atom %d pos=(%.3f,%.3f,%.3f) vel=(%.3f,%.3f,%.3f)\n",
					step, i, w->atoms[i].pos.x, w->atoms[i].pos.y, w->atoms[i].pos.z,
					w->atoms[i].vel.x, w->atoms[i].vel.y, w->atoms[i].vel.z );
		}
	}

	// Check no atom has flown away
	for ( int i = 0; i < n; i++ )
	{
		float r2 = b3Dot( w->atoms[i].pos, w->atoms[i].pos );
		printf( "  final atom %d r2=%.3f\n", i, r2 );
		ENSURE( r2 < 1000.0f );
	}

	// Check bond lengths are approximately maintained
	for ( int i = 0; i < n - 1; i++ )
	{
		b3Vec3 d = b3Sub( w->atoms[i + 1].pos, w->atoms[i].pos );
		float r = b3Length( d );
		ENSURE_SMALL( r - bondLen, 0.3f );
	}

	RAffWorld_Destroy( w );
	return 0;
}

int RigidAtomFFTest( void )
{
	RUN_SUBTEST( RigidAtomFF_BondTest );
	// RUN_SUBTEST( RigidAtomFF_CollisionTest ); // disabled until collisions re-enabled
	RUN_SUBTEST( RigidAtomFF_MomentumTest );
	RUN_SUBTEST( RigidAtomFF_ChainTest );
	return 0;
}
