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

// --- Headless benchmark: sweep solver methods and timesteps ---

static float ComputeKE( RAffWorld* w )
{
	float ke = 0.0f;
	for ( int i = 0; i < w->nAtoms; i++ )
	{
		RAffAtom* a = &w->atoms[i];
		ke += 0.5f * a->mass * b3Dot( a->vel, a->vel );
		b3Matrix3 I = b3InvertMatrix( a->invInertia );
		b3Vec3 Iw = b3MulMV( I, a->omega );
		ke += 0.5f * b3Dot( a->omega, Iw );
	}
	return ke;
}

static float MaxBondError( RAffWorld* w )
{
	float maxErr = 0.0f;
	for ( int i = 0; i < w->nAtoms; i++ )
	{
		for ( int k = 0; k < w->atoms[i].nPorts; k++ )
		{
			RAffPort* port = &w->ports[i * RAFF_MAX_PORTS + k];
			int j = port->neighIdx;
			if ( j < 0 || j >= w->nAtoms ) continue;
			b3Vec3 r = b3RotateVector( w->atoms[i].quat, port->localDir );
			b3Vec3 tip = b3Add( w->atoms[i].pos, r );
			b3Vec3 e = b3Sub( w->atoms[j].pos, tip );
			float C = b3Length( e );
			if ( C > maxErr ) maxErr = C;
		}
	}
	return maxErr;
}

static bool HasNaN( RAffWorld* w )
{
	for ( int i = 0; i < w->nAtoms; i++ )
	{
		if ( isnan( w->atoms[i].pos.x ) || isnan( w->atoms[i].pos.y ) || isnan( w->atoms[i].pos.z ) ) return true;
		if ( isnan( w->atoms[i].vel.x ) || isnan( w->atoms[i].vel.y ) || isnan( w->atoms[i].vel.z ) ) return true;
		if ( isnan( w->atoms[i].omega.x ) || isnan( w->atoms[i].omega.y ) || isnan( w->atoms[i].omega.z ) ) return true;
	}
	return false;
}

static RAffWorld* CreateChainWorld( int n, float bondLen, float mass, float radius, float stretch )
{
	RAffWorld* w = RAffWorld_Create( n );
	for ( int i = 0; i < n; i++ )
	{
		w->atoms[i].pos = (b3Vec3){ i * bondLen * stretch, 0.05f * (float)i, 0.0f };
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
	for ( int i = 0; i < w->nAtoms; i++ )
	{
		w->fsoft[i] = b3Vec3_zero;
		w->tsoft[i] = b3Vec3_zero;
	}
	return w;
}

static const char* SolverName( RAffSolverMethod m )
{
	switch ( m )
	{
		case RAFF_SOLVER_SI:   return "SI";
		case RAFF_SOLVER_XPBD: return "XPBD";
		case RAFF_SOLVER_PD:   return "PD";
		case RAFF_SOLVER_VBD:  return "VBD";
		case RAFF_SOLVER_AVBD: return "AVBD";
		default:               return "???";
	}
}

int RigidAtomFF_Benchmark( void )
{
	int n = 6;
	float bondLen = 1.5f;
	float mass = 12.0f;
	float radius = 0.7f;
	float stretch = 1.3f; // bonds stretched 30%
	int nSteps = 200;

	RAffSolverMethod methods[] = { RAFF_SOLVER_SI, RAFF_SOLVER_XPBD, RAFF_SOLVER_PD, RAFF_SOLVER_VBD, RAFF_SOLVER_AVBD };
	int nMethods = (int)( sizeof( methods ) / sizeof( methods[0] ) );

	float Hs[] = { 0.001f, 0.002f, 0.005f, 0.01f, 0.02f, 0.05f, 0.1f, 0.2f };
	int nHs = (int)( sizeof( Hs ) / sizeof( Hs[0] ) );

	int nsubs[] = { 4, 8, 16, 32 };
	int nNsubs = (int)( sizeof( nsubs ) / sizeof( nsubs[0] ) );

	printf( "\n" );
	printf( "=== RigidAtomFF Stability Benchmark ===\n" );
	printf( "Chain: %d atoms, bondLen=%.2f, stretch=%.2f, mass=%.1f, radius=%.2f\n", n, bondLen, stretch, mass, radius );
	printf( "Steps per run: %d\n", nSteps );
	printf( "\n" );

	// Sweep: for each method x H x nsub, measure residual KE and max bond error
	printf( "%-6s  %8s  %6s  %12s  %12s  %12s  %s\n", "Method", "H", "nsub", "resid_KE", "max_bondErr", "max_pos", "Status" );
	printf( "--------------------------------------------------------------------------------------------------------\n" );

	for ( int mi = 0; mi < nMethods; mi++ )
	{
		for ( int hi = 0; hi < nHs; hi++ )
		{
			float H = Hs[hi];
			for ( int ni = 0; ni < nNsubs; ni++ )
			{
				int nsub = nsubs[ni];
				float h = H / nsub;

				// Skip if h*kBond*wport is absurdly large (would be unfair)
				// wport ~ 2/mass ~ 0.167, h*K*wport for K=5e5, h=0.05 -> 4167 — skip
				float K = 1e6f * 0.5f;
				float wport_approx = 2.0f / mass;
				float a_dim = K * h * h * wport_approx;
				if ( a_dim > 1e4f ) continue;

				RAffWorld* w = CreateChainWorld( n, bondLen, mass, radius, stretch );
				w->solverMethod = methods[mi];
				w->H = H;
				w->nsub = nsub;
				w->kBond = 1e6f;
				w->dBond = 1e4f;
				w->kColl = 5e6f;
				w->dColl = 5e4f;
				w->softCutoff = 0.0f;

				for ( int step = 0; step < nSteps; step++ )
				{
					RAffWorld_Step( w );
					if ( HasNaN( w ) ) break;
				}

				float ke = HasNaN( w ) ? 1e30f : ComputeKE( w );
				float bondErr = HasNaN( w ) ? 1e30f : MaxBondError( w );
				float maxPos = 0.0f;
				for ( int i = 0; i < n; i++ )
				{
					float r2 = b3Dot( w->atoms[i].pos, w->atoms[i].pos );
					if ( r2 > maxPos ) maxPos = r2;
				}
				maxPos = sqrtf( maxPos );

				const char* status;
				if ( HasNaN( w ) ) status = "NaN!";
				else if ( maxPos > 1e4f ) status = "EXPLODED";
				else if ( ke > 1e3f ) status = "AGITATED";
				else if ( ke > 1.0f ) status = "warm";
				else if ( bondErr > 0.5f ) status = "poor_convergence";
				else status = "RELAXED";

				printf( "%-6s  %8.4f  %6d  %12.4e  %12.6f  %12.2f  %s\n",
					SolverName( methods[mi] ), H, nsub, ke, bondErr, maxPos, status );

				RAffWorld_Destroy( w );
			}
		}
		printf( "\n" );
	}

	// Summary: find max stable H per method
	printf( "=== Summary: max stable H (resid_KE < 1.0, no NaN) ===\n" );
	for ( int mi = 0; mi < nMethods; mi++ )
	{
		float bestH = 0.0f;
		int bestNsub = 0;
		for ( int hi = nHs - 1; hi >= 0; hi-- )
		{
			float H = Hs[hi];
			for ( int ni = 0; ni < nNsubs; ni++ )
			{
				int nsub = nsubs[ni];
				float h = H / nsub;
				float K = 1e6f * 0.5f;
				float wport_approx = 2.0f / mass;
				float a_dim = K * h * h * wport_approx;
				if ( a_dim > 1e4f ) continue;

				RAffWorld* w = CreateChainWorld( n, bondLen, mass, radius, stretch );
				w->solverMethod = methods[mi];
				w->H = H;
				w->nsub = nsub;
				w->kBond = 1e6f;
				w->dBond = 1e4f;
				w->kColl = 5e6f;
				w->dColl = 5e4f;
				w->softCutoff = 0.0f;

				for ( int step = 0; step < nSteps; step++ )
				{
					RAffWorld_Step( w );
					if ( HasNaN( w ) ) break;
				}

				if ( !HasNaN( w ) && ComputeKE( w ) < 1.0f )
				{
					if ( H > bestH )
					{
						bestH = H;
						bestNsub = nsub;
					}
				}
				RAffWorld_Destroy( w );
			}
		}
		printf( "  %-6s: max stable H = %.4f (nsub=%d, h=%.5f)\n",
			SolverName( methods[mi] ), bestH, bestNsub, bestH / bestNsub );
	}
	printf( "\n" );

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
