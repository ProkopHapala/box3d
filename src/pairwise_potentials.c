#include "pairwise_potentials.h"
#include <math.h>

void RAffPairPot_Setup( RAffPairPot* pot )
{
	float q = pot->qi_qj;
	float C = pot->C;
	int   n = pot->n;

	// Need attractive electrostatics (q < 0) and repulsion (C > 0) for a minimum.
	if ( q >= 0.0f || C <= 0.0f || n < 2 )
	{
		pot->r0 = 0.0f;
		pot->k_hard = 0.0f;
		pot->E_r0 = 0.0f;
		return;
	}

	// r0^(1-n) = -q / (n*C)  =>  r0 = (n*C / |q|)^(1/(n-1))
	float ratio = n * C / (-q); // -q = |q| since q < 0
	pot->r0 = powf( ratio, 1.0f / (float)( n - 1 ) );

	// k_hard = n*(n-1)*C / r0^(n+2)
	float r0_n2 = powf( pot->r0, (float)( n + 2 ) );
	pot->k_hard = (float)( n * ( n - 1 ) ) * C / r0_n2;

	// E(r0) = q/r0 + C/r0^n
	float r0_n = powf( pot->r0, (float)n );
	pot->E_r0 = q / pot->r0 + C / r0_n;
}

float RAffPairPot_SoftForce( const RAffPairPot* pot, float r )
{
	if ( r < 1e-12f ) return 0.0f;

	float q = pot->qi_qj;
	float C = pot->C;
	int   n = pot->n;

	// Full force: F = -dE/dr = q/r² + n*C/r^(n+1)
	// (positive = repulsive, pushing atoms apart)
	float r2 = r * r;
	float rn = powf( r, (float)n );
	float rn1 = rn * r; // r^(n+1)
	float F_full = q / r2 + (float)n * C / rn1;

	if ( r >= pot->r0 || pot->r0 <= 0.0f )
	{
		// Outside equilibrium: soft = full
		return F_full;
	}
	else
	{
		// Inside r0: subtract hard parabolic force
		// F_hard = k_hard * (r0 - r)  (repulsive)
		float F_hard = pot->k_hard * ( pot->r0 - r );
		return F_full - F_hard;
	}
}

float RAffPairPot_SoftEnergy( const RAffPairPot* pot, float r )
{
	if ( r < 1e-12f ) return 0.0f;

	float q = pot->qi_qj;
	float C = pot->C;
	int   n = pot->n;

	float rn = powf( r, (float)n );
	float E_full = q / r + C / rn;

	if ( r >= pot->r0 || pot->r0 <= 0.0f )
	{
		return E_full;
	}
	else
	{
		// Subtract hard parabolic energy: ½*k*(r0-r)²
		float dr = pot->r0 - r;
		float E_hard = 0.5f * pot->k_hard * dr * dr;
		return E_full - E_hard;
	}
}
