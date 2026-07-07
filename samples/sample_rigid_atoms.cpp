// SPDX-FileCopyrightText: 2025 Erin Catto
// SPDX-License-Identifier: MIT

// RigidAtomFF visualization sample.
// Uses the RAffWorld solver (not Box3D's joint solver) to simulate
// atom chains with port constraints. Atoms are drawn as spheres.
// Mouse picking grabs atoms. A ground plane provides substrate collision.

#include "gfx/draw.h"
#include "gfx/keycodes.h"
#include "imgui.h"
#include "implot.h"
#include "sample.h"

#include "rigid_atom_ff.h"

#include <math.h>
#include <stdio.h>

static const int KE_HISTORY_SIZE = 512;

static float Clampf( float v, float lo, float hi ) { return v < lo ? lo : ( v > hi ? hi : v ); }

// Ray-sphere intersection. Returns t >= 0 if hit, -1 otherwise.
static float RaySphere( b3Vec3 origin, b3Vec3 dir, b3Vec3 center, float radius )
{
	b3Vec3 oc = b3Sub( origin, center );
	float b = b3Dot( oc, dir );
	float c = b3Dot( oc, oc ) - radius * radius;
	float disc = b * b - c;
	if ( disc < 0.0f ) return -1.0f;
	float t = -b - sqrtf( disc );
	if ( t < 0.0f ) t = -b + sqrtf( disc );
	return t;
}

class RigidAtomsSample : public Sample
{
public:
	explicit RigidAtomsSample( SampleContext* context )
		: Sample( context )
	{
		if ( context->restart == false )
		{
			m_camera->SetView( -35.0f, 25.0f, 40.0f, { 0.0f, 5.0f, 0.0f } );
		}

		AddGroundBox( 40.0f );

		m_pickedAtom = -1;
		m_grabActive = false;
		m_grabDepth = 0.0f;

		m_gravity = 9.81f;
		m_restitution = 0.3f;
		m_friction = 0.98f;
		m_grabStiffness = 200.0f;
		m_reportInterval = 60; // frames between terminal reports (0 = off)
		m_keHistoryCount = 0;
		m_keHistoryWrite = 0;
		m_nanDetected = false;

		// Try to load nNonan molecule
		m_mol = {};
		m_molLoaded = false;
		const char* paths[] = {
			"data/mol/nNonan.mol",
			"../data/mol/nNonan.mol",
			"../../data/mol/nNonan.mol",
			"../../../data/mol/nNonan.mol",
		};
		for ( int i = 0; i < 4; i++ )
		{
			if ( b3LoadMolFile( paths[i], &m_mol ) )
			{
				m_molLoaded = true;
				break;
			}
		}

		SpawnMolecule();
	}

	~RigidAtomsSample() override
	{
		if ( m_raff ) RAffWorld_Destroy( m_raff );
	}

	void SpawnMolecule()
	{
		if ( m_raff ) RAffWorld_Destroy( m_raff );

		if ( !m_molLoaded )
		{
			DrawTextLine( "ERROR: could not load nNonan.mol" );
			return;
		}

		int n = m_mol.atomCount;
		m_raff = RAffWorld_Create( n );
		m_raff->H = 0.02f;
		m_raff->nsub = 16;
		m_raff->kBond = 1e6f;
		m_raff->dBond = 1e4f;
		m_raff->softCutoff = 0.0f;

		// Load molecule: scale Angstroms to Box3D units (1 A ≈ 0.1 nm, use 1.0 for direct)
		// Place molecule above ground
		float scale = 1.0f;
		b3Vec3 offset = { 0.0f, 8.0f, 0.0f };
		RAffWorld_LoadMol( m_raff, &m_mol, offset, scale );
	}

	void Step() override
	{
		// Step Box3D world (for ground rendering only — no dynamic bodies)
		Sample::Step();

		if ( m_context->pause == false || m_context->singleStep > 0 )
		{
			// Apply gravity to atoms
			for ( int i = 0; i < m_raff->nAtoms; i++ )
			{
				m_raff->atoms[i].vel.y -= m_gravity * m_raff->H;
			}

			RAffWorld_Step( m_raff );

			// Ground collision: simple position correction + velocity bounce
			float groundY = 0.0f; // ground top surface
			for ( int i = 0; i < m_raff->nAtoms; i++ )
			{
				RAffAtom* a = &m_raff->atoms[i];
				float penetration = ( a->pos.y - a->radius ) - groundY;
				if ( penetration < 0.0f )
				{
					a->pos.y = groundY + a->radius;
					if ( a->vel.y < 0.0f )
						a->vel.y = -a->vel.y * m_restitution;
					a->vel.x *= m_friction;
					a->vel.z *= m_friction;
				}
			}

			// Apply mouse grab force
			if ( m_grabActive && m_pickedAtom >= 0 )
			{
				RAffAtom* a = &m_raff->atoms[m_pickedAtom];
				b3Vec3 target = m_grabTarget;
				b3Vec3 delta = b3Sub( target, a->pos );
				b3Vec3 force = b3MulSV( m_grabStiffness, delta );
				a->vel = b3Add( a->vel, b3MulSV( m_raff->H * a->invMass, force ) );
			}

			m_context->singleStep = b3MaxInt( 0, m_context->singleStep - 1 );
		}

		// Compute diagnostics
		float ke = 0.0f;
		float maxVel = 0.0f;
		bool hasNaN = false;
		if ( m_raff )
		{
			for ( int i = 0; i < m_raff->nAtoms; i++ )
			{
				RAffAtom* a = &m_raff->atoms[i];
				float v2 = b3Dot( a->vel, a->vel );
				float w2 = b3Dot( a->omega, a->omega );
				ke += 0.5f * a->mass * v2 + 0.5f * ( 1.0f / ( a->invMass + 1e-12f ) ) * w2 * 0.4f * a->radius * a->radius;
				if ( v2 > maxVel ) maxVel = v2;
				if ( !isfinite( a->pos.x ) || !isfinite( a->pos.y ) || !isfinite( a->pos.z ) ||
					 !isfinite( a->vel.x ) || !isfinite( a->vel.y ) || !isfinite( a->vel.z ) )
					hasNaN = true;
			}
			maxVel = sqrtf( maxVel );
		}

		if ( hasNaN && !m_nanDetected )
		{
			m_nanDetected = true;
			printf( "[RigidAtomFF] NaN detected at frame %d!\n", m_stepCount );
		}

		// Store KE history
		if ( m_raff && m_keHistoryCount < KE_HISTORY_SIZE )
		{
			m_keHistory[m_keHistoryWrite] = ke;
			m_keHistoryWrite = ( m_keHistoryWrite + 1 ) % KE_HISTORY_SIZE;
			m_keHistoryCount++;
		}
		else if ( m_raff )
		{
			m_keHistory[m_keHistoryWrite] = ke;
			m_keHistoryWrite = ( m_keHistoryWrite + 1 ) % KE_HISTORY_SIZE;
		}

		// Terminal report
		if ( m_raff && m_reportInterval > 0 && m_stepCount % m_reportInterval == 0 )
		{
			printf( "[RigidAtomFF] frame=%d  KE=%.6e  maxVel=%.4f  H=%.4f  nsub=%d  kBond=%.0f  dBond=%.0f%s\n",
				m_stepCount, ke, maxVel, m_raff->H, m_raff->nsub, m_raff->kBond, m_raff->dBond,
				hasNaN ? "  ** NaN **" : "" );
		}

		DrawTextLine( "RigidAtomFF: nNonan (C9H20) with ground collision" );
		DrawTextLine( "Ctrl+Click to grab an atom" );
		if ( m_raff )
		{
			DrawTextLine( "Atoms: %d  Bonds: %d  H=%.3f  nsub=%d", m_raff->nAtoms, m_mol.bondCount, m_raff->H, m_raff->nsub );
			DrawTextLine( "KE=%.3e  maxVel=%.3f%s", ke, maxVel, hasNaN ? "  ** NaN **" : "" );
		}
		else
			DrawTextLine( "ERROR: molecule not loaded" );
	}

	void Render() override
	{
		if ( !m_raff ) return;

		// Draw atoms as spheres, colored by element
		for ( int i = 0; i < m_raff->nAtoms; i++ )
		{
			RAffAtom* a = &m_raff->atoms[i];
			b3WorldTransform xf;
			xf.p = { a->pos.x, a->pos.y, a->pos.z };
			xf.q = a->quat;

			Vec4 color;
			if ( i == m_pickedAtom && m_grabActive )
				color = { 1.0f, 0.5f, 0.2f, 1.0f }; // orange = grabbed
			else if ( i == m_pickedAtom )
				color = { 0.2f, 1.0f, 0.2f, 1.0f }; // green = selected
			else if ( i < m_mol.atomCount && m_mol.atoms[i].element[0] == 'H' )
				color = { 0.9f, 0.9f, 0.9f, 1.0f }; // white = hydrogen
			else
				color = { 0.3f, 0.5f, 0.9f, 1.0f }; // blue = carbon/other

			DrawSphere( xf, a->radius, color );
		}

		// Draw bonds as lines
		Vec4 bondColor = { 0.8f, 0.8f, 0.8f, 1.0f };
		for ( int i = 0; i < m_raff->nAtoms; i++ )
		{
			RAffAtom* ai = &m_raff->atoms[i];
			for ( int k = 0; k < ai->nPorts; k++ )
			{
				RAffPort* port = &m_raff->ports[i * RAFF_MAX_PORTS + k];
				int j = port->neighIdx;
				if ( j < 0 || j >= m_raff->nAtoms ) continue;
				if ( i > j ) continue; // draw each bond once

				RAffAtom* aj = &m_raff->atoms[j];
				b3Vec3 r = b3RotateVector( ai->quat, port->localDir );
				b3Vec3 tip = b3Add( ai->pos, r );
				DrawLine( { tip.x, tip.y, tip.z }, { aj->pos.x, aj->pos.y, aj->pos.z }, bondColor );
			}
		}
	}

	bool DrawControls() override
	{
		if ( !m_raff ) return false;

		ImGui::PushItemWidth( 6.0f * ImGui::GetFontSize() );

		ImGui::Text( "RigidAtomFF Solver" );
		ImGui::Text( "Molecule: nNonan (C9H20)" );
		ImGui::Text( "Atoms: %d  Bonds: %d", m_mol.atomCount, m_mol.bondCount );
		ImGui::Separator();

		ImGui::SliderFloat( "H (macrostep)", &m_raff->H, 0.001f, 0.05f, "%4.3f" );
		ImGui::SliderInt( "nsub", &m_raff->nsub, 1, 128 );
		ImGui::SliderFloat( "kBond", &m_raff->kBond, 1e3f, 1e8f, "%.0f" );
		ImGui::SliderFloat( "dBond", &m_raff->dBond, 1e1f, 1e6f, "%.0f" );

		ImGui::Separator();
		ImGui::SliderFloat( "gravity", &m_gravity, 0.0f, 30.0f, "%.2f" );
		ImGui::SliderFloat( "restitution", &m_restitution, 0.0f, 1.0f, "%.2f" );
		ImGui::SliderFloat( "friction", &m_friction, 0.8f, 1.0f, "%.4f" );
		ImGui::SliderFloat( "grab stiffness", &m_grabStiffness, 10.0f, 2000.0f, "%.0f" );
		ImGui::SliderInt( "report interval", &m_reportInterval, 0, 300 );

		ImGui::Separator();
		if ( ImGui::Button( "Respawn" ) )
		{
			SpawnMolecule();
			m_keHistoryCount = 0;
			m_keHistoryWrite = 0;
			m_nanDetected = false;
		}

		ImGui::PopItemWidth();

		// KE history plot
		if ( m_keHistoryCount > 1 )
		{
			ImGui::Separator();
			ImVec2 plotSize = ImGui::GetContentRegionAvail();
			plotSize.y = 150;
			if ( ImPlot::BeginPlot( "Kinetic Energy", plotSize, ImPlotFlags_NoTitle ) )
			{
				ImPlot::SetupAxes( "frame", "KE" );
				// Plot the ring buffer linearly
				if ( m_keHistoryCount < KE_HISTORY_SIZE )
				{
					ImPlot::PlotLine( "KE", m_keHistory, m_keHistoryCount );
				}
				else
				{
					// Full ring: plot in order
					float ordered[KE_HISTORY_SIZE];
					for ( int i = 0; i < KE_HISTORY_SIZE; i++ )
						ordered[i] = m_keHistory[( m_keHistoryWrite + i ) % KE_HISTORY_SIZE];
					ImPlot::PlotLine( "KE", ordered, KE_HISTORY_SIZE );
				}
				ImPlot::EndPlot();
			}
		}

		return true;
	}

	void MouseDown( b3Vec2 p, int button, int modifiers ) override
	{
		if ( !m_raff ) { Sample::MouseDown( p, button, modifiers ); return; }

		if ( button == 0 && ( modifiers & MOD_CTRL ) )
		{
			// Ctrl+Click: grab nearest atom under cursor
			PickRay pickRay = m_camera->BuildPickRay( p.x, p.y );
			b3Vec3 origin = { (float)pickRay.origin.x, (float)pickRay.origin.y, (float)pickRay.origin.z };
			b3Vec3 dir = b3Normalize( pickRay.translation );

			int bestAtom = -1;
			float bestT = 1e30f;
			for ( int i = 0; i < m_raff->nAtoms; i++ )
			{
				float t = RaySphere( origin, dir, m_raff->atoms[i].pos, m_raff->atoms[i].radius );
				if ( t >= 0.0f && t < bestT )
				{
					bestT = t;
					bestAtom = i;
				}
			}

			if ( bestAtom >= 0 )
			{
				m_pickedAtom = bestAtom;
				m_grabActive = true;
				m_grabDepth = bestT;
				b3Vec3 hitPoint = b3Add( origin, b3MulSV( bestT, dir ) );
				m_grabTarget = hitPoint;
			}
		}
		else
		{
			// Fall through to default behavior (selection, etc.)
			Sample::MouseDown( p, button, modifiers );
		}
	}

	void MouseUp( b3Vec2 p, int button ) override
	{
		if ( m_grabActive )
		{
			m_grabActive = false;
		}
		else
		{
			Sample::MouseUp( p, button );
		}
	}

	void MouseMove( b3Vec2 p ) override
	{
		if ( m_grabActive && m_pickedAtom >= 0 )
		{
			PickRay pickRay = m_camera->BuildPickRay( p.x, p.y );
			b3Vec3 origin = { (float)pickRay.origin.x, (float)pickRay.origin.y, (float)pickRay.origin.z };
			b3Vec3 dir = b3Normalize( pickRay.translation );
			m_grabTarget = b3Add( origin, b3MulSV( m_grabDepth, dir ) );
		}
		else
		{
			Sample::MouseMove( p );
		}
	}

	static Sample* Create( SampleContext* context )
	{
		return new RigidAtomsSample( context );
	}

	RAffWorld* m_raff = nullptr;
	b3MolFile m_mol;
	bool m_molLoaded;
	int m_pickedAtom;
	bool m_grabActive;
	float m_grabDepth;
	b3Vec3 m_grabTarget;

	// Tunable parameters
	float m_gravity;
	float m_restitution;
	float m_friction;
	float m_grabStiffness;
	int m_reportInterval;

	// Diagnostics
	float m_keHistory[KE_HISTORY_SIZE];
	int m_keHistoryCount;
	int m_keHistoryWrite;
	bool m_nanDetected;
};

static int sampleRigidAtoms = RegisterSample( "Molecules", "RigidAtomFF", RigidAtomsSample::Create );
