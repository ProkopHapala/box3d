// SPDX-FileCopyrightText: 2025 Erin Catto
// SPDX-License-Identifier: MIT

// Molecular dynamics demo using rigid ports (spherical joints).
//
// Each atom is a dynamic body with a sphere shape. Each bond is modeled as a
// spherical joint where the anchor on body A is offset from A's center by the
// initial bond vector (the "port"). This simultaneously constrains the bond
// length and generates torque on A to align the port toward B — the "as rigid
// as possible" behavior from RRsp3. When an atom has multiple bonds, the SI
// solver iterates over all its port constraints to find the best compromise
// orientation.

#include "gfx/draw.h"
#include "imgui.h"
#include "mol_loader.h"
#include "sample.h"

#include "box3d/box3d.h"
#include "box3d/math_functions.h"

#include <math.h>

static bool TryLoadMol( b3MolFile* mol )
{
	// Try common paths relative to the working directory and source tree.
	const char* paths[] = {
		"data/mol/nNonan.mol",
		"../data/mol/nNonan.mol",
		"../../data/mol/nNonan.mol",
		"../../../data/mol/nNonan.mol",
	};

	for ( int i = 0; i < (int)( sizeof( paths ) / sizeof( paths[0] ) ); ++i )
	{
		if ( b3LoadMolFile( paths[i], mol ) )
		{
			return true;
		}
	}

	return false;
}

class MoleculesRigidPorts : public Sample
{
public:
	explicit MoleculesRigidPorts( SampleContext* context )
		: Sample( context )
	{
		if ( context->restart == false )
		{
			m_camera->SetView( -35.0f, 25.0f, 40.0f, { 0.0f, 5.0f, 0.0f } );
		}

		AddGroundBox( 40.0f );

		m_mol = {};
		if ( TryLoadMol( &m_mol ) == false )
		{
			DrawTextLine( "ERROR: could not load nNonan.mol" );
			return;
		}

		m_moleculeCount = 10;
		m_jointHertz = 0.0f;
		m_jointDampingRatio = 0.5f;
		m_useSpring = false;

		SpawnMolecules();
	}

	void SpawnMolecules()
	{
		// Place molecules in a grid above the ground
		int cols = 5;
		float spacing = 15.0f;
		float yOffset = 12.0f;

		for ( int m = 0; m < m_moleculeCount; ++m )
		{
			float x = ( m % cols - cols / 2.0f ) * spacing;
			float z = ( m / cols ) * spacing;
			b3Vec3 offset = { x, yOffset, z };

			CreateMolecule( offset, m );
		}
	}

	void CreateMolecule( b3Vec3 offset, int moleculeIndex )
	{
		const b3MolFile& mol = m_mol;

		// Create a body for each atom
		b3BodyId bodyIds[B3_MOL_MAX_ATOMS];
		for ( int i = 0; i < mol.atomCount; ++i )
		{
			const b3MolAtom& atom = mol.atoms[i];

			b3BodyDef bodyDef = b3DefaultBodyDef();
			bodyDef.type = b3_dynamicBody;
			bodyDef.position = { offset.x + atom.x, offset.y + atom.y, offset.z + atom.z };
			bodyDef.linearDamping = 0.1f;
			bodyDef.angularDamping = 0.5f;
			bodyDef.enableSleep = false;
			bodyDef.isAwake = true;

			b3BodyId bodyId = b3CreateBody( m_worldId, &bodyDef );
			bodyIds[i] = bodyId;

			// Sphere shape with radius from covalent radius
			float radius = b3GetCovalentRadius( atom.element );
			float mass = b3GetAtomicMass( atom.element );
			float volume = ( 4.0f / 3.0f ) * B3_PI * radius * radius * radius;
			float density = mass / ( volume + 1e-12f );

			b3ShapeDef shapeDef = b3DefaultShapeDef();
			shapeDef.density = density;
			shapeDef.baseMaterial.friction = 0.2f;
			shapeDef.baseMaterial.restitution = 0.0f;

			b3Sphere sphere = { b3Vec3_zero, radius };
			b3CreateSphereShape( bodyId, &shapeDef, &sphere );
		}

		// Create spherical joint port constraints for each bond.
		// Two joints per bond: one with port on A toward B, one with port on B toward A.
		// This gives both atoms angular constraints from each bond, matching the RRsp3
		// approach where each atom has ports toward all its neighbors. With multiple
		// ports, the SI solver finds the orientation that best aligns all ports —
		// the "as rigid as possible" behavior.
		for ( int i = 0; i < mol.bondCount; ++i )
		{
			const b3MolBond& bond = mol.bonds[i];
			int a = bond.a;
			int b = bond.b;

			b3Vec3 ab = {
				mol.atoms[b].x - mol.atoms[a].x,
				mol.atoms[b].y - mol.atoms[a].y,
				mol.atoms[b].z - mol.atoms[a].z
			};
			b3Vec3 ba = { -ab.x, -ab.y, -ab.z };

			// Port on A toward B: A's port tip must coincide with B's center.
			// This generates torque on A via cross(rA, impulse) to align the port.
			{
				b3SphericalJointDef jointDef = b3DefaultSphericalJointDef();
				jointDef.base.bodyIdA = bodyIds[a];
				jointDef.base.bodyIdB = bodyIds[b];
				jointDef.base.localFrameA.p = ab;
				jointDef.base.localFrameA.q = b3Quat_identity;
				jointDef.base.localFrameB.p = b3Vec3_zero;
				jointDef.base.localFrameB.q = b3Quat_identity;
				jointDef.base.collideConnected = false;

				if ( m_useSpring )
				{
					jointDef.enableSpring = true;
					jointDef.hertz = m_jointHertz;
					jointDef.dampingRatio = m_jointDampingRatio;
				}

				b3CreateSphericalJoint( m_worldId, &jointDef );
			}

			// Port on B toward A: B's port tip must coincide with A's center.
			// This generates torque on B to align its port toward A.
			{
				b3SphericalJointDef jointDef = b3DefaultSphericalJointDef();
				jointDef.base.bodyIdA = bodyIds[b];
				jointDef.base.bodyIdB = bodyIds[a];
				jointDef.base.localFrameA.p = ba;
				jointDef.base.localFrameA.q = b3Quat_identity;
				jointDef.base.localFrameB.p = b3Vec3_zero;
				jointDef.base.localFrameB.q = b3Quat_identity;
				jointDef.base.collideConnected = false;

				if ( m_useSpring )
				{
					jointDef.enableSpring = true;
					jointDef.hertz = m_jointHertz;
					jointDef.dampingRatio = m_jointDampingRatio;
				}

				b3CreateSphericalJoint( m_worldId, &jointDef );
			}
		}
	}

	bool DrawControls() override
	{
		ImGui::PushItemWidth( 6.0f * ImGui::GetFontSize() );

		ImGui::Text( "Molecule: nNonan (C9H20)" );
		ImGui::Text( "Atoms: %d, Bonds: %d", m_mol.atomCount, m_mol.bondCount );
		ImGui::Text( "Molecules: %d", m_moleculeCount );

		ImGui::Separator();

		if ( ImGui::Checkbox( "Soft Spring", &m_useSpring ) )
		{
			// Will apply on next respawn
		}

		if ( m_useSpring )
		{
			ImGui::SliderFloat( "Hertz", &m_jointHertz, 1.0f, 30.0f, "%3.1f" );
			ImGui::SliderFloat( "Damping", &m_jointDampingRatio, 0.0f, 4.0f, "%3.1f" );
		}

		if ( ImGui::Button( "Respawn" ) )
		{
			// Destroy and recreate the world
			CreateWorld( &m_context->capacity );
			AddGroundBox( 40.0f );
			SpawnMolecules();
		}

		ImGui::PopItemWidth();
		return true;
	}

	void Step() override
	{
		Sample::Step();

		DrawTextLine( "Molecules: Rigid Ports (Spherical Joints)" );
		DrawTextLine( "Each bond = spherical joint with port offset on atom A" );
		DrawTextLine( "Atoms: %d  Bonds: %d  Molecules: %d",
			m_mol.atomCount, m_mol.bondCount, m_moleculeCount );
	}

	static Sample* Create( SampleContext* context )
	{
		return new MoleculesRigidPorts( context );
	}

	b3MolFile m_mol;
	int m_moleculeCount;
	float m_jointHertz;
	float m_jointDampingRatio;
	bool m_useSpring;
};

static int sampleMolecules = RegisterSample( "Molecules", "Rigid Ports", MoleculesRigidPorts::Create );
