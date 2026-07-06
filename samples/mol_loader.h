// SPDX-FileCopyrightText: 2025 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>

// Minimal MDL V2000 MOL file parser.
// Parses atom positions, element symbols, and bond connectivity.
// See https://en.wikipedia.org/wiki/Chemical_table_file

#define B3_MOL_MAX_ATOMS 512
#define B3_MOL_MAX_BONDS 512

typedef struct b3MolAtom
{
	float x, y, z;
	char element[4];
} b3MolAtom;

typedef struct b3MolBond
{
	int a;       // first atom index (0-based)
	int b;       // second atom index (0-based)
	int type;    // 1=single, 2=double, 3=triple
} b3MolBond;

typedef struct b3MolFile
{
	b3MolAtom atoms[B3_MOL_MAX_ATOMS];
	b3MolBond bonds[B3_MOL_MAX_BONDS];
	int atomCount;
	int bondCount;
	char name[256];
} b3MolFile;

static inline bool b3LoadMolFile( const char* path, b3MolFile* mol )
{
	FILE* fp = fopen( path, "r" );
	if ( fp == NULL )
	{
		return false;
	}

	memset( mol, 0, sizeof( b3MolFile ) );

	char line[1024];

	// Line 1: molecule name
	if ( fgets( line, sizeof( line ), fp ) == NULL ) { fclose( fp ); return false; }
	// Strip trailing whitespace
	int len = (int)strlen( line );
	while ( len > 0 && isspace( (unsigned char)line[len - 1] ) ) { line[--len] = '\0'; }
	strncpy( mol->name, line, sizeof( mol->name ) - 1 );

	// Line 2: program/timestamp
	if ( fgets( line, sizeof( line ), fp ) == NULL ) { fclose( fp ); return false; }

	// Line 3: comment
	if ( fgets( line, sizeof( line ), fp ) == NULL ) { fclose( fp ); return false; }

	// Line 4: counts line "aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv"
	if ( fgets( line, sizeof( line ), fp ) == NULL ) { fclose( fp ); return false; }

	int atomCount = 0, bondCount = 0;
	sscanf( line, "%3d%3d", &atomCount, &bondCount );

	if ( atomCount < 0 || atomCount > B3_MOL_MAX_ATOMS ||
		 bondCount < 0 || bondCount > B3_MOL_MAX_BONDS )
	{
		fclose( fp );
		return false;
	}

	mol->atomCount = atomCount;
	mol->bondCount = bondCount;

	// Atom block
	for ( int i = 0; i < atomCount; i++ )
	{
		if ( fgets( line, sizeof( line ), fp ) == NULL ) { fclose( fp ); return false; }

		float x, y, z;
		char element[8] = { 0 };

		// Format: xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaadd...
		// Coordinates are 10 chars each, then space, then 3-char element
		int parsed = sscanf( line, "%f %f %f %3s", &x, &y, &z, element );
		if ( parsed < 4 )
		{
			fclose( fp );
			return false;
		}

		mol->atoms[i].x = x;
		mol->atoms[i].y = y;
		mol->atoms[i].z = z;
		strncpy( mol->atoms[i].element, element, 3 );
		mol->atoms[i].element[3] = '\0';
	}

	// Bond block
	for ( int i = 0; i < bondCount; i++ )
	{
		if ( fgets( line, sizeof( line ), fp ) == NULL ) { fclose( fp ); return false; }

		int a, b, type;
		int parsed = sscanf( line, "%3d%3d%3d", &a, &b, &type );
		if ( parsed < 3 )
		{
			fclose( fp );
			return false;
		}

		// MOL files are 1-indexed; convert to 0-indexed
		mol->bonds[i].a = a - 1;
		mol->bonds[i].b = b - 1;
		mol->bonds[i].type = type;
	}

	fclose( fp );
	return true;
}

// Get atomic mass (in amu) for common elements
static inline float b3GetAtomicMass( const char* element )
{
	if ( element[0] == '\0' ) return 1.0f;

	// Check two-character symbols first
	char c0 = (char)toupper( (unsigned char)element[0] );
	char c1 = element[1];

	// Cl
	if ( c0 == 'C' && c1 == 'l' ) return 35.45f;
	// Br
	if ( c0 == 'B' && c1 == 'r' ) return 79.904f;

	switch ( c0 )
	{
		case 'H': return 1.008f;
		case 'C': return 12.011f;
		case 'N': return 14.007f;
		case 'O': return 15.999f;
		case 'F': return 18.998f;
		case 'P': return 30.974f;
		case 'S': return 32.06f;
		case 'B': return 10.81f;
		default:  return 12.0f;
	}
}

// Get van der Waals radius (in Angstroms) for common elements
static inline float b3GetVdwRadius( const char* element )
{
	if ( element[0] == '\0' ) return 1.2f;

	char c = (char)toupper( (unsigned char)element[0] );
	switch ( c )
	{
		case 'H': return 1.2f;
		case 'C': return 1.7f;
		case 'N': return 1.55f;
		case 'O': return 1.52f;
		case 'F': return 1.47f;
		case 'P': return 1.8f;
		case 'S': return 1.8f;
		default:  return 1.5f;
	}
}

// Get covalent radius (in Angstroms) for common elements
static inline float b3GetCovalentRadius( const char* element )
{
	if ( element[0] == '\0' ) return 0.31f;

	char c = (char)toupper( (unsigned char)element[0] );
	switch ( c )
	{
		case 'H': return 0.31f;
		case 'C': return 0.76f;
		case 'N': return 0.71f;
		case 'O': return 0.66f;
		case 'F': return 0.57f;
		case 'P': return 1.07f;
		case 'S': return 1.05f;
		default:  return 0.7f;
	}
}
