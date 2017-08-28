#include "stdafx.h"
#include "PhysicalConstants.h"

using namespace Phys;

// D - mean solar day
// S - Solar mass

// Astronomical system of units
float PhysicalConstants::k = 0.01720209895f; // au^3/2 D^-1 S^-1/2
float PhysicalConstants::kSquared = 0.01720209895f * 0.01720209895f; // au^3 D^-2 S^-1

float PhysicalConstants::au2m = 1.49597870700e11f; // m au^-1
float PhysicalConstants::D2s = 86400.0f; // s D^-1
float PhysicalConstants::S2kg = 1.989e30f; // kg S^-1

// SI
float PhysicalConstants::G = 6.67384e-11f; // m^3 kg^-1 s^-2 (el. N m^2 kg^-2)

// ---
float PhysicalConstants::timeScale = 1.0f;