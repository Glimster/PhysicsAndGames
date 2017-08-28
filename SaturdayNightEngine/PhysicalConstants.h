#pragma once

namespace Phys {

class PhysicalConstants
{
public:

  // Astronomical system of units
  static float k; // Gaussian gravitational constant
  static float kSquared; // Gaussian gravitational constant squared

  static float au2m; // Astronomical unit to meter
  static float D2s; // Mean solar day to seconds
  static float S2kg; // Solar mass to kg

  // SI
  static float G; // Gravitational constant
  
  // ---
  static float timeScale; // Scales the pace of time. Game.maxTimeStep_ shall ensure that physics simulation still holds.

};

}