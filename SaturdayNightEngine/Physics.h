#pragma once
#include "PhysicalConstants.h"

// TODO, ska detta också vara ett namespace, pss som Utilities? Verkar vara bättre c++-standard!
class Physics
{
public:

  // Force of gravity using Gaussian gravitational constant
  inline static Eigen::Vector2f forceOfGravityOverM( float m2, const Eigen::Vector2f& r )
  {
    Eigen::Vector2f F = r;
    F.normalize();
    // TODO, använd softening length!
    const float force = Phys::PhysicalConstants::kSquared * m2 / r.dot( r );
    F *= -force;
    return F;
  }

  inline static Eigen::Vector2f momentum( float m, const Eigen::Vector2f& v )
  {
    return m * v;
  }

  inline static float kineticEnergy( float m, const Eigen::Vector2f& v )
  {
    return m * v.dot( v ) / 2.0f;
  }

  // Note that the potential energy = 0 for r = oo
  inline static float gravitationalPotentialEnergy( float m1, float m2, const Eigen::Vector2f& r )
  {
    // TODO, använd softening length!
    return -Phys::PhysicalConstants::kSquared * m1 * m2 / r.norm();
  }

  // Rotation inertia around z at center of rectangle
  static float computeRotationalInertiaForRectangle( float mass, float a, float b );

  // Computes the velocity of a satellite in circular orbit around a large central mass
  // The velocity will be perpendicular to the relative position vector 
  static Eigen::Vector2f circularOrbitVelocity( const Eigen::Vector2f& relPos, float centralMass );
  
  // Computes the escape velocity (really the speed) of a satellite around a large central mass
  static float escapeVelocity( const Eigen::Vector2f& relPos, float centralMass );

  // Resulting velocities v1p and v2p of two circular objects colliding
  static void elasticCollisionCircularObject( const Eigen::Vector2f& p1,
                                              const Eigen::Vector2f& p2,
                                              const Eigen::Vector2f& v1,
                                              const Eigen::Vector2f& v2,
                                              const float m1,
                                              const float m2,
                                              Eigen::Vector2f& v1p,
                                              Eigen::Vector2f& v2p );

  
  // TODO - ny representation (och precision). Merge:a ihop med ovan!?

  // Nb of particles = stateVector.size() x 4
  // stateVector[i] = px
  // stateVector[i + 1] = py
  // stateVector[i + 2] = vx
  // stateVector[i + 3] = vy
  static double computeTotalEnergy2D( const std::vector< double >& stateVector, 
                                      const std::vector< double >& masses );

};

