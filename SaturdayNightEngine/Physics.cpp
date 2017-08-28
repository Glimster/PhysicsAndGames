#include "stdafx.h"
#include "Physics.h"

using namespace std;
using namespace Phys;
float Physics::computeRotationalInertiaForRectangle( float mass, float a, float b )
{
  return mass * ( a * a + b * b ) / 12.0f;
}

Eigen::Vector2f Physics::circularOrbitVelocity( const Eigen::Vector2f& relPos, float centralMass )
{
  const float r = sqrt( relPos.dot( relPos ) );
  float speed = sqrt( PhysicalConstants::kSquared * centralMass / r );
  Eigen::Vector2f velocity( relPos(1), -relPos(0) );
  velocity.normalize();
  return velocity * speed;
}

float Physics::escapeVelocity( const Eigen::Vector2f& relPos, float centralMass )
{
  const float r = sqrt( relPos.dot( relPos ) );
  return sqrt( 2.0f * PhysicalConstants::kSquared * centralMass / r );
}

void Physics::elasticCollisionCircularObject( const Eigen::Vector2f& p1,
                                              const Eigen::Vector2f& p2, 
                                              const Eigen::Vector2f& v1,
                                              const Eigen::Vector2f& v2,
                                              const float m1,
                                              const float m2,
                                              Eigen::Vector2f& v1p,
                                              Eigen::Vector2f& v2p )
{
  Eigen::Vector2f normal = p1 - p2;
  normal.normalize();

  v1p = v1 - 2.0f * m2 / (m1 + m2) * (v1 - v2).dot( normal ) * normal;
  v2p = v2 - 2.0f * m1 / (m1 + m2) * (v2 - v1).dot( normal ) * normal;
}