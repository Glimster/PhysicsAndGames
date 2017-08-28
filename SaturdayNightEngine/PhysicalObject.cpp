#include "stdafx.h"
#include "PhysicalObject.h"

#include "MathUtil.h"
#include "Utilities.h"

using namespace std;
PhysicalObject::PhysicalObject( float mass ):
  Massive( mass ),
  position_( 0.0f, 0.0f ),
  velocity_( 0.0f, 0.0f ),
  orientation_( 0.0f, 0.0f )
{}

PhysicalObject::~PhysicalObject()
{}

void PhysicalObject::setPosition( const Eigen::Vector2f& position )
{
  position_ = position;
}

void PhysicalObject::setVelocity( const Eigen::Vector2f& velocity )
{
  velocity_ = velocity;
}

void PhysicalObject::setOrientation( const Eigen::Vector2f& orientation )
{
  orientation_ = orientation;
  orientation_.normalize();
}

void PhysicalObject::updateState( const Eigen::Vector2f& dr, const Eigen::Vector2f& dv )
{
  position_ += dr;
  velocity_ += dv;
}

void PhysicalObject::rotate( const float angle )
{
  MathUtil::rotate2D( orientation_, angle );
}

void PhysicalObject::printPOData_( string& string ) const
{
  string.append( ", p = " + Utilities::toString( position_ ) );
  string.append( ", v = " + Utilities::toString( velocity_ ) );
  string.append( ", o = " + Utilities::toString( orientation_ ) );
  string.append( ", mass = " + to_string( getMass() ) );
}