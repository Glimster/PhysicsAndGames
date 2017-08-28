#include "stdafx.h"
#include "CoordinateSystemHandler.h"

using namespace std;

// TODO, inte så snyggt, kanske ska ha en initialize-metod?
CoordinateSystemHandler::CoordinateSystemHandler():
  world2DisplayLength_( 1.0f ),
  world2DisplayLengthV_( Eigen::Vector2f( world2DisplayLength_, -world2DisplayLength_ ) ),
  originInDisplayCS_( Eigen::Vector2f() )
{}

CoordinateSystemHandler::CoordinateSystemHandler( float world2DisplayLength, const Eigen::Vector2f& originInDisplayCS ):
  world2DisplayLength_( world2DisplayLength ),
  world2DisplayLengthV_( Eigen::Vector2f( world2DisplayLength_, -world2DisplayLength_ ) ),
  originInDisplayCS_( originInDisplayCS )
{}

CoordinateSystemHandler::~CoordinateSystemHandler()
{}

sf::Vector2f CoordinateSystemHandler::convertToDisplayCS( const Eigen::Vector2f& position ) const
{
  Eigen::Vector2f posE = position.cwiseProduct( world2DisplayLengthV_ ) + originInDisplayCS_;
  return sf::Vector2f( posE(0), posE(1) );
}

Eigen::Vector2f CoordinateSystemHandler::convertToWorldCS( const sf::Vector2f& position ) const
{
  Eigen::Vector2f posE( position.x, position.y );
  return (posE - originInDisplayCS_).cwiseQuotient( world2DisplayLengthV_ );
}

float CoordinateSystemHandler::computeAngleInDisplayCS( const Eigen::Vector2f& orientationInWorldCS ) const
{
  return atan2( -orientationInWorldCS(1), orientationInWorldCS(0) ) * 180.0f / float( M_PI );
}