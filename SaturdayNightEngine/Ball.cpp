#include "stdafx.h"
#include "Ball.h"

#include "Utilities.h"

Ball::Ball( float radius, float mass, const CoordinateSystemHandler& csHandler ):
  PhysicalObject( mass ),
  SceneNode(),
  csHandler_( csHandler ),
  shape_( radius ),
  radius_( radius )
{
  Utilities::centerOrigin( shape_ );
}


Ball::~Ball()
{}

sf::FloatRect Ball::getBoundingRectDisplayCS() const
{
  return getWorldTransformDisplayCS().transformRect( shape_.getGlobalBounds() );
}

unsigned int Ball::getCategory() const
{
  return Category::Ball;
}

void Ball::drawCurrent_( sf::RenderTarget& target, sf::RenderStates states ) const
{
  target.draw( shape_, states );
}
void Ball::updateCurrent_( sf::Time dt, CommandQueue& commands )
{
  transformable_.setPosition( csHandler_.convertToDisplayCS( position_ ) );
}

void Ball::printCurrent_( std::string& string ) const
{
  string.append( "Ball: " );
  printDefaultData_( string );
  printPOData_( string );
  string.append( "\n" );
}
