#include "stdafx.h"
#include "Planet.h"

#include "CoordinateSystemHandler.h"
#include "ResourceHolder.h"

#include "Utilities.h"

using namespace std;

Planet::Planet( const std::string& name, float mass, float radius,
                int hitPoints,
                const TextureHolder& textureHolder,
                const CoordinateSystemHandler& csHandler ):
  GameEntity( mass, hitPoints, csHandler ),
  name_( name ),
  radius_( radius ),
  sprite_( textureHolder.get( toTextureID_()) )
{
  // TODO, consider making utility function as well
  const float scaleX = radius * csHandler_.world2DisplayLength() * 2.0f / sprite_.getTextureRect().width;
  const float scaleY = radius * csHandler_.world2DisplayLength() * 2.0f / sprite_.getTextureRect().height;
  transformable_.setScale( scaleX, scaleY );

  Utilities::centerOrigin( sprite_ );
}

Planet::~Planet()
{}

sf::FloatRect Planet::getBoundingRectDisplayCS() const
{
  return getWorldTransformDisplayCS().transformRect( sprite_.getGlobalBounds() );
}

unsigned int Planet::getCategory() const
{
  return Category::Planet;
}

void Planet::drawCurrent_(sf::RenderTarget& target, sf::RenderStates states) const
{
  target.draw(sprite_, states);
}

void Planet::updateCurrent_( sf::Time dt, CommandQueue& commands )
{
  transformable_.setPosition( csHandler_.convertToDisplayCS( position_ ) );
}

void Planet::printCurrent_( string& string ) const
{
  string.append( "Planet: category = " + to_string( getCategory() ) );
  printDefaultData_( string );
  string.append( ", " + name_ );
  printPOData_( string );
  string.append( "\n" );
}

Textures::ID Planet::toTextureID_() const
{
  if( name_ == "Mercury" )
  {
    return Textures::ID::Mercury;
  }
  else if( name_ == "Venus" )
  {
    return Textures::ID::Venus;
  }
  else if( name_ == "Earth" )
  {
    return Textures::ID::Earth;
  }
  else if( name_ == "Mars" )
  {
    return Textures::ID::Mars;
  }
  else if( name_ == "Jupiter" )
  {
    return Textures::ID::Jupiter;
  }
  else if( name_ == "Saturn" )
  {
    return Textures::ID::Saturn;
  }
  else if( name_ == "Uranus" )
  {
    return Textures::ID::Uranus;
  }
  else if( name_ == "Neptune" )
  {
    return Textures::ID::Neptune;
  }
  else if( name_ == "Sun" )
  {
    return Textures::ID::Sun;
  }
  else
  {
    throw invalid_argument( "Missing texture " + name_ );
    return Textures::ID::Sun;
  }
}