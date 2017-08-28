#include "stdafx.h"
#include "Utilities.h"

using namespace std;
void Utilities::centerOrigin( sf::Sprite& sprite )
{
  sf::FloatRect bounds = sprite.getLocalBounds();
  sprite.setOrigin( std::floor( bounds.left + bounds.width / 2.0f ), std::floor( bounds.top + bounds.height / 2.0f ) );
}

void Utilities::centerOrigin( sf::CircleShape& circle )
{
  sf::FloatRect bounds = circle.getLocalBounds();
  circle.setOrigin( std::floor( bounds.left + bounds.width / 2.0f ), std::floor( bounds.top + bounds.height / 2.0f ) );
}

void Utilities::centerOrigin( sf::Text& text )
{
  sf::FloatRect bounds = text.getLocalBounds();
  text.setOrigin( std::floor( bounds.left + bounds.width / 2.0f ), std::floor( bounds.top + bounds.height / 2.0f ) );
}

std::string Utilities::toString( const sf::Vector2f& vec )
{
  string str = "(" + to_string( vec.x ) + ", " + to_string( vec.y ) + ")";
  return str;
}

std::string Utilities::toString( const Eigen::Vector2f& vec )
{
  string str = "(" + to_string( vec(0) ) + ", " + to_string( vec(1) ) + ")";
  return str;
}