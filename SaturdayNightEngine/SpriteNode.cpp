#include "stdafx.h"
#include "SpriteNode.h"

#include "Utilities.h"

using namespace std;
SpriteNode::SpriteNode( const sf::Texture& texture, const sf::Vector2f& scaleFactor ):
  sprite_( texture )
{
  Utilities::centerOrigin( sprite_ );
  transformable_.setScale( scaleFactor );
}

SpriteNode::SpriteNode( const sf::Texture& texture, 
                        const sf::IntRect& textureRect, 
                        const sf::Vector2f& scaleFactor ):
  sprite_( texture, textureRect )
{
  Utilities::centerOrigin( sprite_ );
  transformable_.setScale( scaleFactor );
}

void SpriteNode::drawCurrent_( sf::RenderTarget& target, sf::RenderStates states ) const
{
  target.draw( sprite_, states );
}

void SpriteNode::printCurrent_( string& string ) const
{
  string.append( "SpriteNode: category = " + to_string( getCategory() ) + ", " + Utilities::toString( transformable_.getPosition() ) );
  printDefaultData_( string );
  string.append( "\n" );
}