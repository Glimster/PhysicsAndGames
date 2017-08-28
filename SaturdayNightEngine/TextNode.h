#pragma once
#include "SceneNode.h"

#include "ResourceHolder.h"
#include "ResourceIdentifiers.h"

class TextNode : public SceneNode
{
public:
  explicit TextNode( const FontHolder& fonts, const std::string& text );
  ~TextNode();

  inline void setPositionDisplayCS( float x, float y ) { transformable_.setPosition( x, y ); }
  inline void setRotationDisplayCS( float rotation ) { transformable_.setRotation( rotation ); }
  void setString( const std::string& text );

private:
  void drawCurrent_( sf::RenderTarget& target, sf::RenderStates states ) const override;
  virtual void printCurrent_( std::string& string ) const override;

private:
  sf::Text text_;
};

