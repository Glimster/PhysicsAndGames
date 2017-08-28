#include "stdafx.h"
#include "TextNode.h"

#include "Utilities.h"

TextNode::TextNode( const FontHolder& fonts, const std::string& text )
{
  text_.setFont( fonts.get( Fonts::Main ) );
  text_.setCharacterSize( 200 ); // TODO, måste vara såhär stor eftersom skeppet är förminskat (syns inte annars)
  setString( text );
}

TextNode::~TextNode()
{}

void TextNode::setString( const std::string& text )
{
  text_.setString( text );
  Utilities::centerOrigin( text_ );
}

void TextNode::drawCurrent_( sf::RenderTarget& target, sf::RenderStates states ) const
{
  target.draw( text_, states );
}

void TextNode::printCurrent_( std::string& string ) const
{
  string.append( "TextNode: category = " + to_string( getCategory() ) );
  printDefaultData_( string );
  string.append( ", text = " + text_.getString() );
  string.append( "\n" );
}