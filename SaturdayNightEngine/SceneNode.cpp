#include "stdafx.h"
#include "SceneNode.h"

#include "Command.h"
#include "Category.h"

#include "CollisionControl.h"

using namespace std;

SceneNode::SceneNode( Category::Type category ):
  children_(),
  parent_( nullptr ),
  defaultCategory_( category )
{}

void SceneNode::attachChild( Ptr child )
{
  child->parent_ = this;
  children_.push_back( std::move( child ) );
}

SceneNode::Ptr SceneNode::detachChild( const SceneNode& node )
{
  auto found = find_if( children_.begin(), children_.end(), 
    [&]( Ptr& p ) { return p.get() == &node; } );
  
  assert( found != children_.end() );

  Ptr result = std::move( *found );
  result->parent_ = nullptr;
  children_.erase( found );

  return result;
}

void SceneNode::update( sf::Time dt, CommandQueue& commands )
{
  updateCurrent_( dt, commands );

  for( const Ptr& child : children_ )
  {
    child->update( dt, commands );
  }
}

void SceneNode::onCommand( const Command& command, sf::Time dt )
{
  // NB, bitwise and intentional!
  if( getCategory() & command.category )
  {
    command.action( *this, dt );
  }

  for( Ptr& child : children_ )
  {
    child->onCommand( command, dt );
  }
}

unsigned int SceneNode::getCategory() const
{
  return defaultCategory_;
}

sf::FloatRect SceneNode::getBoundingRectDisplayCS() const
{
  return sf::FloatRect();
}

sf::Transform SceneNode::getWorldTransformDisplayCS() const
{
  sf::Transform transform = sf::Transform::Identity;

  for( const SceneNode* node = this; node != nullptr; node = node->parent_ )
    transform = node->getTransformDisplayCS() * transform;

  return transform;
}

//bool SceneNode::collision( const SceneNode& left, const SceneNode& right ) const
//{
//  return left.getBoundingRectDisplayCS().intersects( right.getBoundingRectDisplayCS() );
//}

void SceneNode::checkSceneCollision( SceneNode& sceneGraph, set< SceneNode::Pair >& collisionPairs )
{
  checkNodeCollision( sceneGraph, collisionPairs );

  for( Ptr& child : sceneGraph.children_ )
  {
    checkSceneCollision( *child, collisionPairs );
  }
}

void SceneNode::checkNodeCollision( SceneNode& node, set< SceneNode::Pair >& collisionPairs )
{
  if( this != &node &&
      CollisionControl::isColliding( *this, node ) &&
      //collision( *this, node ) &&
      !isDestroyed() &&
      !node.isDestroyed() )
  {
    // Inserts SceneNode pairs (always sorted the same way)
    collisionPairs.insert( minmax( this, &node ) );
  }
 
  for( Ptr& child : children_ )
  {
    child->checkNodeCollision( node, collisionPairs );
  }
}

void SceneNode::drawBoundingRect( sf::RenderTarget& target, sf::RenderStates ) const
{
  sf::FloatRect rect = getBoundingRectDisplayCS();

  sf::RectangleShape shape;
  shape.setPosition( sf::Vector2f( rect.left, rect.top ) );
  shape.setSize( sf::Vector2f( rect.width, rect.height ) );
  shape.setFillColor( sf::Color::Transparent );
  shape.setOutlineColor( sf::Color::Green );
  shape.setOutlineThickness( 1.f );

  target.draw( shape );
}

bool SceneNode::isDestroyed() const
{
  return false;
}

void SceneNode::print( string& string ) const
{
  this->printCurrent_( string );

  for( const Ptr& child : children_ )
  {
    child->print( string );
  }
}

void SceneNode::draw( sf::RenderTarget& target, sf::RenderStates states ) const
{
  states.transform *= transformable_.getTransform();

  drawCurrent_( target, states );

  for( const Ptr& child : children_ )
  {
    child->draw( target, states );
  }

  //drawBoundingRect( target, states );

}

void SceneNode::drawCurrent_( sf::RenderTarget& target, sf::RenderStates states ) const
{}

void SceneNode::updateCurrent_( sf::Time dt, CommandQueue& commands )
{}

void SceneNode::printCurrent_( string& string ) const
{
  string.append( "SceneNode: category = " + to_string( getCategory() ) );
  printDefaultData_( string );
  string.append( "\n" );
}

void SceneNode::printDefaultData_( string& string ) const
{
  string.append( ", #c = " + to_string( children_.size() ) );
  if( parent_ == nullptr )
    string.append( ". Top node!" );
}