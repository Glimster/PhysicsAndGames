#pragma once
#include "Category.h"

struct Command;
class CommandQueue;

class SceneNode : public sf::Drawable, private sf::NonCopyable
{
public:
  typedef std::unique_ptr< SceneNode > Ptr;
  typedef std::pair< SceneNode*, SceneNode* > Pair;

public:
  SceneNode( Category::Type category = Category::None );

  void attachChild( Ptr child );
  Ptr detachChild( const SceneNode& node );

  void update( sf::Time dt, CommandQueue& commands );

  // Excecutes command if the right category and the forwards to children
  void onCommand( const Command& command, sf::Time dt );
  virtual unsigned int getCategory() const;

  // Kollisionskontrol från SFML-boken. Långsamt!
  virtual sf::FloatRect	getBoundingRectDisplayCS() const;
  sf::Transform	getWorldTransformDisplayCS() const;
  const sf::Transform& getTransformDisplayCS() const { return transformable_.getTransform(); }
  //bool collision( const SceneNode& left, const SceneNode& right ) const; // Använder CollisionControl istället. Kan vara static?
  void checkSceneCollision( SceneNode& sceneGraph, std::set< Pair >& collisionPairs );
  void checkNodeCollision( SceneNode& node, std::set< Pair >& collisionPairs );
  void drawBoundingRect( sf::RenderTarget& target, sf::RenderStates states ) const;
  virtual bool isDestroyed() const;


  void print( std::string& string ) const;

protected:
  void printDefaultData_( std::string& string ) const;

private:
  virtual void draw( sf::RenderTarget& target, sf::RenderStates states ) const final;
  
  virtual void drawCurrent_( sf::RenderTarget& target, sf::RenderStates states ) const;  
  virtual void updateCurrent_( sf::Time dt, CommandQueue& commands );

  virtual void printCurrent_( std::string& string ) const;
protected:
  sf::Transformable transformable_;

private:
  std::vector< Ptr > children_;
  SceneNode* parent_;

  Category::Type defaultCategory_;
};