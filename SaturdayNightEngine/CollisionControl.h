#pragma once
#include "SceneNode.h"

class PhysicalObject;
class Ball;
namespace CollisionControl
{
  // For physics tests
  void handleCollisions( const std::vector< Ball* >& balls );
  void handleCollisionsBalls( SceneNode& sceneGraph );

  // Using the SFML book way
  void handleCollisions( SceneNode& sceneGraph );
  bool matchesCategories( SceneNode::Pair& colliders, Category::Type type1, Category::Type type2 );
  // Istället för att kolla direkt i SceneNode. Hanterar cirkulära objekt.
  bool isColliding( const SceneNode& left, const SceneNode& right );
};

