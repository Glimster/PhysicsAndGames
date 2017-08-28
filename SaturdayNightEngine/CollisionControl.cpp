#include "stdafx.h"
#include "CollisionControl.h"

#include "SceneNode.h"
#include "PhysicalObject.h"
#include "Planet.h"
#include "Projectile.h"
#include "Ball.h"
#include "StarShip.h"

#include "Category.h"
#include "Physics.h"

using namespace std;

void CollisionControl::handleCollisions( const std::vector< Ball* >& balls )
{
  set< SceneNode::Pair > collisionPairs;
  for( size_t i = 0; i < balls.size() - 1; ++i )
  {
    for( size_t j = i + 1; j < balls.size(); ++j )
    {
      float d = (balls[i]->getPosition() - balls[j]->getPosition()).norm();
      if( d <= balls[i]->getRadius() + balls[j]->getRadius() )
      {
        auto& b1 = *balls[i];
        auto& b2 = *balls[j];
        Eigen::Vector2f v1p, v2p;
        Physics::elasticCollisionCircularObject( b1.getPosition(),
                                                 b2.getPosition(),
                                                 b1.getVelocity(),
                                                 b2.getVelocity(),
                                                 b1.getMass(),
                                                 b2.getMass(),
                                                 v1p,
                                                 v2p );

        b1.setVelocity( v1p );
        b2.setVelocity( v2p );
      }
    }
  }
}

void CollisionControl::handleCollisionsBalls( SceneNode& sceneGraph )
{
  set< SceneNode::Pair > collisionPairs;
  sceneGraph.checkSceneCollision( sceneGraph, collisionPairs );

  for( SceneNode::Pair pair : collisionPairs )
  {
    if( matchesCategories( pair, Category::Ball, Category::Ball ) )
    {
      auto& b1 = static_cast<Ball&>(*pair.first);
      auto& b2 = static_cast<Ball&>(*pair.second);
      Eigen::Vector2f v1p, v2p;
      Physics::elasticCollisionCircularObject( b1.getPosition(),
                                               b2.getPosition(),
                                               b1.getVelocity(),
                                               b2.getVelocity(),
                                               b1.getMass(),
                                               b2.getMass(),
                                               v1p,
                                               v2p );

      b1.setVelocity( v1p );
      b2.setVelocity( v2p );
    }
  }
}

void CollisionControl::handleCollisions( SceneNode& sceneGraph )
{
  set< SceneNode::Pair > collisionPairs;
  sceneGraph.checkSceneCollision( sceneGraph, collisionPairs );

  for( SceneNode::Pair pair : collisionPairs )
  {
    // TODO: Hemsk, hemsk design. Använd inte SceneGraph för att kolla kollissioner!!!
    PhysicalObject* o1;
    PhysicalObject* o2;
    if( matchesCategories( pair, Category::Planet, Category::Planet ) )
    {
      o1 = &static_cast<Planet&>(*pair.first);
      o2 = &static_cast<Planet&>(*pair.second);
    }
    else if( matchesCategories( pair, Category::Projectile, Category::Projectile ) )
    {
      o1 = &static_cast<Projectile&>(*pair.first);
      o2 = &static_cast<Projectile&>(*pair.second);
    }
    else if( matchesCategories( pair, Category::Projectile, Category::Planet ) )
    {
      o1 = &static_cast<Projectile&>(*pair.first);
      o2 = &static_cast<Planet&>(*pair.second);
    }
    else if( matchesCategories( pair, Category::Player, Category::Planet ) )
    {
      o1 = &static_cast<StarShip&>(*pair.first);
      o2 = &static_cast<Planet&>(*pair.second);
    }
    else if( matchesCategories( pair, Category::Projectile, Category::Player ) )
    {
      o1 = &static_cast<Projectile&>(*pair.first);
      o2 = &static_cast<StarShip&>(*pair.second);
    }
    else
    {
      o1 = nullptr;
      o2 = nullptr;
    }

    if( o1 != nullptr && o2 != nullptr )
    {
      // TODO, not always circular object
      Eigen::Vector2f v1p, v2p;
      Physics::elasticCollisionCircularObject( o1->getPosition(),
                                               o2->getPosition(),
                                               o1->getVelocity(),
                                               o2->getVelocity(),
                                               o1->getMass(),
                                               o2->getMass(),
                                               v1p,
                                               v2p );

      o1->setVelocity( v1p );
      o2->setVelocity( v2p );
    }
  }
}

bool CollisionControl::matchesCategories( SceneNode::Pair& colliders, Category::Type type1, Category::Type type2 )
{
  auto category1 = colliders.first->getCategory();
  auto category2 = colliders.second->getCategory();

  if( type1 & category1 && type2 & category2 )
  {
    return true;
  }
  else if( type1 & category2 && type2 & category1 )
  {
    swap( colliders.first, colliders.second );
    return true;
  }
  else
  {
    return false;
  }
}

bool CollisionControl::isColliding( const SceneNode& left, const SceneNode& right )
{
  // TODO, jättefult, fixa detta!
  auto catL = left.getCategory();
  auto catR = right.getCategory();
  auto isBallL = catL & Category::Ball;
  auto isPlanetL = catL & Category::Planet;
  auto isProjectileL = catL & Category::Projectile;
  auto isBallR = catR & Category::Ball;
  auto isPlanetR = catR & Category::Planet;
  auto isProjectileR = catR & Category::Projectile;

  if( (isBallL || isPlanetL || isProjectileL) && (isBallR || isPlanetR || isProjectileR) )
  {
    Eigen::Vector2f p1, p2;
    float r1, r2;
    if( isBallL )
    {
      const auto& b1 = static_cast<const Ball&>(left);
      p1 = b1.getPosition();
      r1 = b1.getRadius();
    }
    else if( isPlanetL )
    {
      const auto& b1 = static_cast<const Planet&>(left);
      p1 = b1.getPosition();
      r1 = b1.getRadius();
    }
    else
    {
      const auto& b1 = static_cast<const Projectile&>(left);
      p1 = b1.getPosition();
      r1 = b1.getRadius();
    }

    if( isBallR )
    {
      const auto& b2 = static_cast<const Ball&>(right);
      p2 = b2.getPosition();
      r2 = b2.getRadius();
    }
    else if( isPlanetR )
    {
      const auto& b2 = static_cast<const Planet&>(right);
      p2 = b2.getPosition();
      r2 = b2.getRadius();
    }
    else
    {
      const auto& b2 = static_cast<const Projectile&>(right);
      p2 = b2.getPosition();
      r2 = b2.getRadius();
    }

    // Circle overlapping
    float d = (p1 - p2).norm();
    return d <= r1 + r2;
  }
  else
  {
    // BB overlapping
    return left.getBoundingRectDisplayCS().intersects( right.getBoundingRectDisplayCS() );
  }
}