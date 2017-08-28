#pragma once

#include "ResourceIdentifiers.h"
#include "GameEntity.h"
#include "Command.h"
#include "Projectile.h"

class TextNode;
class SpriteNode;
class CommandQueue;

class StarShip : public GameEntity
{
public:
  enum Type
  {
    Standard,
    TypeCount
  };

public:
  StarShip( Type type,
            const TextureHolder& textureHolder,
            const FontHolder& fonts,
            const CoordinateSystemHandler& csHandler,
            std::vector< PhysicalObject* >& physicalObjects );
  
  ~StarShip();

  virtual float getMass() const override; // Includes StarShip addons
  virtual void computeLinearForceOverM( Eigen::Vector2f& fOverM ) const override;
  virtual float computeTorque() const override;

  virtual sf::FloatRect getBoundingRectDisplayCS() const override;

  void setAftThrusters( bool thrusters ) { aftThrusters_ = thrusters; }
  void setLeftRotationThrusters( bool thrusters ) { leftRotationThrusters_ = thrusters; }
  void setRightRotationThrusters( bool thrusters ) { rightRotationThrusters_ = thrusters; }

  void fire();

  virtual unsigned int getCategory() const override;

  friend class Game; // TODO, bara för debug antar jag
private:
  virtual void drawCurrent_( sf::RenderTarget& target, sf::RenderStates states ) const override;
  virtual void updateCurrent_( sf::Time dt, CommandQueue& commands ) override;
  virtual void printCurrent_( std::string& string ) const override;

  void checkProjectileLaunch_( sf::Time dt, CommandQueue& commands );
  
  void launchCannonBall_( SceneNode& node, 
                          Projectile::Type type, 
                          std::vector< PhysicalObject* >& physicalObjects,
                          const TextureHolder& textureHolder );

private:
  Type type_;
  sf::Sprite sprite_;
  TextNode* healthDisplay_;

  bool aftThrusters_;
  bool leftRotationThrusters_;
  bool rightRotationThrusters_;

  Command fireCommand_;
  bool isFiring_;
  sf::Time fireCountdown_;
  float fireRate_;
  
  float cannonBallLaunchSpeed_;
  float cannonBallMass_;
  SpriteNode* cannonBall_;
};