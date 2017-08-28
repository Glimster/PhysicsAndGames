#pragma once
#include "GameEntity.h"
#include "ResourceIdentifiers.h"

class CommandQueue;
class Projectile : public GameEntity
{
public:
  enum Type
  {
    CannonBall,
    TypeCount
  };

  Projectile( Type type, const TextureHolder& textureHolder, const CoordinateSystemHandler& csHandler );
  ~Projectile();

  float getRadius() const; // TODO, det här funkar bara nu, om det är en CannonBall. Styr upp!

  //void guideTowards( const sf::Vector2f& position ); 
  //bool isGuided() const; 
  virtual unsigned int getCategory() const override;
  virtual sf::FloatRect getBoundingRectDisplayCS() const; 
  //float getMaxSpeed() const; 
  int getDamage() const; 

private: 
  virtual void drawCurrent_( sf::RenderTarget& target, sf::RenderStates states ) const override; 
  virtual void updateCurrent_( sf::Time dt, CommandQueue& commands ) override;
  virtual void printCurrent_( std::string& string ) const override;

private: 
  Type type_;
  sf::Sprite sprite_; 
  //sf::Vector2f targetDirection_;
};

