#pragma once
#include "PhysicalObject.h"
#include "SceneNode.h"
#include "CoordinateSystemHandler.h"

class Ball : public PhysicalObject, public SceneNode
{
public:
  Ball( float radius, float mass, const CoordinateSystemHandler& csHandler_ );
  ~Ball();

  float getRadius() const { return radius_; }

  virtual sf::FloatRect getBoundingRectDisplayCS() const override;

  virtual unsigned int getCategory() const override;
private:
  virtual void drawCurrent_( sf::RenderTarget& target, sf::RenderStates states ) const override;
  virtual void updateCurrent_( sf::Time dt, CommandQueue& commands ) override;
  virtual void printCurrent_( std::string& string ) const override;

private:
  const CoordinateSystemHandler& csHandler_;
  sf::CircleShape shape_;
  float radius_;
};

