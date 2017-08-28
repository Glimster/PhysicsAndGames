#pragma once
#include "Massive.h"

class PhysicalObject : public Massive
{
public:
  PhysicalObject( float mass );
  ~PhysicalObject();

  // TODO, är det dumt att ge tillbaks const& här?
  const Eigen::Vector2f& getPosition() const { return position_; }
  const Eigen::Vector2f& getVelocity() const { return velocity_; }
  const Eigen::Vector2f& getOrientation() const { return orientation_; }
  void setPosition( const Eigen::Vector2f& position );
  void setVelocity( const Eigen::Vector2f& velocity );
  void setOrientation( const Eigen::Vector2f& orientation );

  void updateState( const Eigen::Vector2f& dr, const Eigen::Vector2f& dv );
  void rotate( const float angle );

protected:
  void printPOData_( std::string& string ) const;

protected:
  Eigen::Vector2f position_;
  Eigen::Vector2f velocity_;
  Eigen::Vector2f orientation_;

private:
};

