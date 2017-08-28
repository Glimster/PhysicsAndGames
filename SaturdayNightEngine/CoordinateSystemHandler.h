#pragma once
class CoordinateSystemHandler
{
public:
  CoordinateSystemHandler();
  CoordinateSystemHandler( float world2DisplayLength, const Eigen::Vector2f& originInDisplayCS );
  ~CoordinateSystemHandler();

  inline float world2DisplayLength() const { return world2DisplayLength_; }

  // Display CS:            World CS:
  // ------------->  x            ^ y
  // |                            |
  // |                            |
  // |                     -------|------> x
  // |                            |
  // v  y                         |
  sf::Vector2f convertToDisplayCS( const Eigen::Vector2f& position ) const;
  Eigen::Vector2f convertToWorldCS( const sf::Vector2f& position ) const;

  // Also converts from radians to degrees
  float computeAngleInDisplayCS( const Eigen::Vector2f& orientationInWorldCS ) const;

private:
  float world2DisplayLength_;
  Eigen::Vector2f world2DisplayLengthV_;
  Eigen::Vector2f originInDisplayCS_;
};

