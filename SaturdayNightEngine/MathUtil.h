#pragma once
class PhysicalObject;

// TODO, ska detta också vara ett namespace, pss som Utilities? Verkar vara bättre c++-standard!
class MathUtil
{
public:

  // Ellipses
  // a - semi major axis
  // b - semi minor axis
  
  static float fociCenterDistanceEllipse( float a, float b );
  static float areaEllipse( float a, float b );

  // Ramanujan approximation, 
  // se https://www.mathsisfun.com/geometry/ellipse-perimeter.html
  static float circumferenceEllipse( float a, float b );
  

  // C++11 baby!
  template< typename VecType >
  inline static auto angleBetween( const VecType& v1, VecType& v2 ) -> decltype( v1.norm() )
  {
    return acos( v1.dot( v2 ) / ( v1.norm() * v2.norm() ) );
  }

  // Eigen extensions
  template< typename VecType >
  inline static auto cross2D( const VecType& v1, VecType& v2 ) -> decltype( v1.norm() )
  {
    return abs( v1(0) * v2(1) - v1(1) * v2(0) );
  }

  // Counter clockwise rotation of vector
  template< typename VecType >
  inline static void rotate2D( VecType& vector, const float angle )
  {
    // Eigen for some reason almost 10 times slower!?
#if 0
    const Eigen::Rotation2D< float > rot( angle );
    vector = rot * vector;
#else
    const float cosAngle = cos( angle );
    const float sinAngle = sin( angle );
    const float newX = vector(0) * cosAngle - vector(1) * sinAngle;
    const float newY = vector(0) * sinAngle + vector(1) * cosAngle;
    vector(0) = newX;
    vector(1) = newY;
#endif
  }

  template< typename FuncType >
  inline static auto forwardEulerStep( FuncType& derivative,
                                       const float dt ) -> decltype( derivative() )
  {
    return derivative() * dt;
  };

  template< typename IncrementType >
  inline static auto RK4Step( const IncrementType k1,
                              const IncrementType k2,
                              const IncrementType k3,
                              const IncrementType k4,
                              const float dt ) -> decltype( k1 )
  {
    return oneOverSix_ * (k1 + 2.0f * (k2 + k3) + k4) * dt;
  }

private:

  static const float oneOverSix_;
};
