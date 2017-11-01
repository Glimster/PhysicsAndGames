#pragma once
#include "IStepper.h"

template< typename T >
struct RK4 : IStepper
{
public:
  RK4( double& x,
       std::vector< double >& y,
       std::vector< double >& dydx );

  void dy( const double dx, T& derivatives );
  void step( const double dx, T& derivatives );

  typedef T DType;

private:
  const double oneOverSix_ = 1.0 / 6.0;
};

template< typename T >
RK4< T >::RK4( double& x,
               std::vector< double > &y,
               std::vector< double > &dydx ) :
  IStepper( x, y, dydx )
{}

template< typename T >
void RK4< T >::dy( const double dx, T& derivatives )
{
  const double dxHalf = dx / 2.0;
  const size_t nbOfVariables = y_.size();
  // TODO, lägg in k som privata variabler?
  vector< double > k1( nbOfVariables ), k2( nbOfVariables ), k3( nbOfVariables ), k4( nbOfVariables );
  vector< double > yTemp( nbOfVariables );

  derivatives( x_, y_, k1 );

  std::transform( y_.begin(), y_.end(), k1.begin(), yTemp.begin(), [dxHalf]( const double y, const double dydx ) { return y + dxHalf * dydx; } );
  derivatives( x_ + dxHalf, yTemp, k2 );

  std::transform( y_.begin(), y_.end(), k2.begin(), yTemp.begin(), [dxHalf]( const double y, const double dydx ) { return y + dxHalf * dydx; } );
  derivatives( x_ + dxHalf, yTemp, k3 );

  std::transform( y_.begin(), y_.end(), k3.begin(), yTemp.begin(), [dx]( const double y, const double dydx ) { return y + dx * dydx; } );
  derivatives( x_ + dx, yTemp, k4 );

  for( size_t i = 0; i < nbOfVariables; ++i )
    y_[i] += oneOverSix_ * (k1[i] + 2.0 * (k2[i] + k3[i]) + k4[i]) * dx;
}

template< typename T>
void RK4< T >::step( const double dx, T& derivatives )
{
  dy( dx, derivatives );
  x_ += dx;
}