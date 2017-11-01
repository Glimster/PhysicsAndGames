#pragma once
#include "IStepper.h"

// Forward (explicit) Euler method
template< typename T > 
struct EulerForward : IStepper 
{
  EulerForward( double& x,
                std::vector< double > &y,
                std::vector< double > &dydx );

  void dy( const double dx, T& derivatives );
  void step( const double dx, T& derivatives );

  typedef T DType;
};

template< typename T >
EulerForward< T >::EulerForward( double& x,
                                 std::vector< double >& y,
                                 std::vector< double >& dydx ) :
  IStepper( x, y, dydx )
{}

template< typename T >
void EulerForward< T >::dy( const double dx, T& derivatives )
{
  derivatives( x_, y_, dydx_ );
  std::transform( dydx_.begin(), dydx_.end(), y_.begin(), y_.begin(), [dx]( const double dydx, const double y ) { return y + dx * dydx; } );
}

template< typename T >
void EulerForward< T >::step( const double dx, T& derivatives )
{
  dy( dx, derivatives );
  x_ += dx;
}


// Symplectic (semi-implicit) Euler method
// TODO, kan man baka ihop 1D och 2D lätt?
template< typename T >
struct SymplecticEuler1D : IStepper 
{
  SymplecticEuler1D( double& x,
                     std::vector< double > &y,
                     std::vector< double > &dydx );

  void dy( const double dx, T& derivatives );
  void step( const double dx, T& derivatives );

  typedef T DType;
};

template< typename T >
SymplecticEuler1D< T >::SymplecticEuler1D( double& x,
                                           std::vector< double > &y,
                                           std::vector< double > &dydx ) :
  IStepper( x, y, dydx )
{
  if( y.size() % 2 != 0 || y.size() != dydx.size() )
    throw runtime_error( "Symplectic Euler assumes coupled ODE:s (e.g. position and velocity) in 1 dimension" );
}

template< typename T >
void SymplecticEuler1D< T >::dy( const double dx, T& derivatives )
{
  derivatives( x_, y_, dydx_ );
  for( size_t i = 0; i < dydx_.size(); i += 2 )
  {
    y_[i] += dydx_[i] * dx;
  }

  derivatives( x_, y_, dydx_ );
  for( size_t i = 0; i < dydx_.size(); i += 2 )
  {
    y_[i + 1] += dydx_[i + 1] * dx;
  }
}

template< typename T >
void SymplecticEuler1D< T >::step( const double dx, T& derivatives )
{
  dy( dx, derivatives );
  x_ += dx;
}

template< typename T >
struct SymplecticEuler2D : IStepper {

  SymplecticEuler2D( double& x,
                     std::vector< double > &y,
                     std::vector< double > &dydx );

  void dy( const double dx, T& derivatives );
  void step( const double dx, T& derivatives );

  typedef T DType;
};

template< typename T >
SymplecticEuler2D< T >::SymplecticEuler2D( double& x,
                                           std::vector< double > &y,
                                           std::vector< double > &dydx ) :
  IStepper( x, y, dydx )
{
  if( y.size() % 4 != 0 || y.size() != dydx.size() )
    throw runtime_error( "Symplectic Euler assumes coupled ODE:s (e.g. position and velocity) in 2 dimensions" );
}

template< typename T >
void SymplecticEuler2D< T >::dy( const double dx, T& derivatives )
{
  derivatives( x_, y_, dydx_ );
  for( size_t i = 0; i < dydx_.size(); i += 4 )
  {
    y_[i] += dydx_[i] * dx;
    y_[i + 1] += dydx_[i + 1] * dx;
  }
  
  derivatives( x_, y_, dydx_ );
  for( size_t i = 0; i < dydx_.size(); i += 4 )
  {
    y_[i + 2] += dydx_[i + 2] * dx;
    y_[i + 3] += dydx_[i + 3] * dx;
  }
}

template < typename T>
void SymplecticEuler2D< T >::step( const double dx, T& derivatives )
{
  dy( dx, derivatives );
  x_ += dx;
}