#pragma once
#include <vector>

#include "NR_OdeOutput.h"


template< typename Stepper >
class OdeIntegrator
{
public:
  OdeIntegrator( double& x0,
                 std::vector< double >& y0,
                 const double dx,
                 const bool denseOutput,
                 typename Stepper::DType& derivatives );

  void integrate( const double deltaX );

  template< typename Type >
  void integrateUntil( Type& criterion );

  const std::vector< double >& getXDense() const { return xDense_; }
  const std::vector< std::vector< double > >& getYDense() const { return yDense_; }

private:

  //double minStep_;
  //size_t nVariables_;
  //std::vector< double > yStart_;
  //int nGood_, nBad_;
  //double t1_, t2_;
  double& x_;
  std::vector< double >& y_;
  std::vector< double > dydx_;
  const double dx_;

  // TODO: so far only stored at steps dx_, support interpolation
  bool denseOutput_;
  std::vector< double > xDense_;
  std::vector< std::vector< double > > yDense_;

  //NR_OdeOutput output_;

  typename Stepper::DType& derivatives_;
  Stepper stepper_;
};

// Implementation ----------------------
template< typename Stepper >
OdeIntegrator< Stepper >::OdeIntegrator( double& x0,
                                         std::vector< double >& y0,
                                         const double dx,
                                         const bool denseOutput,
                                         typename Stepper::DType& derivatives ) :
  x_( x0 ),
  y_( y0 ),
  dx_( dx ),
  dydx_( y_.size() ), // Behövs denna i denna klass?
  denseOutput_( denseOutput ),
  xDense_(),
  yDense_(),
  derivatives_( derivatives ),
  stepper_( x_, y_, dydx_ )
{}

template< typename Stepper >
void OdeIntegrator< Stepper >::integrate( const double deltaX )
{
  // TODO: inte så numeriskt stabilt:
  const size_t nbOfWholeSteps = size_t( floor( deltaX / dx_ ) );
  const double remainingStep = deltaX - nbOfWholeSteps * dx_;
  assert( remainingStep >= 0.0 );
  size_t step = 0;
  if( denseOutput_ )
  {
    size_t outputSize = remainingStep > 1.0e-10 ? nbOfWholeSteps + 2 : nbOfWholeSteps + 1;
    xDense_.resize( outputSize );
    yDense_.resize( outputSize );
    xDense_[0] = x_;
    yDense_[0] = y_;

    while( step < nbOfWholeSteps )
    {
      stepper_.step( dx_, derivatives_ );
      xDense_[step + 1] = x_;
      yDense_[step + 1] = y_;
      ++step;
    }

    if( remainingStep > 0.0 )
    {
      stepper_.step( remainingStep, derivatives_ );
      xDense_[step + 1] = x_;
      yDense_[step + 1] = y_;
    }
  }
  else {
    while( step < nbOfWholeSteps )
    {
      stepper_.step( dx_, derivatives_ );
      ++step;
    }

    if( remainingStep > 1.0e-10 )
      stepper_.step( remainingStep, derivatives_ );
  }
}

template< typename Stepper >
template< typename Type >
void OdeIntegrator< Stepper >::integrateUntil( Type& criterion )
{
  if( denseOutput_ )
  {
    size_t outputSize = 1000; // TODODODO!!!!
    xDense_.reserve( outputSize );
    yDense_.reserve( outputSize );

    xDense_.push_back( x_ );
    yDense_.push_back( y_ );

    do {
      stepper_.step( dx_, derivatives_ );
      xDense_.push_back( x_ );
      yDense_.push_back( y_ );
    } while( criterion( *(xDense_.rbegin() + 1), *(yDense_.rbegin() + 1), x_, y_ ) );
  }
  else {
    double xPrev = x_;
    vector< double > yPrev = y_;
    do {
      stepper_.step( dx_, derivatives_ );
    } while( criterion( xPrev, yPrev, x_, y_ ) );
  }
}