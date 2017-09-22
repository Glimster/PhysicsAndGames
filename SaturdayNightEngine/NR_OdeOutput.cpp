#include "stdafx.h"
#include "NR_OdeOutput.h"

using namespace std;

NR_OdeOutput::NR_OdeOutput():
  kmax_( -1 ),
  dense_( false ),
  count_( 0 )
{}

NR_OdeOutput::NR_OdeOutput( const int nSaved ) :
  kmax_( 500 ),
  nSaved_( nSaved ),
  count_( 0 ),
  xsave_( kmax_ )
{
  dense_ = nSaved_ > 0 ? true : false;
}

void NR_OdeOutput::init( const int neqn,
                      const double x1,
                      const double x2 )
{
  nvar_ = neqn;
  if( kmax_ == -1 )
    return; // Why?

  ysave_.resize( nvar_, kmax_ );
  if( dense_ )
  {
    x1_ = x1;
    x2_ = x2;
    xout_ = x1;
    dxout_ = (x2 - x1) / nSaved_;
  } 
}

void NR_OdeOutput::resize()
{
  // TODO, borde kunna göras med stl
  int kold = kmax_;
  kmax_ *= 2;
  vector< double > tempvec( xsave_ );
  xsave_.resize( kmax_ );
  for( size_t k = 0; k < kold; ++k )
    xsave_[k] = tempvec[k];
 
  // Och det här med Eigen
  Eigen::MatrixXd tempmat( ysave_ );
  ysave_.resize( nvar_, kmax_ );
  for( size_t i = 0; i < nvar_; ++i )
    for( size_t k = 0; k < kold; ++k )
      ysave_( i, k ) = tempmat( i, k );
}

//template< typename Stepper >
//void NR_OdeOutput::saveDense( Stepper& s, const double xout, const double h )
//{
//  if( count_ == kmax_ )
//    resize();
//  for( size_t i = 0; i < nvar_; ++i )
//    ysave_( i, count_ ) = s.denseOut( i, xout, h );
//  xsave_[count_++] = xout;
//}

void NR_OdeOutput::save( const double x, const std::vector< double >& y )
{
  if( kmax_ < 0.0 )
    return;
  if( count_ == kmax_ )
    resize();

  for( size_t i = 0; i < nvar_; ++i )
    ysave_( i, count_ ) = y[i];
  xsave_[count_++] = x;
}

//template< typename Stepper >
//void NR_OdeOutput::out( const int nstp, const double x, const vector< double >& y, Stepper& s, const double h )
//{
//  if( !dense_ )
//    throw runtime_error( "Dense output not set in NR_OdeOutput" );
//  if( nstp == -1 )
//  {
//    save( x, y );
//    xout_ += dxout_;
//  }
//  else
//  {
//    while( (x - xout_) * (x2_ - x1_) > 0.0 )
//    {
//      saveDense( s, xout_, h );
//      xout_ += dxout_;
//    }
//  }
//}