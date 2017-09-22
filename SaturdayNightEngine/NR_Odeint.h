#pragma once
#include "NR_OdeOutput.h"

template<class Stepper>
struct NR_Odeint
{
  static const int MAXSTP = 50000;
  double EPS;
  int nok;
  int nbad;
  int nvar;
  double x1, x2, hmin;
  bool dense;
  std::vector< double > y, dydx;
  std::vector< double >& ystart;
  NR_OdeOutput& out;
  typename Stepper::Dtype& derivs;
  Stepper s;
  int nstp;
  double x, h;
  NR_Odeint( std::vector< double > &ystartt, 
             const double xx1, 
             const double xx2,
             const double atol, 
             const double rtol, 
             const double h1,
             const double hminn, 
             NR_OdeOutput &outt, 
             typename Stepper::Dtype &derivss );
  void integrate();

  template<class T>
  inline T SIGN( const T &a, const T &b )
  {
    return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
  }
};

template<class Stepper>
NR_Odeint<Stepper>::NR_Odeint( std::vector< double > &ystartt,
                               const double xx1, 
                               const double xx2,
                               const double atol, 
                               const double rtol,
                               const double h1, 
                               const double hminn,
                               NR_OdeOutput &outt, 
                               typename Stepper::Dtype &derivss ) :
  nvar( int( ystartt.size() ) ),
  y( nvar ), 
  dydx( nvar ), 
  ystart( ystartt ), 
  x( xx1 ), 
  nok( 0 ), 
  nbad( 0 ),
  x1( xx1 ), 
  x2( xx2 ), 
  hmin( hminn ), 
  dense( outt.dense_ ), 
  out( outt ), 
  derivs( derivss ),
  s( y, dydx, x, atol, rtol, dense )
{
  EPS = numeric_limits<double>::epsilon();
  h = SIGN( h1, x2 - x1 );
  for( int i = 0; i < nvar; i++ ) y[i] = ystart[i];

  out.init( s.neqn_, x1, x2 );
}

template<class Stepper>
void NR_Odeint<Stepper>::integrate()
{
  derivs( x, y, dydx );
  if( dense )
    out.out( -1, x, y, s, h );
  else
    out.save( x, y );
  for( nstp = 0; nstp < MAXSTP; nstp++ )
  {
    if( (x + h * 1.0001 - x2) * (x2 - x1) > 0.0 )
      h = x2 - x;
    s.step( h, derivs );
    if( s.hdid_ == h ) ++nok; else ++nbad;
    if( dense )
      out.out( nstp, x, y, s, s.hdid_ );
    else
      out.save( x, y );
    if( (x - x2) * (x2 - x1) >= 0.0 )
    {
      for( int i = 0; i < nvar; i++ ) ystart[i] = y[i];
      if( out.kmax_ > 0 && abs( out.xsave_[out.count_ - 1] - x2 ) > 100.0 * abs( x2 )*EPS )
        out.save( x, y );
      return;
    }
    if( abs( s.hnext_ ) <= hmin ) throw("Step size too small in NR_Odeint");
    h = s.hnext_;
  }
  throw ("Too many steps in routine NR_Odeint");
}