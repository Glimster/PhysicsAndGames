#pragma once
#include "NR_StepperBase.h"
#include "MathUtil.h"


template <class D>
struct NR_StepperDopr5 : NR_StepperBase {
  typedef D Dtype;
  std::vector< double > k2, k3, k4, k5, k6;
  std::vector< double > rcont1, rcont2, rcont3, rcont4, rcont5;
  std::vector< double > dydxnew;
  NR_StepperDopr5( std::vector< double > &yy, 
                   std::vector< double > &dydxx, 
                   double &xx,
                   const double atoll, 
                   const double rtoll, 
                   bool dens );
  void step( const double htry, D &derivs );
  void dy( const double h, D &derivs );
  void prepare_dense( const double h, D &derivs );
  double dense_out( const int i, const double x, const double h );
  double error();
  struct Controller {
    double hnext, errold;
    bool reject;
    Controller();
    bool success( const double err, double &h );
  };
  Controller con;
};

template <class D>
NR_StepperDopr5<D>::NR_StepperDopr5( std::vector< double > &yy, std::vector< double > &dydxx, double &xx,
                               const double atoll, const double rtoll, bool dens ) :
  NR_StepperBase( yy, dydxx, xx, atoll, rtoll, dens ), 
  k2( n_ ), 
  k3( n_ ), 
  k4( n_ ), 
  k5( n_ ), 
  k6( n_ ),
  rcont1( n_ ), 
  rcont2( n_ ), 
  rcont3( n_ ), 
  rcont4( n_ ), 
  rcont5( n_ ), 
  dydxnew( n_ ) 
{
  //EPS = numeric_limits<double>::epsilon();
}

template <class D>
NR_StepperDopr5<D>::Controller::Controller() : reject( false ), errold( 1.0e-4 ) {}
template <class D>
bool NR_StepperDopr5<D>::Controller::success( const double err, double &h ) {
  static const double beta = 0.0, alpha = 0.2 - beta*0.75, safe = 0.9, minscale = 0.2,
    maxscale = 10.0;
  double scale;
  if( err <= 1.0 ) {
    if( err == 0.0 )
      scale = maxscale;
    else {
      scale = safe*pow( err, -alpha )*pow( errold, beta );
      if( scale<minscale ) scale = minscale;
      if( scale>maxscale ) scale = maxscale;
    }
    if( reject )
      hnext = h*MathUtil::min( scale, 1.0 );
    else
      hnext = h*scale;
    errold = MathUtil::max( err, 1.0e-4 );
    reject = false;
    return true;
  }
  else {
    scale = MathUtil::max( safe*pow( err, -alpha ), minscale );
    h *= scale;
    reject = true;
    return false;
  }
}
template <class D>
double NR_StepperDopr5<D>::dense_out( const int i, const double x, const double h ) {
  double s = (x - xOld_) / h;
  double s1 = 1.0 - s;
  return rcont1[i] + s*(rcont2[i] + s1*(rcont3[i] + s*(rcont4[i] + s1*rcont5[i])));
}
template <class D>
void NR_StepperDopr5<D>::dy( const double h, D &derivs ) {
  static const double c2 = 0.2, c3 = 0.3, c4 = 0.8, c5 = 8.0 / 9.0, a21 = 0.2, a31 = 3.0 / 40.0,
    a32 = 9.0 / 40.0, a41 = 44.0 / 45.0, a42 = -56.0 / 15.0, a43 = 32.0 / 9.0, a51 = 19372.0 / 6561.0,
    a52 = -25360.0 / 2187.0, a53 = 64448.0 / 6561.0, a54 = -212.0 / 729.0, a61 = 9017.0 / 3168.0,
    a62 = -355.0 / 33.0, a63 = 46732.0 / 5247.0, a64 = 49.0 / 176.0, a65 = -5103.0 / 18656.0,
    a71 = 35.0 / 384.0, a73 = 500.0 / 1113.0, a74 = 125.0 / 192.0, a75 = -2187.0 / 6784.0,
    a76 = 11.0 / 84.0, e1 = 71.0 / 57600.0, e3 = -71.0 / 16695.0, e4 = 71.0 / 1920.0,
    e5 = -17253.0 / 339200.0, e6 = 22.0 / 525.0, e7 = -1.0 / 40.0;
  std::vector< double > ytemp( n_ );
  int i;
  for( i = 0; i<n_; i++ )
    ytemp[i] = y_[i] + h*a21*dydx_[i];
  derivs( x_ + c2*h, ytemp, k2 );
  for( i = 0; i<n_; i++ )
    ytemp[i] = y_[i] + h*(a31*dydx_[i] + a32*k2[i]);
  derivs( x_ + c3*h, ytemp, k3 );
  for( i = 0; i<n_; i++ )
    ytemp[i] = y_[i] + h*(a41*dydx_[i] + a42*k2[i] + a43*k3[i]);
  derivs( x_ + c4*h, ytemp, k4 );
  for( i = 0; i<n_; i++ )
    ytemp[i] = y_[i] + h*(a51*dydx_[i] + a52*k2[i] + a53*k3[i] + a54*k4[i]);
  derivs( x_ + c5*h, ytemp, k5 );
  for( i = 0; i<n_; i++ )
    ytemp[i] = y_[i] + h*(a61*dydx_[i] + a62*k2[i] + a63*k3[i] + a64*k4[i] + a65*k5[i]);
  double xph = x_ + h;
  derivs( xph, ytemp, k6 );
  for( i = 0; i<n_; i++ )
    yout_[i] = y_[i] + h*(a71*dydx_[i] + a73*k3[i] + a74*k4[i] + a75*k5[i] + a76*k6[i]);
  derivs( xph, yout_, dydxnew );
  for( i = 0; i<n_; i++ ) {
    yerr_[i] = h*(e1*dydx_[i] + e3*k3[i] + e4*k4[i] + e5*k5[i] + e6*k6[i] + e7*dydxnew[i]);
  }
}
template <class D>
double NR_StepperDopr5<D>::error() {
  double err = 0.0, sk;
  for( int i = 0; i<n_; i++ ) {
    sk = absoluteTolerance_ + relativeTolerance_*MathUtil::max( abs( y_[i] ), abs( yout_[i] ) );
    err += SQR( yerr_[i] / sk );
  }
  return sqrt( err / n_ );
}
template <class D>
void NR_StepperDopr5<D>::prepare_dense( const double h, D &derivs ) {
  std::vector< double > ytemp( n_ );
  static const double d1 = -12715105075.0 / 11282082432.0,
    d3 = 87487479700.0 / 32700410799.0, d4 = -10690763975.0 / 1880347072.0,
    d5 = 701980252875.0 / 199316789632.0, d6 = -1453857185.0 / 822651844.0,
    d7 = 69997945.0 / 29380423.0;
  for( int i = 0; i<n_; i++ ) {
    rcont1[i] = y_[i];
    double ydiff = yout_[i] - y_[i];
    rcont2[i] = ydiff;
    double bspl = h*dydx_[i] - ydiff;
    rcont3[i] = bspl;
    rcont4[i] = ydiff - h*dydxnew[i] - bspl;
    rcont5[i] = h*(d1*dydx_[i] + d3*k3[i] + d4*k4[i] + d5*k5[i] + d6*k6[i] +
                    d7*dydxnew[i]);
  }
}
template <class D>
void NR_StepperDopr5<D>::step( const double htry, D &derivs ) {
  double h = htry;
  for( ;;) {
    dy( h, derivs );
    double err = error();
    if( con.success( err, h ) ) break;
    if( abs( h ) <= abs( x_ )*numeric_limits<double>::epsilon() )
      throw("stepsize underflow in NR_StepperDopr5");
  }
  if( dense_ )
    prepare_dense( h, derivs );
  dydx_ = dydxnew;
  y_ = yout_;
  xOld_ = x_;
  x_ += (hdid_ = h);
  hnext_ = con.hnext;
}
