#include "stdafx.h"
#include "NR_StepperBase.h"

using namespace std;
NR_StepperBase::NR_StepperBase( vector< double >& y,
                                vector< double >& dydx,
                                double& x,
                                const double absoluteTolerance,
                                const double relativeTolerance,
                                bool dense ) :
  x_( x ),
  xOld_(),
  y_( y ),
  dydx_( dydx ),
  absoluteTolerance_( absoluteTolerance ),
  relativeTolerance_( relativeTolerance ),
  dense_( dense ),
  hdid_(),
  hnext_(),
  n_( int( y_.size() ) ),
  neqn_( n_ ),
  yout_( n_ ),
  yerr_( n_ )
{}