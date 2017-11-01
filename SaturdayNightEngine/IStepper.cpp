#include "stdafx.h"
#include "IStepper.h"

using namespace std;
IStepper::IStepper( double& x,
                    vector< double >& y,
                    vector< double >& dydx ):
  x_( x ),
  y_( y ),
  dydx_( dydx )
{}