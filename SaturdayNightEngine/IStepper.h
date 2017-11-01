#pragma once
#include< vector >

// Interface class for ODE integrators
class IStepper
{
public:
  IStepper( double& x,
            std::vector< double >& y,
            std::vector< double >& dydx );

protected:
  double& x_;
  std::vector< double >& y_;
  std::vector< double >& dydx_;
};

