#pragma once
#include< vector >

class NR_StepperBase
{
public:
  NR_StepperBase( std::vector< double >& y,
                  std::vector< double >& dydx,
                  double& x,
                  const double absoluteTolerance,
                  const double relativeTolerance,
                  bool dense );

  template<class T>
  inline T SQR( const T a ) { return a*a; }

public:
  double& x_;
  double xOld_;
  std::vector< double >& y_;
  std::vector< double >& dydx_;
  double absoluteTolerance_, relativeTolerance_;
  bool dense_;
  double hdid_;
  double hnext_;
  int n_, neqn_;
  std::vector< double > yout_, yerr_;
};

