#pragma once
#include "stdafx.h"

#include "MathUtil.h"
#include "Physics.h"
#include "PhysicalConstants.h"

#include <gtest/gtest.h>

enum class IntegrationMethod { EulerForward, SymplecticEuler, RK4 };

class PhysicsTests : public ::testing::TestWithParam< IntegrationMethod >
{
public:
  //static void compileMe()
  //{
  //  std::cout << "Mumrik" << std::endl;
  //}

  //static void test();
private:
};

//------- Functors for stop criteria

struct AngleIsLessThanTwoPi {
  AngleIsLessThanTwoPi():
    totalAngle_( 0.0 )
  {}

  bool operator() ( const double xPrev, 
                    const std::vector< double >& yPrev, 
                    const double x, 
                    const std::vector< double >& y )
  {
    // TODO, varför funkar inte detta? Räknar fel på någe sätt...
    //const auto previousSunPosition = Eigen::Vector2d( yPrev[0], yPrev[1] );
    //const auto sunPosition = Eigen::Vector2d( y[0], y[1] );
    //const auto previousPlanetPosition = Eigen::Vector2d( yPrev[4], yPrev[5] );
    //const auto planetPosition = Eigen::Vector2d( y[4], y[5] );

    //const auto rPrev = previousPlanetPosition - previousSunPosition;
    //const auto r = planetPosition - sunPosition;
    //const double dAngle = MathUtil::angleBetween( r, rPrev );
    //totalAngle_ += dAngle;

    const double rx = y[0] - y[4];
    const double ry = y[1] - y[5];
    const double rxPrev = yPrev[0] - yPrev[4];
    const double ryPrev = yPrev[1] - yPrev[5];

    const double dotProduct = rx * rxPrev + ry * ryPrev;
    const double normPrev = sqrt( rxPrev * rxPrev + ryPrev * ryPrev );
    const double norm = sqrt( rx * rx + ry * ry );
    const double dAngle = acos( dotProduct / (norm * normPrev) );
    totalAngle_ += dAngle;

    return totalAngle_ < 2.0 * M_PI;
  }

private:
  double totalAngle_;
};

struct TotalEnergyIsConstant {
  TotalEnergyIsConstant( const std::vector< double >& masses,
                         const std::vector< double >& y,
                         const double fraction,
                         const double xMax ):
    masses_( masses ),
    fraction_( fraction ),
    xMax_( xMax )
  {
    initialEnergy_ = Physics::computeTotalEnergy2D( y, masses_ );
  }

  bool operator() ( const double xPrev,
                    const std::vector< double >& yPrev,
                    const double x,
                    const std::vector< double >& y )
  {
    return x < xMax_ && abs( Physics::computeTotalEnergy2D( y, masses_ ) - initialEnergy_ ) < abs( initialEnergy_ * fraction_ );
  }

private:
  const std::vector< double >& masses_; // TODO
  double initialEnergy_;
  double fraction_;
  double xMax_;
};