#include "stdafx.h"
#include "PhysicsTests.h"

#include "ODEs.h"
#include "OdeIntegrator.h"
#include "EulerMethods.h"
#include "RK4.h"

#include "PhysicalData.h"

#include "NR_Odeint.h"
#include "NR_StepperDopr5.h"

#include "Utilities.h"

// TODO, gammalt, ta bort?
#include "MotionManager.h"
#include "PhysicalObject.h"

using namespace std;

INSTANTIATE_TEST_CASE_P( IntegrationTechniques,
                         PhysicsTests,
                         ::testing::Values( IntegrationMethod::EulerForward, IntegrationMethod::SymplecticEuler, IntegrationMethod::RK4 ) );

TEST_F( PhysicsTests, NR_VanDerPolOscillator )
{
  const int nVariables = 2; // y0 and y1 = y0_dot
  const double absoluteTolerance = 1.0e-3;
  const double relativeTolerance = absoluteTolerance;
  const double initialStep = 0.01; // h1
  const double minStep = 0.0; // hmin

  // Integration limits
  const double x1 = 0.0;
  const double x2 = 2.0;

  // Initial values
  vector< double > yStart = { 2.0, 0.0 };
  NR_OdeOutput output( 20 ); // Dense output at 20 points (plus x1)
  VanDerPol derivatives( 1.0e-3 );
  NR_Odeint< NR_StepperDopr5< VanDerPol > > integrator( yStart,
                                                        x1, x2,
                                                        absoluteTolerance,
                                                        relativeTolerance,
                                                        initialStep,
                                                        minStep,
                                                        output,
                                                        derivatives );
  integrator.integrate();

  for( size_t i = 0; i < output.count_; ++i )
    cout << "x[i] = " << output.xsave_[i] << ", y[0, i] = " << output.ysave_( 0, i ) << ", y[1, i] = " << output.ysave_( 1, i ) << endl;

  // TODO, ska det vara nåt test av detta (bara taget från NR)?
}

TEST_F( PhysicsTests, FallingObject )
{
  // Initial state
  const double px = 0.0, py = 0.0;
  const double vx = 10.0, vy = 5.0;
  const vector< double > y0 = { px, py, vx, vy };
  const double t0 = 0.0, t1 = 5.0;

  const double absoluteTolerance = 1.0e-3;

  const double g = 9.82;
  FixedGravitationalField derivatives( g );

  // Analytical solution
  const double pxAnalytic = vx * (t1 - t0) + px;
  const double pyAnalytic = 0.5 * (-g) * (t1 - t0) * (t1 - t0) + vy * (t1 - t0) + py;
  const double vxAnalytic = vx;
  const double vyAnalytic = vy + (-g) * (t1 - t0);
  cout << "Analytic solution:" << endl;
  cout << "Final position = (" << pxAnalytic << ", " << pyAnalytic << ")" << endl;
  cout << "Final velocity = (" << vxAnalytic << ", " << vyAnalytic << ")" << endl;

  {
    cout << endl << "--------StepperDopr5" << endl;
    vector< double > yStart( y0 );
    const double relativeTolerance = absoluteTolerance;
    const double initialStep = 0.1;
    const double minStep = 0.0;
    NR_OdeOutput output( -1 );
    NR_Odeint< NR_StepperDopr5< FixedGravitationalField > > integrator( yStart,
                                                                        t0, t1,
                                                                        absoluteTolerance,
                                                                        relativeTolerance,
                                                                        initialStep,
                                                                        minStep,
                                                                        output,
                                                                        derivatives );
    integrator.integrate();

    cout << "Final" << endl;
    cout << "Pos = " << output.ysave_( 0, output.count_ - 1 ) << ", " << output.ysave_( 1, output.count_ - 1 ) << endl;
    cout << "Vel = " << output.ysave_( 2, output.count_ - 1 ) << ", " << output.ysave_( 3, output.count_ - 1 ) << endl;

    EXPECT_NEAR( pxAnalytic, output.ysave_( 0, output.count_ - 1 ), absoluteTolerance );
    EXPECT_NEAR( pyAnalytic, output.ysave_( 1, output.count_ - 1 ), absoluteTolerance );
    EXPECT_NEAR( vxAnalytic, output.ysave_( 2, output.count_ - 1 ), absoluteTolerance );
    EXPECT_NEAR( vyAnalytic, output.ysave_( 3, output.count_ - 1 ), absoluteTolerance );
  }

  {
    cout << endl << "----------Euler forward" << endl;
    const double dt = 0.1; // Should fail miserably
    double time = t0;
    vector< double > y( y0 );
    OdeIntegrator< EulerForward< FixedGravitationalField > > integrator( time,
                                                                         y,
                                                                         dt,
                                                                         false,
                                                                         derivatives );
    integrator.integrate( t1 - time );

    cout << "Time = " << time << endl;
    cout << "Final position = " << y[0] << ", " << y[1] << endl;
    cout << "Final velocity = " << y[2] << ", " << y[3] << endl;

    EXPECT_NEAR( time, t1, 1.0e-10 );
    EXPECT_NEAR( pxAnalytic, y[0], absoluteTolerance );
    EXPECT_GT( fabs( pyAnalytic - y[1] ), absoluteTolerance );
    EXPECT_NEAR( vxAnalytic, y[2], absoluteTolerance );
    EXPECT_NEAR( vyAnalytic, y[3], absoluteTolerance );
  }

  { // RK4
    cout << endl << "----------RK4" << endl;
    const double dt = 1.0;
    double time = t0;
    vector< double > y( y0 );
    OdeIntegrator< RK4< FixedGravitationalField > > integrator( time,
                                                                y,
                                                                dt,
                                                                false,
                                                                derivatives );
    integrator.integrate( t1 - time );
    cout << "Time = " << time << endl;
    cout << "Final position = " << y[0] << ", " << y[1] << endl;
    cout << "Final velocity = " << y[2] << ", " << y[3] << endl;

    EXPECT_NEAR( time, t1, 1.0e-10 );
    EXPECT_NEAR( pxAnalytic, y[0], absoluteTolerance );
    EXPECT_NEAR( pyAnalytic, y[1], absoluteTolerance );
    EXPECT_NEAR( vxAnalytic, y[2], absoluteTolerance );
    EXPECT_NEAR( vyAnalytic, y[3], absoluteTolerance );
  }
}

TEST_F( PhysicsTests, SimpleHarmonicOscillator )
{
  const double k = 27.0; // Spring constant
  const double m = 2.1; // kg

  // Initial conditions
  const double x0 = -5.0;
  const double v0 = 3.0;
  const double t0 = 0.0;

  const double absoluteTolerance = 1.0e-3;

  // Analytic solution
  const double frequency = 1 / (2.0 * M_PI) * sqrt( k / m );
  const double period = 1.0 / frequency;
  auto hoAmplitude = []( double f, double x, double v ) {
    if( v == 0.0 )
      return fabs( x );
    const double angularFrequence = 2.0 * M_PI * f;
    const double phase = atan( angularFrequence * x / v );
    return x / sin( phase );
  };

  cout << "Analytic solution:" << endl;
  cout << "Frequency = " << frequency << endl;
  cout << "Period = " << 1.0 / frequency << endl;
  cout << "Amplitude = " << hoAmplitude( frequency, x0, v0 ) << endl;

  HookesLaw derivatives( k, m );
  {
    cout << endl << "----------Euler forward" << endl;
    const double dt = 0.00001; // TODO: Så mycket krävs...
    double time = t0;
    vector< double > y = { x0, v0 };
    OdeIntegrator< EulerForward< HookesLaw > > integrator( time,
                                                           y,
                                                           dt,
                                                           true,
                                                           derivatives );
    integrator.integrate( period );

    const auto& xDense = integrator.getXDense();
    const auto& yDense = integrator.getYDense();
    auto minmax = minmax_element( yDense.begin(), yDense.end(), []( const vector< double >& a, const vector< double >& b ) { return a[0] < b[0]; } );
    cout << "Time = " << time << endl;
    cout << "Final position = " << y[0] << endl;
    cout << "Final velocity = " << y[1] << endl;
    cout << "Min x = " << (*minmax.first)[0] << endl;
    cout << "Max x = " << (*minmax.second)[0] << endl;

    EXPECT_NEAR( time, period, 1.0e-10 );
    EXPECT_NEAR( y[0], x0, absoluteTolerance );
    EXPECT_NEAR( y[1], v0, absoluteTolerance );
    EXPECT_NEAR( fabs( (*minmax.first)[0] ), fabs( (*minmax.second)[0] ), absoluteTolerance );
    EXPECT_NEAR( fabs( (*minmax.first)[0] ), hoAmplitude( frequency, x0, v0 ), absoluteTolerance );
  }
  {
    cout << endl << "----------Symplectic Euler" << endl;
    const double dt = 0.0001; // TODO: Så mycket krävs...
    double time = t0;
    vector< double > y = { x0, v0 };
    OdeIntegrator< SymplecticEuler1D< HookesLaw > > integrator( time,
                                                                y,
                                                                dt,
                                                                true,
                                                                derivatives );
    integrator.integrate( period );

    const auto& xDense = integrator.getXDense();
    const auto& yDense = integrator.getYDense();
    auto minmax = minmax_element( yDense.begin(), yDense.end(), []( const vector< double >& a, const vector< double >& b ) { return a[0] < b[0]; } );
    cout << "Time = " << time << endl;
    cout << "Final position = " << y[0] << endl;
    cout << "Final velocity = " << y[1] << endl;
    cout << "Min x = " << (*minmax.first)[0] << endl;
    cout << "Max x = " << (*minmax.second)[0] << endl;

    EXPECT_NEAR( time, period, 1.0e-10 );
    EXPECT_NEAR( y[0], x0, absoluteTolerance );
    EXPECT_NEAR( y[1], v0, absoluteTolerance );
    EXPECT_NEAR( fabs( (*minmax.first)[0] ), fabs( (*minmax.second)[0] ), absoluteTolerance );
    EXPECT_NEAR( fabs( (*minmax.first)[0] ), hoAmplitude( frequency, x0, v0 ), absoluteTolerance );
  }
  {
    cout << endl << "----------RK4" << endl;
    const double dt = 0.01; // TODO: Så mycket krävs...
    double time = t0;
    vector< double > y = { x0, v0 };
    OdeIntegrator< RK4< HookesLaw > > integrator( time,
                                                  y,
                                                  dt,
                                                  true,
                                                  derivatives );
    integrator.integrate( period );

    const auto& xDense = integrator.getXDense();
    const auto& yDense = integrator.getYDense();
    auto minmax = minmax_element( yDense.begin(), yDense.end(), []( const vector< double >& a, const vector< double >& b ) { return a[0] < b[0]; } );
    cout << "Time = " << time << endl;
    cout << "Final position = " << y[0] << endl;
    cout << "Final velocity = " << y[1] << endl;
    cout << "Min x = " << (*minmax.first)[0] << endl;
    cout << "Max x = " << (*minmax.second)[0] << endl;

    EXPECT_NEAR( time, period, 1.0e-10 );
    EXPECT_NEAR( y[0], x0, absoluteTolerance );
    EXPECT_NEAR( y[1], v0, absoluteTolerance );
    EXPECT_NEAR( fabs( (*minmax.first)[0] ), fabs( (*minmax.second)[0] ), absoluteTolerance );
    EXPECT_NEAR( fabs( (*minmax.first)[0] ), hoAmplitude( frequency, x0, v0 ), absoluteTolerance );
  }

  {
    cout << endl << "--------StepperDopr5" << endl;
    vector< double > yStart = { x0, v0 };
    const double relativeTolerance = absoluteTolerance;
    const double initialStep = 0.01;
    const double minStep = 0.0;
    NR_OdeOutput output( 500 );
    // TODO, why is ~1e-2 better tolerance needed to satisfy the tolerance?
    NR_Odeint< NR_StepperDopr5< HookesLaw > > integrator( yStart,
                                                          t0,
                                                          t0 + period,
                                                          absoluteTolerance * 1e-2,
                                                          relativeTolerance * 1e-2,
                                                          initialStep,
                                                          minStep,
                                                          output,
                                                          derivatives );
    integrator.integrate();

    double min = numeric_limits< double >::max();
    double max = numeric_limits< double >::min();
    for( size_t i = 0; i < output.count_; ++i ) {
      //cout << "Pos = " << output.ysave_( 0, i ) << ", " << output.ysave_( 1, i ) << endl;
      min = MathUtil::min( min, output.ysave_( 0, i ) );
      max = MathUtil::max( max, output.ysave_( 0, i ) );
    }
    cout << "Time = " << output.xsave_[output.count_ - 1] << endl;
    cout << "Nok = " << integrator.nok << endl;
    cout << "Nbad = " << integrator.nbad << endl;
    cout << "Final position = " << output.ysave_( 0, output.count_ - 1 ) << endl;
    cout << "Final velocity = " << output.ysave_( 1, output.count_ - 1 ) << endl;
    cout << "Min x = " << min << endl;
    cout << "Max x = " << max << endl;
    cout << "*Eigen* Min x (why is this value not correct?) = " << output.ysave_.rowwise().minCoeff()[0] << endl;
    cout << "*Eigen* Max x = " << output.ysave_.rowwise().maxCoeff()[0] << endl;

    EXPECT_NEAR( output.xsave_[output.count_ - 1], period, 1.0e-10 );
    EXPECT_NEAR( output.ysave_( 0, output.count_ - 1 ), x0, absoluteTolerance );
    EXPECT_NEAR( output.ysave_( 1, output.count_ - 1 ), v0, absoluteTolerance );
    EXPECT_NEAR( fabs( min ), fabs( max ), absoluteTolerance );
    EXPECT_NEAR( fabs( min ), hoAmplitude( frequency, x0, v0 ), absoluteTolerance );
  }
}

// Tests that a (light) planet orbiting around a (heavy) sun follows Kepler's
// laws of planetary motion
// https://en.wikipedia.org/wiki/Kepler%27s_laws_of_planetary_motion
//
// Fails for forward Euler
// Succeeds for RK4
TEST_P( PhysicsTests, KeplersLaws )
{
  // Setup
  auto spaceObjectData = PhysicalData::setupHeavySunLightPlanet();
  assert( spaceObjectData.size() == 2 );
  const auto& sun = spaceObjectData[0];
  const auto& planet = spaceObjectData[1];
  
  vector< double > masses = { sun.mass, planet.mass };
  GravitationalNBody derivatives( masses );

  const double t0 = 0.0;
  const double dt = 0.001;
  double time = t0;
  vector< double > y = {
    sun.position(0),
    sun.position(1),
    sun.velocity(0),
    sun.velocity(1),
    planet.position(0),
    planet.position(1),
    planet.velocity(0),
    planet.velocity(1),
  };

  const double tolerance = 1.0e-3;

  // Act
  vector< double > timeDense;
  vector< vector< double > > yDense;

  const auto& integrationMethod = GetParam();
  switch( integrationMethod ) {
  case IntegrationMethod::EulerForward:
  {
    cout << endl << "----------Euler Forward" << endl;
    OdeIntegrator< EulerForward< GravitationalNBody > > integrator( time,
                                                                    y,
                                                                    dt,
                                                                    true,
                                                                    derivatives );
    AngleIsLessThanTwoPi angleIsLessThanTwoPi;
    integrator.integrateUntil( angleIsLessThanTwoPi );
    timeDense = integrator.getXDense();
    yDense = integrator.getYDense();
    break;
  }
  case IntegrationMethod::SymplecticEuler:
  {
    cout << endl << "----------Symplectic Euler" << endl;
    OdeIntegrator< SymplecticEuler2D< GravitationalNBody > > integrator( time,
                                                                         y,
                                                                         dt,
                                                                         true,
                                                                         derivatives );
    AngleIsLessThanTwoPi angleIsLessThanTwoPi;
    integrator.integrateUntil( angleIsLessThanTwoPi );
    timeDense = integrator.getXDense();
    yDense = integrator.getYDense();
    break;
  }
  case IntegrationMethod::RK4:
  {
    cout << endl << "----------Runge Kutta 4" << endl;
    OdeIntegrator< RK4< GravitationalNBody > > integrator( time,
                                                           y,
                                                           dt,
                                                           true,
                                                           derivatives );
    AngleIsLessThanTwoPi angleIsLessThanTwoPi;
    integrator.integrateUntil( angleIsLessThanTwoPi );
    timeDense = integrator.getXDense();
    yDense = integrator.getYDense();
    break;
  }
  default:
    throw runtime_error( "Unknown integration technique" );
  }

  const auto finalSunPosition = Eigen::Vector2d( y[0], y[1] );
  const auto finalPlanetPosition = Eigen::Vector2d( y[4], y[5] );
  cout << "Final position planet = " << Utilities::toString( finalPlanetPosition ) << endl;
  cout << "Final position sun = " << Utilities::toString( finalSunPosition ) << endl;

  if( integrationMethod == IntegrationMethod::EulerForward ) {
    cout << "Expected to fail!" << endl;
    EXPECT_GT( (planet.position.cast< double >() - finalPlanetPosition).norm(), tolerance );
    return;
  }

  EXPECT_NEAR( (sun.position.cast< double >() - finalSunPosition).norm(), 0.0, tolerance );
  EXPECT_NEAR( (planet.position.cast< double >() - finalPlanetPosition).norm(), 0.0, tolerance );

  // Analyze
  double distanceTravelled = 0.0;
  double totalArea = 0.0;
  double totalAngle = 0.0;

  vector< double > sweptAreas;

  Eigen::Vector2d minPosition( planet.position.cast< double >() ); // TODO, ta bort floats!
  Eigen::Vector2d maxPosition( planet.position.cast< double >() );
  // Distances to sun
  double dPericenter = (planet.position - sun.position).norm();
  double dApocenter = (planet.position - sun.position).norm();
  double vPericenter = planet.velocity.norm();
  double vApocenter = planet.velocity.norm();
  Eigen::Vector2d pericenter( planet.position.cast< double >() ); // TODO, ta bort floats!
  Eigen::Vector2d apocenter( planet.position.cast< double >() );

  double sweptArea = 0.0;
  const size_t segmentAreaStepInterval = timeDense.size() / 10;
  for( size_t i = 1; i < timeDense.size(); ++i )
  {
    const auto previousSunPosition = Eigen::Vector2d( yDense[i - 1][0], yDense[i - 1][1] );
    const auto sunPosition = Eigen::Vector2d( yDense[i][0], yDense[i][1] );
    const auto previousPlanetPosition = Eigen::Vector2d( yDense[i - 1][4], yDense[i - 1][5] );
    const auto planetPosition = Eigen::Vector2d( yDense[i][4], yDense[i][5] );
    const auto planetVelocity = Eigen::Vector2d( yDense[i][6], yDense[i][7] );

    const auto rPrev = previousPlanetPosition - previousSunPosition;
    const auto r = planetPosition - sunPosition;
    const auto d = (planetPosition - previousPlanetPosition).norm();

    distanceTravelled += d;
    totalArea += MathUtil::cross2D( r, rPrev ) / 2.0; // Area of triangle
    sweptArea += MathUtil::cross2D( r, rPrev ) / 2.0; // Area of triangle
    if( i % segmentAreaStepInterval == 0 )
    {
      sweptAreas.push_back( sweptArea );
      sweptArea = 0.0f;
    }
    const double dAngle = MathUtil::angleBetween( r, rPrev );
    totalAngle += dAngle;

    minPosition(0) = min( planetPosition(0), minPosition(0) );
    minPosition(1) = min( planetPosition(1), minPosition(1) );
    maxPosition(0) = max( planetPosition(0), maxPosition(0) );
    maxPosition(1) = max( planetPosition(1), maxPosition(1) );

    if( r.norm() < dPericenter )
    {
      pericenter = planetPosition;
      vPericenter = planetVelocity.norm();
      dPericenter = r.norm();
    }
    if( r.norm() > dApocenter )
    {
      apocenter = planetPosition;
      vApocenter = planetVelocity.norm();
      dApocenter = r.norm();
    }
  }

  // -------------------------------------------
  cout << endl << "Kepler 1: The orbit of a planet is an ellipse with the Sun at one of the two foci." << endl;
  const auto centerPosition = (maxPosition + minPosition) / 2.0f;
  const auto axes = maxPosition - minPosition;
  const double semiMajorAxis = axes.maxCoeff() / 2.0f;
  const double semiMinorAxis = axes.minCoeff() / 2.0f;
  cout << "Semi major axis = " << semiMajorAxis << endl;
  cout << "Semi minor axis = " << semiMinorAxis << endl;

  cout << "pericenter =  " << Utilities::toString( pericenter ) << endl;
  cout << "apocenter =  " << Utilities::toString( apocenter ) << endl;

  // According to https://en.wikipedia.org/wiki/Apsis
  const double eccentricity = (dApocenter - dPericenter) / (dApocenter + dPericenter);
  const double standardGravitationalParameter = Phys::PhysicalConstants::kSquared * sun.mass;
  cout << "Distance at pericenter = " << dPericenter << endl;
  cout << "Speed at pericenter = " << vPericenter << endl;
  cout << "Distance at apocenter = " << dApocenter << endl;
  cout << "Speed at apocenter = " << vApocenter << endl;
  EXPECT_NEAR( dPericenter, (1 - eccentricity) * semiMajorAxis, tolerance );
  EXPECT_NEAR( vPericenter, sqrt( (1 + eccentricity) * standardGravitationalParameter / ((1 - eccentricity) * semiMajorAxis) ), tolerance );
  EXPECT_NEAR( dApocenter, (1 + eccentricity) * semiMajorAxis, tolerance );
  EXPECT_NEAR( vApocenter, sqrt( (1 - eccentricity) * standardGravitationalParameter / ((1 + eccentricity) * semiMajorAxis) ), tolerance );

  const double f = double( MathUtil::fociCenterDistanceEllipse( float( semiMajorAxis ), float( semiMinorAxis ) ) ); // TODO, styr upp float/double
  const Eigen::Vector2d F1( centerPosition(0) - f, centerPosition(1) );
  const Eigen::Vector2d F2( centerPosition(0) + f, centerPosition(1) );
  cout << "Focal point #1 =  " << Utilities::toString( F1 ) << endl;
  cout << "Focal point #2 (sun) =  " << Utilities::toString( F2 ) << endl;
  
  // TODO, styr upp float/double
  EXPECT_NEAR( (finalSunPosition - sun.position.cast< double >()).norm(), 0.0, tolerance ); // Sun "infinitely" heavy
  EXPECT_NEAR( (F2 - finalSunPosition).norm(), 0.0, tolerance ); // Sun at focal point two

  cout << "Distance travelled = " << distanceTravelled << endl;
  EXPECT_NEAR( distanceTravelled, MathUtil::circumferenceEllipse( float( semiMajorAxis ), float( semiMinorAxis ) ), tolerance );

  cout << "Total swept area = " << totalArea << endl;
  EXPECT_NEAR( totalArea, MathUtil::areaEllipse( float( semiMajorAxis ), float( semiMinorAxis ) ), tolerance );

  cout << "Total angle = " << totalAngle << endl;
  EXPECT_NEAR( totalAngle, 2 * M_PI, tolerance );

  // -------------------------------------------
  cout << endl << "Kepler 2: A line joining a planet and the Sun sweeps out equal areas during equal intervals of time." << endl;
  cout << "Swept area in time interval = " << sweptAreas[0] << endl;
  EXPECT_GT( sweptAreas.size(), 1 );
  for( size_t i = 1; i < sweptAreas.size(); ++i )
  {
    EXPECT_NEAR( sweptAreas[i], sweptAreas[0], tolerance );
  }

  // -------------------------------------------
  cout << endl << "Kepler 3: The square of the orbital period of a planet is directly proportional to the cube of the semi-major axis of its orbit." << endl;
  const double constant = 4.0 * M_PI * M_PI / standardGravitationalParameter;
  const double T = sqrt( pow( semiMajorAxis, 3 ) * constant );
  cout << "Orbital period = " << timeDense.back() << endl;
  EXPECT_NEAR( timeDense.back(), T, tolerance );
}

// Tests that the total energy remains constant for a long time
TEST_P( PhysicsTests, StabililtyWRTEnergy )
{
  // Setup
  auto spaceObjectData = PhysicalData::setupHeavySunLightPlanet();
  assert( spaceObjectData.size() == 2 );
  const auto& sun = spaceObjectData[0];
  const auto& planet = spaceObjectData[1];

  vector< double > masses = { sun.mass, planet.mass };
  GravitationalNBody derivatives( masses );

  const double t0 = 0.0;
  const double dt = 0.001;
  const double tMax = 50.0;

  double time = t0;
  vector< double > y = {
    sun.position(0),
    sun.position(1),
    sun.velocity(0),
    sun.velocity(1),
    planet.position(0),
    planet.position(1),
    planet.velocity(0),
    planet.velocity(1),
  };

  const double tolerance = 1.0e-3;

  // Act
  const double fraction = 0.05; // 5% of initial energy
  TotalEnergyIsConstant totalEnergyIsConstant( masses, y, fraction, tMax );
  const double initialEnergy = Physics::computeTotalEnergy2D( y, masses );

  const auto& integrationMethod = GetParam();
  switch( integrationMethod ) {
  case IntegrationMethod::EulerForward:
  {
    // Inexact an unstable (energy increases)
    cout << endl << "----------Euler Forward" << endl;
    OdeIntegrator< EulerForward< GravitationalNBody > > integrator( time,
                                                                    y,
                                                                    dt,
                                                                    false,
                                                                    derivatives );

    integrator.integrateUntil( totalEnergyIsConstant );
    break;
  }
  case IntegrationMethod::SymplecticEuler:
  {
    // Stable, but inexact
    cout << endl << "----------Symplectic Euler" << endl;
    OdeIntegrator< SymplecticEuler2D< GravitationalNBody > > integrator( time,
                                                                         y,
                                                                         dt,
                                                                         true,
                                                                         derivatives );
    integrator.integrateUntil( totalEnergyIsConstant );
    break;
  }
  case IntegrationMethod::RK4:
  {
    // TODO, RK4 actually not stable, energy slowly decreases
    cout << endl << "----------Runge Kutta 4" << endl;
    OdeIntegrator< RK4< GravitationalNBody > > integrator( time,
                                                           y,
                                                           dt,
                                                           true,
                                                           derivatives );
    integrator.integrateUntil( totalEnergyIsConstant );
    break;
  }
  default:
    throw runtime_error( "Unknown integration technique" );
  }

  cout << "Initial energy = " << initialEnergy << endl;
  const double finalEnergy = Physics::computeTotalEnergy2D( y, masses );
  cout << "Time = " << time << endl;
  cout << "Final energy = " << finalEnergy << endl;

  if( integrationMethod == IntegrationMethod::EulerForward )
  {
    cout << "Expected to fail!" << endl;
    EXPECT_LT( time, tMax - dt );
    return;
  }

  EXPECT_NEAR( time, tMax, dt ); // TODO, bör det vara mer exakt än dt?
  EXPECT_NEAR( finalEnergy, initialEnergy, abs( initialEnergy * fraction ) );
}

// TODO: for some reason, the old integration is much faster than the new!?
TEST_F( PhysicsTests, Performance_aLotOfPlanets_oldIntegration )
{
  cout << endl << "--------------------------------" << endl;
  cout << "Test performance for a lot of planets" << endl;
  // Setup
  MotionManager motionManager( MotionManager::Integration::RK4,
                               MotionManager::Technique::standard );

  vector< unique_ptr< PhysicalObject > > celestialBodies;
  vector< PhysicalObject* > physicalObjects;

  auto spaceObjectData = PhysicalData::setupAlotOfPlanets();
  for( auto& sod : spaceObjectData )
  {
    const PhysicalData::PlanetData& data = sod;
    unique_ptr< PhysicalObject > object( new PhysicalObject( data.mass ) );
    object->setPosition( data.position );
    object->setVelocity( data.velocity );

    physicalObjects.push_back( object.get() );
    celestialBodies.push_back( move( object ) );
  }

  auto computeTotalKE = [&]()
  {
    float totalKE = 0.0f;
    for( const auto& sp : physicalObjects )
    {
      totalKE += Physics::kineticEnergy( sp->getMass(), sp->getVelocity() );
    }
    return totalKE;
  };

  auto computeTotalPE = [&]()
  {
    float totalPE = 0.0f;
    for( size_t i = 0; i < physicalObjects.size() - 1; ++i )
    {
      for( size_t j = i + 1; j < physicalObjects.size(); ++j )
      {
        const Eigen::Vector2f r = physicalObjects[i]->getPosition() - physicalObjects[j]->getPosition();
        totalPE += Physics::gravitationalPotentialEnergy( physicalObjects[i]->getMass(), physicalObjects[j]->getMass(), r );
      }
    }
    return totalPE;
  };

  auto initialEnergy = computeTotalPE() + computeTotalKE();
  cout << "Initial energy = " << initialEnergy << endl;

  // Act
  sf::Clock clock;
  float time = 0.0f;
  float dt = 1.0f / 100.0f;
  size_t i = 0;
  while( i < 1000 )
  {
    motionManager.computeLinearDynamics( dt, physicalObjects );
    ++i;
    time += dt;
  }

  cout << "---------------------------------------" << endl;
  cout << "Simulated time = " << time << endl;
  cout << "Simulation time = " << clock.getElapsedTime().asMilliseconds() << " ms" << endl;
  auto finalEnergy = computeTotalPE() + computeTotalKE();
  cout << "Final energy = " << finalEnergy << endl;
  EXPECT_NEAR( initialEnergy, finalEnergy, abs( initialEnergy * 1e-5 ) );
}

// TODO: for some reason, the old integration is much faster than the new!?
TEST_P( PhysicsTests, Performance_aLotOfPlanets )
{
  // Setup
  auto spaceObjectData = PhysicalData::setupAlotOfPlanets();

  const double t0 = 0.0;
  const double dt = 0.01;
  const double t1 = t0 + dt * double( 1000 );

  const size_t nbOfBodies = spaceObjectData.size();
  vector< double > masses( nbOfBodies );
  vector< double > y( nbOfBodies * 4 );
  for( size_t i = 0; i < nbOfBodies; ++i ) {
    const auto& sod = spaceObjectData[i];
    masses[i] = sod.mass;

    size_t j = i * 4;
    y[j] = sod.position(0);
    y[j + 1] = sod.position(1);
    y[j + 2] = sod.velocity(0);
    y[j + 3] = sod.velocity(1);
  }
  auto initialEnergy = Physics::computeTotalEnergy2D( y, masses );

  GravitationalNBody derivatives( masses );

  // Act
  double time = t0;
  sf::Clock clock;
  const auto& integrationMethod = GetParam();
  switch( integrationMethod ) {
  case IntegrationMethod::EulerForward:
  {
    cout << endl << "----------Euler Forward" << endl;
    OdeIntegrator< EulerForward< GravitationalNBody > > integrator( time,
                                                                    y,
                                                                    dt,
                                                                    false,
                                                                    derivatives );
    integrator.integrate( t1 - time );
    break;
  }
  case IntegrationMethod::SymplecticEuler:
  {
    cout << endl << "----------Symplectic Euler" << endl;
    OdeIntegrator< SymplecticEuler2D< GravitationalNBody > > integrator( time,
                                                                         y,
                                                                         dt,
                                                                         false,
                                                                         derivatives );
    integrator.integrate( t1 - time );
    break;
  }
  case IntegrationMethod::RK4:
  {
    cout << endl << "----------Runge Kutta 4" << endl;
    OdeIntegrator< RK4< GravitationalNBody > > integrator( time,
                                                           y,
                                                           dt,
                                                           false,
                                                           derivatives );
    integrator.integrate( t1 - time );
    break;
  }
  default:
    throw runtime_error( "Unknown integration technique" );
  }
  cout << "Initial energy = " << initialEnergy << endl;

  cout << "Number of bodies = " << nbOfBodies << endl;
  cout << "Simulated time = " << time << endl;
  cout << "Simulation time = " << clock.getElapsedTime().asMilliseconds() << " ms" << endl;
  auto finalEnergy = Physics::computeTotalEnergy2D( y, masses );
  cout << "Final energy = " << finalEnergy << endl;
  EXPECT_NEAR( initialEnergy, finalEnergy, abs( initialEnergy * 1e-5 ) );

  // TODO, energy conservation not even close to happen. Investigate!
}

//void PhysicsTests::test() // For performance profiling
TEST_F( PhysicsTests, OldNewIntegrationBenchmark_RK4 )
{
  // Setup
  const double t0 = 0.0;
  const double t1 = 10000.0;
  const double dt = 0.01;

  auto spaceObjectData = PhysicalData::setupRealisticSolarSystem();

  cout << "-------------Old integration" << endl;
  double initialEnergyOld, finalEnergyOld;
  {
    MotionManager motionManager( MotionManager::Integration::RK4,
                                 MotionManager::Technique::functor );

    vector< unique_ptr< PhysicalObject > > celestialBodies;
    vector< PhysicalObject* > physicalObjects;

    for( auto& sod : spaceObjectData )
    {
      const PhysicalData::PlanetData& data = sod;
      unique_ptr< PhysicalObject > object( new PhysicalObject( data.mass ) );
      object->setPosition( data.position );
      object->setVelocity( data.velocity );

      physicalObjects.push_back( object.get() );
      celestialBodies.push_back( move( object ) );
    }

    auto computeTotalKE = [&]()
    {
      float totalKE = 0.0f;
      for( const auto& sp : physicalObjects )
      {
        totalKE += Physics::kineticEnergy( sp->getMass(), sp->getVelocity() );
      }
      return totalKE;
    };

    auto computeTotalPE = [&]()
    {
      float totalPE = 0.0f;
      for( size_t i = 0; i < physicalObjects.size() - 1; ++i )
      {
        for( size_t j = i + 1; j < physicalObjects.size(); ++j )
        {
          const Eigen::Vector2f r = physicalObjects[i]->getPosition() - physicalObjects[j]->getPosition();
          totalPE += Physics::gravitationalPotentialEnergy( physicalObjects[i]->getMass(), physicalObjects[j]->getMass(), r );
        }
      }
      return totalPE;
    };

    initialEnergyOld = computeTotalPE() + computeTotalKE();
    cout << "Initial energy = " << initialEnergyOld << endl;

    // Act
    sf::Clock clock1;
    float timef = float( t0 );
    float t1f = float( t1 );
    float dtf = float( dt );
    while( timef <= t1f )
    {
      motionManager.computeLinearDynamics( dtf, physicalObjects );
      timef += dtf;
    }

    cout << "Simulated time = " << timef << endl;
    cout << "Simulation time = " << clock1.getElapsedTime().asMilliseconds() << " ms" << endl;
    finalEnergyOld = computeTotalPE() + computeTotalKE();
    cout << "Final energy = " << finalEnergyOld << endl;
  }

  cout << endl << "-------------New integration" << endl;
  double initialEnergyNew, finalEnergyNew;
  {
    const size_t nbOfBodies = spaceObjectData.size();
    vector< double > masses( nbOfBodies );
    vector< double > y( nbOfBodies * 4 );
    for( size_t i = 0; i < nbOfBodies; ++i ) {
      const auto& sod = spaceObjectData[i];
      masses[i] = sod.mass;

      size_t j = i * 4;
      y[j] = sod.position( 0 );
      y[j + 1] = sod.position( 1 );
      y[j + 2] = sod.velocity( 0 );
      y[j + 3] = sod.velocity( 1 );
    }
    initialEnergyNew = Physics::computeTotalEnergy2D( y, masses );
    cout << "Initial energy = " << initialEnergyNew << endl;

    GravitationalNBody derivatives( masses );

    // Act
    double time = t0;
    sf::Clock clock2;
    OdeIntegrator< RK4< GravitationalNBody > > integrator( time,
                                                           y,
                                                           dt,
                                                           false,
                                                           derivatives );
    integrator.integrate( t1 - time );

    cout << "Simulated time = " << time << endl;
    cout << "Simulation time = " << clock2.getElapsedTime().asMilliseconds() << " ms" << endl;
    finalEnergyNew = Physics::computeTotalEnergy2D( y, masses );
    cout << "Final energy = " << finalEnergyNew << endl;
  }

  EXPECT_NEAR( initialEnergyNew, initialEnergyOld, abs( initialEnergyNew * 1e-5 ) );
  EXPECT_NEAR( finalEnergyNew, finalEnergyOld, abs( finalEnergyNew * 1e-5 ) );
}

TEST( Performance, Rotations )
{
  size_t nbOfTimes = size_t( 1e8 );
  const float angle = float( M_PI / 7.6f );

  sf::Clock clock;
  Eigen::Vector2f orientation( 1.0f, 0.0f );
  for( size_t i = 0; i < nbOfTimes; ++i )
  {
    MathUtil::rotate2D( orientation, angle );
  }

  cout << "Simulation time = " << clock.getElapsedTime().asMilliseconds() << " ms" << endl;
  cout << "Final orientation = " << Utilities::toString( orientation ) << endl;
  cout << "Normalized: " << orientation.norm() << endl;
}