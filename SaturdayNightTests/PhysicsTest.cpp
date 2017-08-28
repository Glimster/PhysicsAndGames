#include "stdafx.h"
#include "PhysicsTest.h"

//#include <gtest/gtest.h>

#include "PhysicalObject.h"
#include "Ball.h"
#include "Physics.h"
#include "PhysicalData.h"
#include "PhysicalConstants.h"
#include "MathUtil.h"
//#include "MotionManager.h"
#include "Utilities.h"
#include "CollisionControl.h"

using namespace std;

//PhysicsTest::PhysicsTest():
//{}

INSTANTIATE_TEST_CASE_P( IntegrationTechniques, 
                         PhysicsTest,
                         ::testing::Values( MotionManager::Integration::EulerForward, MotionManager::Integration::RK4 ) );

// Tests that a (light) planet orbiting around a (heavy) sun follows Kepler's
// laws of planetary motion
// https://en.wikipedia.org/wiki/Kepler%27s_laws_of_planetary_motion
//
// Fails for forward Euler
// Succeeds for RK4
TEST_P( PhysicsTest, KeplersLaws )
{
  // Setup
  auto integrationMethod = GetParam();
  MotionManager motionManager( integrationMethod,
                               MotionManager::Technique::functor );
  cout << "Integration technique: " << motionManager.toString() << endl;

  vector< unique_ptr< PhysicalObject > > celestialBodies;
  vector< PhysicalObject* > physicalObjects;

  vector< PhysicalData::PlanetData > spaceObjectData;
  PhysicalData::setupHeavySunLightPlanet( spaceObjectData );
  const auto& sunData = spaceObjectData[0];
  for( auto& sod : spaceObjectData )
  {
    const PhysicalData::PlanetData& data = sod;
    unique_ptr< PhysicalObject > object( new PhysicalObject( data.mass ) );
    object->setPosition( data.position );
    object->setVelocity( data.velocity );

    physicalObjects.push_back( object.get() );

    celestialBodies.push_back( std::move( object ) );
  }
  const auto& sun = *celestialBodies[0];
  const auto& planet = *celestialBodies[1];

  // Act
  sf::Clock clock;
  size_t updates = 0;
  float totalTime = 0.0f;
  float distanceTravelled = 0.0f;
  float totalArea = 0.0f;
  float totalAngle = 0.0f;
  
  vector< float > sweptAreas;

  auto startingPosition = planet.getPosition();
  auto previousPosition = startingPosition;
  auto minPosition = startingPosition;
  auto maxPosition = startingPosition;
  // Distances to sun
  float dPericenter = ( startingPosition - sun.getPosition() ).norm();
  float dApocenter = ( startingPosition - sun.getPosition() ).norm();
  float vPericenter = planet.getVelocity().norm();
  float vApocenter = planet.getVelocity().norm();
  auto pericenter = startingPosition;
  auto apocenter = startingPosition;

  float sweptArea = 0.0f;
  float dt = 1.0f / 1000.0f;
  size_t segmentAreaStepInterval = size_t( 0.5f / dt );
  while( totalAngle <= 2.0f * float( M_PI ) )
  {
    updates++;

    motionManager.computeLinearDynamics( dt, physicalObjects );

    auto r_prev = previousPosition - sun.getPosition();
    auto r = planet.getPosition() - sun.getPosition();
    auto d = ( planet.getPosition() - previousPosition ).norm();
    
    totalTime += dt;
    distanceTravelled += d;
    totalArea += MathUtil::cross2D( r, r_prev ) / 2.0f; // Area of triangle
    sweptArea += MathUtil::cross2D( r, r_prev ) / 2.0f; // Area of triangle
    if( updates % segmentAreaStepInterval == 0 )
    {
      sweptAreas.push_back( sweptArea );
      sweptArea = 0.0f;
    }
    { // Floating point precision not enough, use double
      Eigen::Vector2d rd = r.cast< double >();
      Eigen::Vector2d r_prevd = r_prev.cast< double >();
      double dAngle = MathUtil::angleBetween( rd, r_prevd );
      totalAngle += float( dAngle );
    }

    minPosition(0) = min( planet.getPosition()(0), minPosition(0) );
    minPosition(1) = min( planet.getPosition()(1), minPosition(1) );
    maxPosition(0) = max( planet.getPosition()(0), maxPosition(0) );
    maxPosition(1) = max( planet.getPosition()(1), maxPosition(1) );

    if( r.norm() < dPericenter )
    {
      pericenter = planet.getPosition();
      vPericenter = planet.getVelocity().norm();
      dPericenter = r.norm();
    }
    if( r.norm() > dApocenter )
    {
      apocenter = planet.getPosition();
      vApocenter = planet.getVelocity().norm();
      dApocenter = r.norm();
    }

    previousPosition = planet.getPosition();
  }

  cout << "Simulation time = " << clock.getElapsedTime().asMilliseconds() << " ms" << endl;
  cout << "Number of steps = " << updates << endl;
  cout << "Final position = " << Utilities::toString( planet.getPosition() ) << endl;

  if( integrationMethod == MotionManager::Integration::EulerForward )
  {
    cout << "Expected to fail!" << endl;
    EXPECT_GT( (startingPosition - planet.getPosition()).norm(), 1e-3 );
    return;
  }

  EXPECT_NEAR( ( startingPosition - planet.getPosition() ).norm(), 0.0, 1e-3 );

  // -------------------------------------------
  cout << endl << "Kepler 1: The orbit of a planet is an ellipse with the Sun at one of the two foci." << endl;
  auto centerPosition = ( maxPosition + minPosition ) / 2.0f;
  auto axes = maxPosition - minPosition;
  float semiMajorAxis = axes.maxCoeff() / 2.0f;
  float semiMinorAxis = axes.minCoeff() / 2.0f;
  cout << "Semi major axis = " << semiMajorAxis << endl;
  cout << "Semi minor axis = " << semiMinorAxis << endl;

  cout << "pericenter =  " << Utilities::toString( pericenter ) << endl;
  cout << "apocenter =  " << Utilities::toString( apocenter ) << endl;

  // According to https://en.wikipedia.org/wiki/Apsis
  float eccentricity = ( dApocenter - dPericenter ) / ( dApocenter + dPericenter );
  float standardGravitationalParameter = Phys::PhysicalConstants::kSquared * ( sun.getMass() );
  cout << "Distance at pericenter = " << dPericenter << endl;
  cout << "Speed at pericenter = " << vPericenter << endl;
  cout << "Distance at apocenter = " << dApocenter << endl;
  cout << "Speed at apocenter = " << vApocenter << endl;
  EXPECT_NEAR( dPericenter, ( 1 - eccentricity ) * semiMajorAxis, 1e-3 );
  EXPECT_NEAR( vPericenter, sqrt( ( 1 + eccentricity ) * standardGravitationalParameter / ( ( 1 - eccentricity ) * semiMajorAxis ) ), 1e-3 );
  EXPECT_NEAR( dApocenter, ( 1 + eccentricity ) * semiMajorAxis, 1e-3 );
  EXPECT_NEAR( vApocenter, sqrt( ( 1 - eccentricity ) * standardGravitationalParameter / ( ( 1 + eccentricity ) * semiMajorAxis ) ), 1e-3 );

  float f = MathUtil::fociCenterDistanceEllipse( semiMajorAxis, semiMinorAxis );
  Eigen::Vector2f F1( centerPosition( 0 ) - f, centerPosition( 1 ) );
  Eigen::Vector2f F2( centerPosition( 0 ) + f, centerPosition( 1 ) );
  cout << "Focal point #1 =  " << Utilities::toString( F1 ) << endl;
  cout << "Focal point #2 (sun) =  " << Utilities::toString( F2 ) << endl;
  EXPECT_NEAR( ( sun.getPosition() - sunData.position ).norm(), 0.0, 1e-4 ); // Sun "infinitely" heavy
  EXPECT_NEAR( ( F2 - sun.getPosition() ).norm(), 0.0, 1e-4 ); // Sun at focal point two
  
  cout << "Distance travelled = " << distanceTravelled << endl;
  EXPECT_NEAR( distanceTravelled, MathUtil::circumferenceEllipse( semiMajorAxis, semiMinorAxis ), 1e-3 );

  cout << "Total swept area = " << totalArea << endl;
  EXPECT_NEAR( totalArea, MathUtil::areaEllipse( semiMajorAxis, semiMinorAxis ), 1e-3 );

  cout << "Total angle = " << totalAngle << endl;
  EXPECT_NEAR( totalAngle, 2 * M_PI, 1e-3 );

  // -------------------------------------------
  cout << endl << "Kepler 2: A line joining a planet and the Sun sweeps out equal areas during equal intervals of time." << endl;
  cout << "Swept area in time interval = " << sweptAreas[0] << endl;
  for( size_t i = 1; i < sweptAreas.size(); ++i )
  {
    EXPECT_NEAR( sweptAreas[i], sweptAreas[0], 1e-4 );
  }

  // -------------------------------------------
  cout << endl << "Kepler 3: The square of the orbital period of a planet is directly proportional to the cube of the semi-major axis of its orbit." << endl;
  const float constant = 4.0f * float( M_PI * M_PI ) / standardGravitationalParameter;
  float T = sqrt( pow( semiMajorAxis, 3 ) * constant );
  cout << "Orbital period = " << totalTime << endl;
  EXPECT_NEAR( totalTime, T, 1e-2 );
}

// Tests that the total energy remains constant for 1000 turns
// Fails for forward Euler
// Succeeds for RK4
TEST_P( PhysicsTest, StabililtyWRTEnergy )
{
  cout << endl << "--------------------------------" << endl;
  cout << "Test stability" << endl;
  // Setup
  auto integrationMethod = GetParam();
  MotionManager motionManager( integrationMethod,
                               MotionManager::Technique::functor );
  cout << "Integration technique: " << motionManager.toString() << endl;

  vector< unique_ptr< PhysicalObject > > celestialObjects;
  vector< PhysicalObject* > physicalObjects;

  vector< PhysicalData::PlanetData > spaceObjectData;
  PhysicalData::setupHeavySunLightPlanet( spaceObjectData );
  for( auto& sod : spaceObjectData )
  {
    const PhysicalData::PlanetData& data = sod;
    unique_ptr< PhysicalObject > object( new PhysicalObject( data.mass ) );
    object->setPosition( data.position );
    object->setVelocity( data.velocity );

    physicalObjects.push_back( object.get() );
    celestialObjects.push_back( std::move( object ) );
  }
  const auto& sun = *celestialObjects[0];
  const auto& planet = *celestialObjects[1];

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
        totalPE = Physics::gravitationalPotentialEnergy( physicalObjects[i]->getMass(), physicalObjects[j]->getMass(), r );
      }
    }
    return totalPE;
  };

  cout << "Initial position =  " << Utilities::toString( planet.getPosition() ) << endl;
  cout << "Initial distance = " << ( planet.getPosition() - sun.getPosition() ).norm() << endl;
  const float initialEnergy = computeTotalPE() + computeTotalKE();
  cout << "Initial total energy = " << initialEnergy << endl;

  // Act
  sf::Clock clock;
  float totalAngle = 0.0f;
  auto startingPosition = planet.getPosition();
  auto previousPosition = startingPosition;
  float dt = 1.0f / 1000.0f;

  int stable = true;
  int numberOfTurns = 1000;
  while( int( totalAngle / ( 2.0f * float( M_PI ) ) ) < numberOfTurns )
  {
    motionManager.computeLinearDynamics( dt, physicalObjects );

    auto r_prev = previousPosition - sun.getPosition();
    auto r = planet.getPosition() - sun.getPosition();

    { // Floating point computation not sufficient, use double precision
      Eigen::Vector2d rd = r.cast< double >();
      Eigen::Vector2d r_prevd = r_prev.cast< double >();
      double dAngle = MathUtil::angleBetween( rd, r_prevd );
      totalAngle += float( dAngle );
    }

    const float totalEnergy = computeTotalPE() + computeTotalKE();
    if( fabs( totalEnergy - initialEnergy ) > fabs( initialEnergy * 0.001f ) )
    {
      stable = false;
      cout << "Total energy of system not constant, terminating!!!" << endl;
      break;
    }
    previousPosition = planet.getPosition();
  }

  if( !stable )
  {
    EXPECT_EQ( MotionManager::Integration::EulerForward, integrationMethod );
    return;
  }

  cout << "---------------------------------------" << endl;
  cout << "Simulation time = " << clock.getElapsedTime().asMilliseconds() << " ms" << endl;
  cout << "Final position =  " << Utilities::toString( planet.getPosition() ) << endl;

  cout << "Total angle = " << totalAngle << endl;
  cout << "Number of turns = " << totalAngle / ( 2.0f * float( M_PI ) ) << endl;
  EXPECT_NEAR( totalAngle, numberOfTurns * 2 * M_PI, 1e-3 );

  float finalEnergy = computeTotalPE() + computeTotalKE();
  cout << "Final total energy = " << finalEnergy << endl;
  EXPECT_NEAR( initialEnergy, finalEnergy, 1e-3 );
}

TEST( Performance, aLotOfPlanets )
{
  cout << endl << "--------------------------------" << endl;
  cout << "Test performance for a lot of planets" << endl;
  // Setup
  MotionManager motionManager( MotionManager::Integration::RK4,
                               MotionManager::Technique::standard );

  vector< PhysicalObject* > physicalObjects;

  vector< PhysicalData::PlanetData > spaceObjectData;
  PhysicalData::setupAlotOfPlanets( spaceObjectData );
  for( auto& sod : spaceObjectData )
  {
    const PhysicalData::PlanetData& data = sod;
    unique_ptr< PhysicalObject > object( new PhysicalObject( data.mass ) );
    object->setPosition( data.position );
    object->setVelocity( data.velocity );

    physicalObjects.push_back( object.get() );
  }

  // Act
  sf::Clock clock;
  float dt = 1.0f / 100.0f;
  size_t i = 0;
  while( i < 1000 )
  {
    motionManager.computeLinearDynamics( dt, physicalObjects );
    ++i;
  }

  cout << "---------------------------------------" << endl;
  cout << "Simulation time = " << clock.getElapsedTime().asMilliseconds() << " ms" << endl;
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

TEST( ElasticCollisions, ConservationOfEnergyAndMomentum )
{
  // Setup
  CoordinateSystemHandler csHandler; // TODO, should not be necessary, bad design
  MotionManager motionManager( MotionManager::Integration::RK4,
                               MotionManager::Technique::functor );

  vector< string > testNames( { "Moving and stationary ball", 
                                "Two colliding balls", 
                                "A lot of balls" } ); // TODO, failar just nu test för bevarande av energi och rörelsemängd
  vector< vector< PhysicalData::BallData > > ballDataV( 3 );
  PhysicalData::setupMovingAndStationaryBall( ballDataV[0] );
  PhysicalData::setupTwoBalls( ballDataV[1] );
  PhysicalData::setupALotOfBalls( ballDataV[2] );
  for( size_t i = 0; i < 3; ++i )
  {
    cout << endl << testNames[i] << endl;
    const auto& ballData = ballDataV[i];
    // TODO, stoppa i gemensamt ställe för GUITests också?
    vector< unique_ptr< Ball > > objects;
    vector< PhysicalObject* > physicalObjects;
    vector< Ball* > balls;
    for( auto data : ballData )
    {
      float mass = float( M_PI ) * data.radius * data.radius;
      std::unique_ptr< Ball > ball( new Ball( data.radius, mass, csHandler ) );
      ball->setPosition( data.position );
      ball->setVelocity( data.velocity );
      physicalObjects.push_back( ball.get() );
      balls.push_back( ball.get() );
      objects.push_back( std::move( ball ) );
    }

    const auto& b0 = *objects[0];
    const auto& b1 = *objects[1];

    auto v0 = b0.getVelocity();
    auto v1 = b1.getVelocity();
    cout << "Ball 0 velocity =  " << Utilities::toString( v0 ) << endl;
    cout << "Ball 1 velocity =  " << Utilities::toString( v1 ) << endl;

    auto computeTotalKE = [&]()
    {
      float totalKE = 0.0f;
      for( const auto& sp : physicalObjects )
      {
        totalKE += Physics::kineticEnergy( sp->getMass(), sp->getVelocity() );
      }
      return totalKE;
    };

    auto computeTotalMomentum = [&]()
    {
      Eigen::Vector2f totalMomentum( 0.0f, 0.0f );
      for( const auto& sp : physicalObjects )
      {
        totalMomentum += Physics::momentum( sp->getMass(), sp->getVelocity() );
      }
      return totalMomentum;
    };

    auto initialKE = computeTotalKE();
    auto initialMomentum = computeTotalMomentum();
    cout << "Initial total energy = " << initialKE << endl;
    cout << "Initial total momentum = " << Utilities::toString( initialMomentum ) << endl;

    // Act
    cout << "Simulate collision" << endl;
    float totalSimulationTime = 0.0f;
    const float dt = 1.0f / 1000.0f;
    while( totalSimulationTime < 4.0f )
    {
      motionManager.computeLinearKinematics( dt, physicalObjects );
      CollisionControl::handleCollisions( balls );
      totalSimulationTime += dt;
    }

    // Verify
    auto v0p = b0.getVelocity();
    auto v1p = b1.getVelocity();

    cout << "Ball 0 velocity =  " << Utilities::toString( v0p ) << endl;
    cout << "Ball 1 velocity =  " << Utilities::toString( v1p ) << endl;
    if( i == 0 )
    {
      cout << "Balls of equal masses should exchange velocities" << endl;
      EXPECT_NEAR( v0p( 0 ), v1( 0 ), 1e-6 );
      EXPECT_NEAR( v0p( 1 ), v1( 1 ), 1e-6 );
      EXPECT_NEAR( v1p( 0 ), v0( 0 ), 1e-6 );
      EXPECT_NEAR( v1p( 1 ), v0( 1 ), 1e-6 );
    }
    auto finalKE = computeTotalKE();
    auto finalMomentum = computeTotalMomentum();
    cout << "Final total energy = " << finalKE << endl;
    cout << "Final total momentum = " << Utilities::toString( finalMomentum ) << endl;
    EXPECT_NEAR( finalKE, initialKE, initialKE * 1e-6 );
    EXPECT_NEAR( finalMomentum( 0 ), initialMomentum( 0 ), initialMomentum( 0 ) * 1e-6 );
    EXPECT_NEAR( finalMomentum( 1 ), initialMomentum( 1 ), initialMomentum( 0 ) * 1e-6 );
  }
}