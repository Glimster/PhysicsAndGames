#include "stdafx.h"
#include "CollissionTests.h"

#include "CoordinateSystemHandler.h"
#include "PhysicalData.h"
#include "Ball.h"
#include "Utilities.h"
#include "Physics.h"
#include "CollisionControl.h"

using namespace std;

//CollissionTests::CollissionTests():
//{}

//INSTANTIATE_TEST_CASE_P( IntegrationTechniques, 
//                         CollissionTests,
//                         ::testing::Values( MotionManager::Integration::EulerForward, 
//                                            MotionManager::Integration::RK4 ) );

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
  ballDataV[0] = PhysicalData::setupMovingAndStationaryBall();
  ballDataV[1] = PhysicalData::setupTwoBalls();
  ballDataV[2] = PhysicalData::setupALotOfBalls();
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