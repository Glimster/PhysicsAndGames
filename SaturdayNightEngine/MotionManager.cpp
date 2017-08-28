#include "stdafx.h"
#include "MotionManager.h"

#include "Integrator.h"
#include "MathUtil.h"
#include "PhysicalObject.h"

using namespace std;

MotionManager::MotionManager( Integration integration,
                              Technique technique ):
integration_( integration ),
technique_( technique )
{}

MotionManager::~MotionManager() {}

void MotionManager::computeLinearDynamics( const float dt,
                                           const vector< PhysicalObject* >& physicalObjects )
{
  vector< Eigen::Vector2f > drs, dvs;
  switch( technique_ )
  {
  case( Technique::standard ) :
  {
    if( integration_ == Integration::EulerForward ) {
      Integrator::N2EulerStep( physicalObjects, dt, drs, dvs );
    }
    else {
      assert( integration_ == Integration::RK4 );
      Integrator::N2RK4Step( physicalObjects, dt, drs, dvs );
    }
    break;
  }
  case( Technique::functor ) :
  {
    if( integration_ == Integration::EulerForward ) {
      Integrator::N2EulerStepFunctors( physicalObjects, dt, drs, dvs );
    }
    else {
      assert( integration_ == Integration::RK4 );
      Integrator::N2RK4StepFunctors( physicalObjects, dt, drs, dvs );
    }
    break;
  }
  default:
    assert( false );
    break;
  }

  for( size_t i = 0; i != physicalObjects.size(); ++i )
  {
    physicalObjects[i]->updateState( drs[i], dvs[i] );
  }
}

void MotionManager::computeRotationalDynamics( const float dt,
                                               PhysicalObject& physicalObject )
{
  float dL, dTheta;
  switch( technique_ )
  {
  case( Technique::standard ):
  {
    assert( integration_ != MotionManager::Integration::RK4 );
    const float torque = physicalObject.computeTorque();
    dL = torque * dt;

    const float angularVelocity = physicalObject.getAngularMomentum() / physicalObject.getMomentOfInertia();
    dTheta = angularVelocity * dt;
    break;
  }
  case( Technique::functor ):
  {
    if( integration_ == MotionManager::Integration::EulerForward ) {
      Integrator::N2RotEulerStepFunctor( physicalObject, dt, dL, dTheta );
    }
    else {
      assert( integration_ == MotionManager::Integration::RK4 );
      Integrator::N2RotRK4StepFunctor( physicalObject, dt, dL, dTheta );
    }
    break;
  }
  default:
    assert( false );
    break;
  }

  physicalObject.updateAngularMomentum( dL );
  physicalObject.rotate( dTheta );
}

void MotionManager::computeLinearKinematics( const float dt,
                                             vector< PhysicalObject* > physicalObjects )
{
  for( size_t i = 0; i != physicalObjects.size(); ++i )
  {
    physicalObjects[i]->updateState( physicalObjects[i]->getVelocity() * dt, 
                                     Eigen::Vector2f( 0.0f, 0.0f ) );
  }
}