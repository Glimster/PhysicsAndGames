#include "stdafx.h"
#include "Integrator.h"

#include "PhysicalObject.h"
#include "MathUtil.h"

using namespace std;

void Integrator::N2EulerStep( const vector< PhysicalObject* >& objects,
                            const float dt,
                            vector< Eigen::Vector2f >& drs,
                            vector< Eigen::Vector2f >& dvs )
{
  drs.resize( objects.size() );
  dvs.resize( objects.size() );
  for( size_t i = 0; i != objects.size(); ++i )
  {
    const Eigen::Vector2f& r1 = objects[i]->getPosition();

    Eigen::Vector2f FOverMTot( 0.0f, 0.0f );
    for( size_t j = 0; j != objects.size(); ++j ) {
      if( i == j )
        continue;

      FOverMTot += Physics::forceOfGravityOverM( objects[j]->getMass(), r1 - objects[j]->getPosition() );
    }

    // Thrusters
    Eigen::Vector2f FOverM2( 0.0f, 0.0f );
    objects[i]->computeLinearForceOverM( FOverM2 );

    FOverMTot += FOverM2;

    drs[i] = objects[i]->getVelocity() * dt;
    dvs[i] = FOverMTot * dt;
  }
}

void Integrator::N2RK4Step( const vector< PhysicalObject* >& objects,
                            const float dt,
                            vector< Eigen::Vector2f >& drs,
                            vector< Eigen::Vector2f >& dvs )
{
  const size_t nbOfObjects = objects.size();

  vector< Eigen::Vector2f > k1_r( nbOfObjects ), k1_v( nbOfObjects );
  vector< Eigen::Vector2f > k2_r( nbOfObjects ), k2_v( nbOfObjects );
  vector< Eigen::Vector2f > k3_r( nbOfObjects ), k3_v( nbOfObjects );
  vector< Eigen::Vector2f > k4_r( nbOfObjects ), k4_v( nbOfObjects );

  computeRK4Derivatives_( objects, 0.0f,
                          vector< Eigen::Vector2f >( nbOfObjects ), vector< Eigen::Vector2f >( nbOfObjects ),
                          k1_r, k1_v );

  computeRK4Derivatives_( objects, 0.5f * dt,
                          k1_r, k1_v,
                          k2_r, k2_v );

  computeRK4Derivatives_( objects, 0.5f * dt,
                          k2_r, k2_v,
                          k3_r, k3_v );

  computeRK4Derivatives_( objects, dt,
                          k3_r, k3_v,
                          k4_r, k4_v );

  drs.resize( objects.size() );
  dvs.resize( objects.size() );
  for( size_t i = 0; i < drs.size(); ++i ) {
    drs[i] = 1.0f / 6.0f * (k1_r[i] + 2.0f * (k2_r[i] + k3_r[i]) + k4_r[i]) * dt;
    dvs[i] = 1.0f / 6.0f * (k1_v[i] + 2.0f * (k2_v[i] + k3_v[i]) + k4_v[i]) * dt;
  }
}

void Integrator::N2EulerStepFunctors( const vector< PhysicalObject* >& objects,
                                      const float dt,
                                      vector< Eigen::Vector2f >& drs, 
                                      vector< Eigen::Vector2f >& dvs )
{
  const size_t nbOfObjects = objects.size();
  vector< StateLin > initialStates( nbOfObjects );
  for( size_t i = 0; i < nbOfObjects; ++i ) {
    initialStates[i].mass = objects[i]->getMass();
    initialStates[i].position = objects[i]->getPosition();
    initialStates[i].velocity = objects[i]->getVelocity();
    objects[i]->computeLinearForceOverM( initialStates[i].linearForceOverM );
  }

  drs.resize( nbOfObjects );
  dvs.resize( nbOfObjects );
  for( size_t i = 0; i < nbOfObjects; ++i )
  {
    // dv/dt = a = F / m
    LinearForceOverM fOverM( initialStates, i );
    dvs[i] = MathUtil::forwardEulerStep( fOverM, dt );

    // dr/dt = v
    drs[i] = MathUtil::forwardEulerStep( [&]() { return initialStates[i].velocity; }, dt );
  }
}

void Integrator::N2RK4StepFunctors( const vector< PhysicalObject* >& objects,
                                    const float dt,
                                    vector< Eigen::Vector2f >& drs,
                                    vector< Eigen::Vector2f >& dvs )
{
  const size_t nbOfObjects = objects.size();
  vector< StateLin > initialStates( nbOfObjects );
  for( size_t i = 0; i < nbOfObjects; ++i ) {
    initialStates[i].mass = objects[i]->getMass();
    initialStates[i].position = objects[i]->getPosition();
    initialStates[i].velocity = objects[i]->getVelocity();
    objects[i]->computeLinearForceOverM( initialStates[i].linearForceOverM );
  }

  vector< Eigen::Vector2f > k1_r( nbOfObjects ), k1_v( nbOfObjects );
  vector< Eigen::Vector2f > k2_r( nbOfObjects ), k2_v( nbOfObjects );
  vector< Eigen::Vector2f > k3_r( nbOfObjects ), k3_v( nbOfObjects );
  vector< Eigen::Vector2f > k4_r( nbOfObjects ), k4_v( nbOfObjects );

  computeRK4Derivatives_( initialStates, 0.0f, 
                          vector< Eigen::Vector2f >( nbOfObjects, Eigen::Vector2f( 0.0f, 0.0f ) ), vector< Eigen::Vector2f >( nbOfObjects, Eigen::Vector2f( 0.0f, 0.0f ) ),
                          k1_r, k1_v );

  computeRK4Derivatives_( initialStates, 0.5f * dt, 
                          k1_r, k1_v,
                          k2_r, k2_v );

  computeRK4Derivatives_( initialStates, 0.5f * dt, 
                          k2_r, k2_v,
                          k3_r, k3_v );

  computeRK4Derivatives_( initialStates, dt, 
                          k3_r, k3_v,
                          k4_r, k4_v );

  drs.resize( nbOfObjects );
  dvs.resize( nbOfObjects );
  for( size_t i = 0; i < drs.size(); ++i ) {
    drs[i] = MathUtil::RK4Step( k1_r[i], k2_r[i], k3_r[i], k4_r[i], dt );
    dvs[i] = MathUtil::RK4Step( k1_v[i], k2_v[i], k3_v[i], k4_v[i], dt );
  }
}

void Integrator::N2RotEulerStepFunctor( const PhysicalObject& object,
                                        const float dt,
                                        float& dL,
                                        float& dTheta )
{
  StateRot initialState = { object.getAngularMomentum(), object.getMomentOfInertia(), object.computeTorque() };

  dL = MathUtil::forwardEulerStep( [&]() { return initialState.torque; }, dt );

  AngularVelocity omega( initialState );
  dTheta = MathUtil::forwardEulerStep( omega, dt );
}

void Integrator::N2RotRK4StepFunctor( const PhysicalObject& object,
                                      const float dt,
                                      float& dL,
                                      float& dTheta )
{
  StateRot initialState = { object.getAngularMomentum(), object.getMomentOfInertia(), object.computeTorque() };

  float k1_L = 0.0f, k2_L = 0.0f, k3_L = 0.0f, k4_L = 0.0f;
  float k1_theta = 0.0f, k2_theta = 0.0f, k3_theta = 0.0f, k4_theta = 0.0f;

  computeRK4Derivatives_( initialState, 0.0f, 
                          0.0f, 0.0f,
                          k1_L, k1_theta );

  computeRK4Derivatives_( initialState, 0.5f * dt, 
                          k1_L, k1_theta,
                          k2_L, k2_theta );

  computeRK4Derivatives_( initialState, 0.5f * dt, 
                          k2_L, k2_theta,
                          k3_L, k3_theta );

  computeRK4Derivatives_( initialState, dt, 
                          k3_L, k3_theta,
                          k4_L, k4_theta );

  dL = MathUtil::RK4Step( k1_L, k2_L, k3_L, k4_L, dt );
  dTheta = MathUtil::RK4Step( k1_theta, k2_theta, k3_theta, k4_theta, dt );
}

// ---------------------------------------------------------------------------
void Integrator::computeRK4Derivatives_( const vector< PhysicalObject* >& objects,
                                         const float dt,
                                         const vector< Eigen::Vector2f >& drdtIn,
                                         const vector< Eigen::Vector2f >& dvdtIn,
                                         vector< Eigen::Vector2f >& drdtOut,
                                         vector< Eigen::Vector2f >& dvdtOut )
{
  vector< Eigen::Vector2f > rs( objects.size() );
  vector< Eigen::Vector2f > vs( objects.size() );
  for( size_t i = 0; i < objects.size(); ++i ) {
    rs[i] = objects[i]->getPosition() + drdtIn[i] * dt;
    vs[i] = objects[i]->getVelocity() + dvdtIn[i] * dt;
  }

  for( size_t i = 0; i < objects.size(); ++i ) {
    Eigen::Vector2f FOverMTot( 0.0f, 0.0f );
    for( size_t j = 0; j != objects.size(); ++j ) {
      if( i == j )
        continue;

      FOverMTot += Physics::forceOfGravityOverM( objects[j]->getMass(), rs[i] - rs[j] );
    }

    // Thrusters
    Eigen::Vector2f FOverM2;
    objects[i]->computeLinearForceOverM( FOverM2 );

    FOverMTot += FOverM2;

    drdtOut[i] = vs[i];
    dvdtOut[i] = FOverMTot;
  }
}

void Integrator::computeRK4Derivatives_( const StateRot& initialState,
                                         const float dt,
                                         const float dLdtIn, 
                                         const float dThetadtIn,
                                         float& dLdtOut, 
                                         float& dThetadtOut )
{
  StateRot intermediateState( initialState ); // Make sure all values are copied
  intermediateState.angularMomentum += dLdtIn * dt;

  dLdtOut = intermediateState.torque;

  AngularVelocity w( intermediateState );
  dThetadtOut = w();
}

void Integrator::computeRK4Derivatives_( const vector< StateLin >& initialStates, 
                                         const float dt,
                                         const vector< Eigen::Vector2f >& drdtIn, 
                                         const vector< Eigen::Vector2f >& dvdtIn,
                                         vector< Eigen::Vector2f >& drdtOut, 
                                         vector< Eigen::Vector2f >& dvdtOut )
{
  size_t nbOfStates = initialStates.size();
  vector< StateLin > intermediateStates( initialStates ); // Make sure all values are copied
  for( size_t i = 0; i < nbOfStates; ++i ) {
    intermediateStates[i].position += drdtIn[i] * dt;
    intermediateStates[i].velocity += dvdtIn[i] * dt;
  }

  for( size_t i = 0; i < nbOfStates; ++i ) {
    LinearForceOverM linearForceOverM( intermediateStates, i );
    Eigen::Vector2f fOverM = linearForceOverM();

    drdtOut[i] = intermediateStates[i].velocity;
    dvdtOut[i] = fOverM;
  }
}