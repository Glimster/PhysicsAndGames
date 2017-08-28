#pragma once

#include "Physics.h"
class PhysicalObject;

struct StateLin
{
  Eigen::Vector2f position;
  Eigen::Vector2f velocity;
  float mass;
  Eigen::Vector2f linearForceOverM;
};

struct StateRot
{
  float angularMomentum;
  float momentOfInertia;
  float torque;
};

class Integrator
{
public:

  // Newton's second law
  //
  // Second order ODE
  static void N2EulerStep( const std::vector< PhysicalObject* >& objects,
                           const float dt,
                           std::vector< Eigen::Vector2f >& drs,
                           std::vector< Eigen::Vector2f >& dvs );
  static void N2RK4Step( const std::vector< PhysicalObject* >& objects,
                         const float dt,
                         std::vector< Eigen::Vector2f >& drs,
                         std::vector< Eigen::Vector2f >& dvs );

  // Använder funktorer med states
  static void N2EulerStepFunctors( const std::vector< PhysicalObject* >& objects, 
                                   const float dt,
                                   std::vector< Eigen::Vector2f >& drs, 
                                   std::vector< Eigen::Vector2f >& dvs );

  static void N2RK4StepFunctors( const std::vector< PhysicalObject* >& objects,
                                 const float dt,
                                 std::vector< Eigen::Vector2f >& drs,
                                 std::vector< Eigen::Vector2f >& dvs );
  
  // Newton's second law for rotation
  //
  // Second order ODE
  static void N2RotEulerStepFunctor( const PhysicalObject& object, 
                                     const float dt,
                                     float& dL,
                                     float& dTheta );

  static void N2RotRK4StepFunctor( const PhysicalObject& object,
                                   const float dt,
                                   float& dL,
                                   float& dTheta );
private:

  // Computes dr/dt (=v) and dv/dt (=a) at time t + dt
  // Helper class for the standard fourth order Runge-Kutta integrator
  static void computeRK4Derivatives_( const std::vector< PhysicalObject* >& objects,
                                      const float dt,
                                      const std::vector< Eigen::Vector2f >& drdtIn,
                                      const std::vector< Eigen::Vector2f >& dvdtIn,
                                      std::vector< Eigen::Vector2f >& drdtOut,
                                      std::vector< Eigen::Vector2f >& dvdtOut );

  // Computes dr/dt (=v) and dv/dt (=a) at time t + dt
  // Helper class for the fourth order Runge-Kutta integrator
  static void computeRK4Derivatives_( const StateRot& initialState,
                                      const float dt,
                                      const float dLdtIn, 
                                      const float dThetadtIn,
                                      float& dLdtOut, 
                                      float& dThetadtOut );

  // Computes dr/dt (=v) and dv/dt (=a) at time t + dt
  // Helper class for the fourth order Runge-Kutta integrator
  static void computeRK4Derivatives_( const std::vector< StateLin >& initialStates,
                                      const float dt,
                                      const std::vector< Eigen::Vector2f >& drdtIn, 
                                      const std::vector< Eigen::Vector2f >& dvdtIn,
                                      std::vector< Eigen::Vector2f >& drdtOut, 
                                      std::vector< Eigen::Vector2f >& dvdtOut );
};

//--- Functors

// Computes all linear forces on bodies_(index_)
//
// - Gravity from N - 1 bodies
// - Forward thrusters
class LinearForceOverM
{
public:
  LinearForceOverM( const std::vector< StateLin >& bodies, size_t index ):
    bodies_( bodies ),
    index_( index )
  {}

public:
  Eigen::Vector2f operator() ()
  {
    Eigen::Vector2f fOverM( 0.0f, 0.0f );
    for( size_t j = 0; j != bodies_.size(); ++j ) {
      if ( index_ == j )
        continue;

      fOverM += Physics::forceOfGravityOverM( bodies_[j].mass,
                                              bodies_[index_].position - bodies_[j].position );
    }

    // Thrusters
    fOverM += bodies_[index_].linearForceOverM;
    
    return fOverM;
  }

private:
  const size_t index_; // Index of body in question
  const std::vector< StateLin >& bodies_;
};

class AngularVelocity
{
public:
  AngularVelocity( const StateRot& state ): state_( state ) {}
  float operator() ()
  {
    return state_.angularMomentum / state_.momentOfInertia;
  }
private:
  const StateRot& state_;
};