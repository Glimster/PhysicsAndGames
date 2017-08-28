#pragma once
class PhysicalObject;

// TODO, hitta på bättre namn?
class MotionManager
{
public:
  enum class Integration { EulerForward, RK4 };
  enum class Technique { standard, functor };

  MotionManager( Integration integration, 
                 Technique technique );
  ~MotionManager();

  void computeLinearDynamics( const float dt,
                              const std::vector< PhysicalObject* >& physicalObjects );

  void computeRotationalDynamics( const float dt,
                                  PhysicalObject& physicalObject );
  
  // No forces
  void computeLinearKinematics( const float dt,
                                std::vector< PhysicalObject* > physicalObjects );

  std::string toString() 
  {
    switch( integration_ ) {
      case Integration::EulerForward:
        return "Euler forward";
      case Integration::RK4:
        return "Runge Kutta 4";
      default:
        return "Undefined";
    }
  }

private:
  Integration integration_;
  Technique technique_;
};

