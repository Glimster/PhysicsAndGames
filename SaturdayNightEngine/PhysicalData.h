#pragma once

// TODO, ska man göra samma lösning här som för GameDataTables?
class PhysicalData
{
public:

  struct PlanetData
  {
    std::string name;
    float mass;
    float radius;
    Eigen::Vector2f velocity;
    Eigen::Vector2f position;
  };

  struct BallData
  {
    float radius;
    Eigen::Vector2f velocity;
    Eigen::Vector2f position;
  };

  static std::vector< PlanetData > setupHeavySunLightPlanet();
  static std::vector< PlanetData > setupTwoBodySystem();
  static std::vector< PlanetData > setupThreeBodySystem();
  static std::vector< PlanetData > setupPlanetarySystem();
  static std::vector< PlanetData > setupRealisticSolarSystem();
  static std::vector< PlanetData > setupAlotOfPlanets();

  static std::vector< BallData > setupMovingAndStationaryBall();
  static std::vector< BallData > setupTwoBalls();
  static std::vector< BallData > setupALotOfBalls();
};

