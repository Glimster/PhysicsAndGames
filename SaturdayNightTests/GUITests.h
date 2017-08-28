#pragma once
#include <array>
#include "ResourceHolder.h"
#include "ResourceIdentifiers.h"
#include "CoordinateSystemHandler.h"
#include "PhysicalData.h"
#include "SceneNode.h"

class PhysicalObject;
class PhysicalData;
class Planet;
class Ball;

class GUITests
{
public:
  enum class Setup
  {
    Space,
    RollingBalls
  };

public:

	GUITests();
  void init( Setup setup );
  ~GUITests();

	void runSpaceSimulation();
  void runCollisionSimulation();

private:

	enum Layer
	{
		Background,
		ActionLayer,
		LayerCount
	};

  void loadResources_();
	void buildScene_( const std::vector< PhysicalData::PlanetData >& data );
	void render_();

  void renderDebugPrintouts_();
  void setupDebug_();
  void updateFPSStats_( const sf::Time& elapsedTime );

private:
	sf::RenderWindow mainWindow_;
	TextureHolder textureHolder_;
  FontHolder fontHolder_;
	sf::View view_;
	CoordinateSystemHandler csHandler_;

	std::array< SceneNode*, LayerCount > sceneLayers_;
	SceneNode sceneGraph_;

	std::vector< PhysicalObject* > physicalObjects_;
  std::vector< Ball* > balls_;

  sf::Time dt_;
  float timeScale_;

  // Debug
  sf::Text debugText_;
  sf::Text fpsText_;
  sf::Time fpsUpdateTime_;
  std::size_t fpsNumFrames_;
};