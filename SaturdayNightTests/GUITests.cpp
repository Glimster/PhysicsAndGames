#include "stdafx.h"
#include "GUITests.h"

#include "MotionManager.h"
#include "Planet.h"
#include "Ball.h"
#include "CollisionControl.h"

#include "CommandQueue.h"

GUITests::GUITests() :
  mainWindow_(),
  textureHolder_(),
  view_(),
  csHandler_(),
  physicalObjects_(),
  dt_( sf::seconds( 1.0f / 60.0f ) ),
  timeScale_( 1.0f )
{
  mainWindow_.create( sf::VideoMode::getDesktopMode(), "GUI tests" );
  view_ = mainWindow_.getDefaultView();

  loadResources_();
  setupDebug_();
}

void GUITests::init( Setup setup )
{
  switch( setup ) {
  case GUITests::Setup::Space:
  {
    dt_ = sf::seconds( 1.0f / 50.0f );
    timeScale_ = 5.0f;
    vector< PhysicalData::PlanetData > planetarySystemData;
    //PhysicalData::setupHeavySunLightPlanet( planetarySystemData );
    //PhysicalData::setupTwoBodySystem( planetarySystemData );
    //PhysicalData::setupThreeBodySystem( planetarySystemData );
    //PhysicalData::setupAlotOfPlanets( planetarySystemData );
    PhysicalData::setupRealisticSolarSystem( planetarySystemData );

    auto planetData = max_element( planetarySystemData.begin(), planetarySystemData.end(), []
    ( const PhysicalData::PlanetData& pd1, const PhysicalData::PlanetData& pd2 )
    {
      return pd1.position.norm() < pd2.position.norm(); // TODO, assumes sun in the middle
    } );
    float spaceLength = planetData->position.norm() * 2.0f; // a.u.
    //spaceLength = 50.0f;

    //auto lowerLeft = Eigen::Vector2f( -spaceLength / 2.0f, -spaceLength / 2.0f );
    //auto upperRight = Eigen::Vector2f( spaceLength / 2.0f, spaceLength / 2.0f );
    //spaceBoundingBox_ = Eigen::AlignedBox2f( lowerLeft, upperRight );

    const float screenWidth = mainWindow_.getDefaultView().getSize().x;
    const float screenHeight = mainWindow_.getDefaultView().getSize().y;
    csHandler_ = CoordinateSystemHandler( min( screenWidth, screenHeight ) / spaceLength,
                                          Eigen::Vector2f( screenWidth / 2.0f, screenHeight / 2.0f ) );

    const float height = 10.0f * csHandler_.world2DisplayLength();
    const float width = height * screenWidth / screenHeight;
    view_.setSize( width, height );
    
    buildScene_( planetarySystemData );
    break;
  }
  case GUITests::Setup::RollingBalls:
  {
    dt_ = sf::seconds( 1.0f / 50.0f );
    timeScale_ = 1.0f;

    const float screenWidth = mainWindow_.getDefaultView().getSize().x;
    const float screenHeight = mainWindow_.getDefaultView().getSize().y;
    csHandler_ = CoordinateSystemHandler( 1.0f, Eigen::Vector2f( screenWidth / 2.0f, screenHeight / 2.0f ) );

    // Add one empty SceneNode for each scene layer
    for( size_t i = 0; i < LayerCount; ++i )
    {
      SceneNode::Ptr layer( new SceneNode() );
      sceneLayers_[i] = layer.get();
      sceneGraph_.attachChild( std::move( layer ) );
    }

    vector< PhysicalData::BallData > ballData;
    //PhysicalData::setupMovingAndStationaryBall( ballData );
    //PhysicalData::setupTwoBalls( ballData );
    PhysicalData::setupALotOfBalls( ballData );
    for( auto data : ballData )
    {
      float mass = float( M_PI ) * data.radius * data.radius;
      std::unique_ptr< Ball > ball( new Ball( data.radius, mass, csHandler_ ) );
      ball->setPosition( data.position );
      ball->setVelocity( data.velocity );
      physicalObjects_.push_back( ball.get() );
      balls_.push_back( ball.get() );
      sceneLayers_[ActionLayer]->attachChild( std::move( ball ) );
    }
    break;
  }
  default:
  {
    assert( false );
    break;
  }
  }
}

GUITests::~GUITests()
{}

void GUITests::runSpaceSimulation()
{
  MotionManager motionManager( MotionManager::Integration::RK4,
                               MotionManager::Technique::functor );

  sf::Clock clock;
  sf::Time totalTime;

  CommandQueue commandQueue; // Not used
  sf::Time timeSinceLastUpdate = sf::Time::Zero;
  while( mainWindow_.isOpen() )
  {
    sf::Time elapsedAbsoluteTime = clock.restart();
    timeSinceLastUpdate += elapsedAbsoluteTime * timeScale_;
    while( timeSinceLastUpdate >= dt_ )
    {
      motionManager.computeLinearDynamics( dt_.asSeconds(), physicalObjects_ );

      sceneGraph_.update( dt_, commandQueue );

      timeSinceLastUpdate -= dt_;
      totalTime += dt_;
    }

    sf::Event event;
    while( mainWindow_.pollEvent( event ) )
    {
      switch( event.type ) {
      case sf::Event::Closed:
      {
        mainWindow_.close();
        break;
      }
      }
    }

    //sceneGraph_.update( dt_, commandQueue );

    // Debug
    updateFPSStats_( elapsedAbsoluteTime );

    view_.setCenter( csHandler_.convertToDisplayCS( physicalObjects_[2]->getPosition() ) );
    render_();
  }
}

void GUITests::runCollisionSimulation()
{
  MotionManager motionManager( MotionManager::Integration::RK4,
                               MotionManager::Technique::functor );

  sf::Clock clock;
  sf::Time totalTime;

  CommandQueue commandQueue; // Not used
  sf::Time timeSinceLastUpdate = sf::Time::Zero;
  while( mainWindow_.isOpen() )
  {
    sf::Time elapsedAbsoluteTime = clock.restart();
    timeSinceLastUpdate += elapsedAbsoluteTime * timeScale_;
    while( timeSinceLastUpdate >= dt_ )
    {
      motionManager.computeLinearKinematics( dt_.asSeconds(), physicalObjects_ );
      sceneGraph_.update( dt_, commandQueue );
      //CollisionControl::handleCollisionsBalls( sceneGraph_ );
      CollisionControl::handleCollisions( balls_ );

      timeSinceLastUpdate -= dt_;
      totalTime += dt_;
    }

    sf::Event event;
    while( mainWindow_.pollEvent( event ) )
    {
      switch( event.type ) {
      case sf::Event::Closed:
      {
        mainWindow_.close();
        break;
      }
      }
    }

    // Debug
    updateFPSStats_( elapsedAbsoluteTime );

    render_();
  }
}

void GUITests::loadResources_()
{
  textureHolder_.load( Textures::ID::Mercury, "Resources/Mercury.jpg" );
  textureHolder_.load( Textures::ID::Venus, "Resources/Venus.jpg" );
  textureHolder_.load( Textures::ID::Earth, "Resources/Earth.jpg" );
  textureHolder_.load( Textures::ID::Mars, "Resources/Mars.jpg" );
  textureHolder_.load( Textures::ID::Jupiter, "Resources/Jupiter.jpg" );
  textureHolder_.load( Textures::ID::Saturn, "Resources/Saturn.jpg" );
  textureHolder_.load( Textures::ID::Uranus, "Resources/Uranus.jpg" );
  textureHolder_.load( Textures::ID::Neptune, "Resources/Neptune.jpg" );
  textureHolder_.load( Textures::ID::Sun, "Resources/Sun.jpg" );
  textureHolder_.load( Textures::ID::StarShip, "Resources/StarShip.jpg" );
  textureHolder_.load( Textures::ID::Space, "Resources/Space.jpg" );
  textureHolder_.load( Textures::ID::Rock, "Resources/Rock.png" );

  fontHolder_.load( Fonts::Main, "Resources/garreg.ttf" );
}

void GUITests::buildScene_( const std::vector< PhysicalData::PlanetData >& planetData )
{
  // Add one empty SceneNode for each scene layer
  for( size_t i = 0; i < LayerCount; ++i )
  {
    SceneNode::Ptr layer( new SceneNode() );
    sceneLayers_[i] = layer.get();
    sceneGraph_.attachChild( std::move( layer ) );
  }

  physicalObjects_.clear();
  for( size_t i = 0; i < planetData.size(); ++i )
  {
    const PhysicalData::PlanetData& data = planetData[i];
    unique_ptr< Planet > planet( new Planet( data.name, data.mass, data.radius, 0 /*TODO, borde inte behöva sätta detta*/, textureHolder_, csHandler_ ) );
    planet->setPosition( data.position );
    planet->setVelocity( data.velocity );

    physicalObjects_.push_back( planet.get() );

    sceneLayers_[ActionLayer]->attachChild( std::move( planet ) );
  }
}

void GUITests::render_()
{
  mainWindow_.clear( sf::Color::Black );
  mainWindow_.setView( view_ );

  mainWindow_.draw( sceneGraph_ );

  renderDebugPrintouts_();
  mainWindow_.display();
}

void GUITests::setupDebug_()
{
  debugText_.setFont( fontHolder_.get( Fonts::Main ) );
  debugText_.setCharacterSize( 20 ); // in pixels, not points!
  debugText_.setColor( sf::Color::White );

  fpsText_.setFont( fontHolder_.get( Fonts::Main ) );
  fpsText_.setPosition( float( mainWindow_.getSize().y - 100 ), 0.0f );
  fpsText_.setCharacterSize( 20 );
}

void GUITests::renderDebugPrintouts_()
{
  mainWindow_.setView( view_ );

  if( false )
  {
    // Time
    //float days = totalTime_.asSeconds();
    //string timeString = "Elapsed time = " + to_string( days ) + " days\n";

    string sceneGraphStr = "SceneGraph:\n";
    sceneGraph_.print( sceneGraphStr );

    string debugString = /*timeString +*/ sceneGraphStr;
    debugText_.setString( debugString );
    mainWindow_.draw( debugText_ );
  }

  mainWindow_.draw( fpsText_ );
}

void GUITests::updateFPSStats_( const sf::Time& elapsedTime )
{
  fpsUpdateTime_ += elapsedTime;
  fpsNumFrames_ += 1;

  if( fpsUpdateTime_ >= sf::seconds( 1.0f ) )
  {
    fpsText_.setString(
      "Frames / Second = " + to_string( fpsNumFrames_ ) + "\n" +
      "Time / Update = " + to_string( fpsUpdateTime_.asMicroseconds() / fpsNumFrames_ ) + "us" );

    fpsUpdateTime_ -= sf::seconds( 1.0f );
    fpsNumFrames_ = 0;
  }
}