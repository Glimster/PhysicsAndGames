#include "stdafx.h"
#include "Game.h"

#include "PhysicalObject.h"
#include "Planet.h"
#include "StarShip.h"

#include "PhysicalConstants.h"

using namespace std;

Game::Game(void):
  mainWindow_( sf::VideoMode::getDesktopMode(), "Solar System" ),
  GUIView_( mainWindow_.getDefaultView() ),
  fontHolder_(),
  space_( mainWindow_, fontHolder_ ),
  isPaused_( false ),
  clock_(),
  totalTime_( sf::Time::Zero ),
  maxTimeStep_( sf::seconds( 1.0f / 60.0f ) ), // TODO: how high can this value be set?
  // Debug
  debugText_(),
  fpsText_(),
  fpsUpdateTime_(),
  fpsNumFrames_( 0 )
{
  Phys::PhysicalConstants::timeScale = 10.0f;

  // TODO, egentligen borde fonten läsas in här, inte i Space?
  //fontHolder_.load( Fonts::Main, "Resources/garreg.ttf" );

  // Debug
  setupDebug_();
}

Game::~Game(void)
{}

void Game::run()
{
  clock_ = sf::Clock(); // Start the clock

  // Game logic (and physics simulation) time step
  sf::Time dt = sf::seconds( 1.0f / 60.0f ) * Phys::PhysicalConstants::timeScale;
  dt = dt > maxTimeStep_ ? maxTimeStep_ : dt; // Ensures that physics simulation will not explode
  
  sf::Time timeSinceLastUpdate = sf::Time::Zero;
  while( mainWindow_.isOpen() )
  {
    const sf::Time elapsedAbsoluteTime = clock_.restart();
    timeSinceLastUpdate += elapsedAbsoluteTime * Phys::PhysicalConstants::timeScale;
    
    while( timeSinceLastUpdate >= dt ) {
      //processEvents(); // Behöver denna verkligen göras så ofta? Det här har ju att göra med physik juh!
      if( !isPaused_ ) // Ska man inte göra så att tiden också stannar när det är pausat?
        update( dt ); 

      timeSinceLastUpdate -= dt;
      totalTime_ += dt;
    }
    
    processEvents(); // Är det ok att den här ligger utanför innersta loopen?


    // Debug
    updateFPSStats_( elapsedAbsoluteTime );

    // TODO, add interpolation so that the rendered state always matches the current time!

    render();
  }
}

void Game::processEvents()
{
  CommandQueue& commands = space_.getCommandQueue();

  sf::Event event;
  while( mainWindow_.pollEvent( event ) )
  {
    player_.handleEvent( event, commands );

    switch( event.type )
    {
      case sf::Event::LostFocus:
      {
        isPaused_ = true;
        break;
      }
      case sf::Event::GainedFocus:
      {
        isPaused_ = false;
        break;
      }
      case sf::Event::Closed:
      {
        mainWindow_.close();
        break;
      }
    }
  }

  player_.handleRealtimeInput( commands );
}

void Game::update( const sf::Time& timeStep )
{
  space_.update( timeStep );
}

void Game::render()
{
  mainWindow_.clear( sf::Color::Black );
  space_.render();
  renderDebugPrintouts_( true );
  mainWindow_.display();
}

void Game::setupDebug_()
{
  debugText_.setFont( fontHolder_.get( Fonts::Main ) );
  debugText_.setCharacterSize( 20 ); // in pixels, not points!
  debugText_.setColor( sf::Color::White );

  fpsText_.setFont( fontHolder_.get( Fonts::Main ) );
  fpsText_.setPosition( float( mainWindow_.getSize().y - 100 ), 0.0f );
  fpsText_.setCharacterSize( 20 );
}

void Game::renderDebugPrintouts_( bool onlyFPS )
{
  mainWindow_.setView( GUIView_ );
  
  if( !onlyFPS )
  {
    // Time
    float days = totalTime_.asSeconds();
    string timeString = "Elapsed time = " + to_string( days ) + " days\n";

    const SceneNode& sceneGraph = space_.getSceneGraph();
    string sceneGraphStr = "SceneGraph:\n";
    sceneGraph.print( sceneGraphStr );

    string debugString = timeString + sceneGraphStr;
    debugText_.setString( debugString );
    mainWindow_.draw( debugText_ );
  }

  mainWindow_.draw( fpsText_ );
}

void Game::updateFPSStats_( const sf::Time& elapsedTime )
{
	fpsUpdateTime_ += elapsedTime;
	fpsNumFrames_ += 1;

	if( fpsUpdateTime_ >= sf::seconds( 1.0f ) )
	{
		fpsText_.setString(
			"Frames / Second = " + to_string( fpsNumFrames_ ) + "\n" +
			"Time / Update = " + to_string( fpsUpdateTime_.asMicroseconds() / fpsNumFrames_ ) + "us");
							 
		fpsUpdateTime_ -= sf::seconds( 1.0f );
		fpsNumFrames_ = 0;
	}
}