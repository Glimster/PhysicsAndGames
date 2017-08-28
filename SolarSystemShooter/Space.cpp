#include "stdafx.h"
#include "Space.h"

#include "PhysicalObject.h"
#include "SpriteNode.h"
#include "Planet.h"
#include "StarShip.h"

#include "PhysicalData.h"
#include "PhysicalConstants.h"
#include "MathUtil.h"
#include "Integrator.h"
#include "CollisionControl.h"
using namespace std;

Space::Space( sf::RenderWindow& mainWindow, FontHolder& fontHolder ):
  mainWindow_( mainWindow ),
  textureHolder_(),
  fontHolder_( fontHolder ),
  spaceBoundingBox_(),
  csHandler_(),
  playerView_(),
  motionManager_( MotionManager::Integration::RK4, MotionManager::Technique::functor ),
  commandQueue_(),
  sceneLayers_(),
  sceneGraph_(),
  physicalObjects_(),
  player_()
{
  vector< PhysicalData::PlanetData > planetarySystemData;
  PhysicalData::setupRealisticSolarSystem( planetarySystemData );
  //PhysicalData::setupPlanetarySystem( planetarySystemData );
  
  {
    auto planetData = max_element( planetarySystemData.begin(), planetarySystemData.end(), []
    ( const PhysicalData::PlanetData& pd1, const PhysicalData::PlanetData& pd2 ) 
    {
      return pd1.position.norm() < pd2.position.norm(); // TODO, assumes sun in the middle
    } );
    const float spaceLength = planetData->position.norm() * 2.0f; // a.u.
    
    auto lowerLeft = Eigen::Vector2f( -spaceLength / 2.0f, -spaceLength / 2.0f );
    auto upperRight = Eigen::Vector2f( spaceLength / 2.0f, spaceLength / 2.0f );
    spaceBoundingBox_ = Eigen::AlignedBox2f( lowerLeft, upperRight );

    const float screenWidth = mainWindow_.getDefaultView().getSize().x;
    const float screenHeight = mainWindow_.getDefaultView().getSize().y;
    csHandler_ = CoordinateSystemHandler( min( screenWidth, screenHeight ) / spaceLength,
                                          Eigen::Vector2f( screenWidth / 2.0f, screenHeight / 2.0f ) );
    // Setup initial view
    playerView_ = mainWindow_.getDefaultView();
    // TODO, which size?
    const float height = 2.0f * csHandler_.world2DisplayLength();
    const float width = height * screenWidth / screenHeight;
    playerView_.setSize( width, height );
  }

  loadTextures_();
  fontHolder_.load( Fonts::Main, "Resources/garreg.ttf" ); // TODO, borde egentligen läsas in i Game, men då kraschar det eftersom den inte finns nu isf...
  buildScene_( planetarySystemData );
}

Space::~Space()
{}

void Space::update( const sf::Time& timeStep )
{
  const float dt = timeStep.asSeconds();

  // TODO, frågan är om jag kan stoppa den här nån annanstans, t.ex. utanför inre fysik-loopen,
  // eftersom jag inte använder dt!?!?!?
  while( !commandQueue_.isEmpty() )
  {
    sceneGraph_.onCommand( commandQueue_.pop(), timeStep );
  }

  sceneGraph_.update( timeStep, commandQueue_ );

  CollisionControl::handleCollisions( sceneGraph_ );

  motionManager_.computeLinearDynamics( dt, physicalObjects_ );
  motionManager_.computeRotationalDynamics( dt, *player_ );
}

void Space::render()
{
  mainWindow_.setView( playerView_ );
  mainWindow_.draw( sceneGraph_ );
}

void Space::loadTextures_()
{
  textureHolder_.load( Textures::ID::Mercury,   "Resources/Mercury.jpg" );
  textureHolder_.load( Textures::ID::Venus,     "Resources/Venus.jpg" );
  textureHolder_.load( Textures::ID::Earth,     "Resources/Earth.jpg" );
  textureHolder_.load( Textures::ID::Mars,      "Resources/Mars.jpg" );
  textureHolder_.load( Textures::ID::Jupiter,   "Resources/Jupiter.jpg" );
  textureHolder_.load( Textures::ID::Saturn,    "Resources/Saturn.jpg" );
  textureHolder_.load( Textures::ID::Uranus,    "Resources/Uranus.jpg" );
  textureHolder_.load( Textures::ID::Neptune,   "Resources/Neptune.jpg" );

  textureHolder_.load( Textures::ID::Sun,       "Resources/Sun.jpg" );

  textureHolder_.load( Textures::ID::StarShip,  "Resources/StarShip.jpg" );

  textureHolder_.load(Textures::ID::Space,      "Resources/Space.jpg");

  textureHolder_.load( Textures::ID::Rock,      "Resources/Rock.png" );
}

void Space::buildScene_( const std::vector< PhysicalData::PlanetData >& planetarySystemData )
{
  // Add one empty SceneNode for each scene layer
  for( size_t i = 0; i < LayerCount; ++i )
  {
    Category::Type category = (i == ActionLayer) ? Category::Scene : Category::None;

    SceneNode::Ptr layer( new SceneNode( category ) );
    sceneLayers_[i] = layer.get();
    sceneGraph_.attachChild( std::move( layer ) );
  }

  { // Background
    sf::Texture& texture = textureHolder_.get( Textures::ID::Space );
    
    // Relate to current view. TODO, support zooming etc.
    // TODO, kommer tryckas ihop lite...
    auto scaleFactor = sf::Vector2f( playerView_.getSize().x / texture.getSize().x,
                                     playerView_.getSize().y / texture.getSize().y );
    // TODO, support repeated texture
    //texture.setRepeated( true );
    //sf::IntRect textureRect( 0, 0, ?, ? );
    //unique_ptr< SpriteNode > background( new SpriteNode( texture, textureRect, scaleFactor ) );
    unique_ptr< SpriteNode > background( new SpriteNode( texture, scaleFactor ) );

    Eigen::Vector2f origo( 0.0f, 0.0f );
    background->setPositionDisplayCS( csHandler_.convertToDisplayCS( origo ) );
    sceneLayers_[Background]->attachChild( std::move( background ) );
  }
  
  for( size_t i = 0; i < planetarySystemData.size(); ++i )
  {
    const PhysicalData::PlanetData& data = planetarySystemData[i];
    unique_ptr< Planet > planet( new Planet( data.name, data.mass, data.radius, 1000 /*TODO, ska planeter ha hit points?*/, textureHolder_, csHandler_ ) );
    planet->setPosition( data.position );
    planet->setVelocity( data.velocity );
    
    physicalObjects_.push_back( planet.get() );

    sceneLayers_[ActionLayer]->attachChild( std::move( planet ) );
  }

  { // Create starship
    unique_ptr< StarShip > player( new StarShip( StarShip::Standard, textureHolder_, fontHolder_, csHandler_, physicalObjects_ ) );
    player->setPosition( Eigen::Vector2f( 1.0f, 1.0f ) );
    player->setVelocity( Eigen::Vector2f( 0.0f, 0.0f ) );
    player->setOrientation( Eigen::Vector2f( 1.0f, 0.0f ) );

    player_ = player.get();

    physicalObjects_.push_back( player.get() );
    sceneLayers_[ActionLayer]->attachChild( std::move( player ) );
  }
}