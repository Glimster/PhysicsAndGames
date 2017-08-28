#include "stdafx.h"
#include "StarShip.h"

#include "GameDataTables.h"
#include "GameEntity.h"
#include "Category.h"
#include "CoordinateSystemHandler.h"
#include "ResourceHolder.h"
#include "Physics.h"
#include "Utilities.h"
#include "TextNode.h"
#include "SpriteNode.h"

#include "CommandQueue.h"

using namespace std;

// Anonymous namespace to avoid name collision
namespace
{
  const std::vector< StarShipData > Table = setupStarShipData();
  const std::vector< ProjectileData > PTable = setupProjectileData();
}

StarShip::StarShip( Type type,
                    const TextureHolder& textureHolder,
                    const FontHolder& fonts,
                    const CoordinateSystemHandler& csHandler,
                    vector< PhysicalObject* >& physicalObjects ):
  GameEntity( Table[type].mass, Table[type].hitPoints, csHandler ),
  type_( type ),
  sprite_( textureHolder.get( Table[type].texture ) ),
  aftThrusters_( false ),
  leftRotationThrusters_( false ),
  rightRotationThrusters_( false ),
  fireCommand_(),
  isFiring_( false ),
  fireCountdown_( sf::Time::Zero ),
  fireRate_( 10.0f ),
  cannonBallLaunchSpeed_( Table[type].cannonBallLaunchSpeed ),
  cannonBallMass_( PTable[Projectile::CannonBall].mass ),
  cannonBall_( nullptr )
{
  const float shipWidth = Table[type].length * sprite_.getTextureRect().height / sprite_.getTextureRect().width;
  const float scaleX = Table[type].length * csHandler_.world2DisplayLength() / sprite_.getTextureRect().width;
  const float scaleY = shipWidth * csHandler_.world2DisplayLength() / sprite_.getTextureRect().height;
  transformable_.setScale( scaleX, scaleY );

  Utilities::centerOrigin( sprite_ );

  // Approximate ship with rectangle
  momentOfInertia_ = Physics::computeRotationalInertiaForRectangle( Table[type].mass, Table[type].length, shipWidth );

  unique_ptr< TextNode > healthDisplay( new TextNode( fonts, "" ) );
  healthDisplay_ = healthDisplay.get();
  attachChild( std::move( healthDisplay ) );

  
  { // Create ready-to-fire cannon ball texture (not the actual object being launched)
    // TODO, få in skalningen i SpritNode-klassen! Skicka in vilken faktiskt storlek den ska ha istället
    const auto& ballTexture = textureHolder.get( PTable[Projectile::CannonBall].texture ); 
    const float scale = 2.0f * PTable[Projectile::CannonBall].radius * csHandler_.world2DisplayLength() / ballTexture.getSize().x;
    sf::Vector2f scaleFactor( scale / scaleX, scale / scaleY ); // Scaling should not be done twice (since it is a child of the StarShip)
    unique_ptr< SpriteNode > ball( new SpriteNode( ballTexture, scaleFactor ) );
    // Set position relative star ship
    float pos = Table[type_].length * csHandler_.world2DisplayLength() / 2.0f;
    pos /= scaleX; // Scaling should not be done twice (since it is a child of the StarShip)
    ball->setPositionDisplayCS( sf::Vector2f( pos, 0.0f ) );

    cannonBall_ = ball.get();
    attachChild( std::move( ball ) );
  }

  fireCommand_.category = Category::Scene;
  fireCommand_.action = [this, &textureHolder, &physicalObjects] ( SceneNode& node, sf::Time ) 
  { 
    launchCannonBall_( node, Projectile::CannonBall, physicalObjects, textureHolder ); 
  };
}

StarShip::~StarShip()
{}

float StarShip::getMass() const
{ 
  return mass_ + cannonBallMass_ * (fireRate_ - fireCountdown_.asSeconds()) / fireRate_;
}

void StarShip::computeLinearForceOverM( Eigen::Vector2f& fOverM ) const
{
  // TODO, lägg in variabel kraft?
  if( aftThrusters_ )
  {
    const float magnitude = /*getMass() **/ 0.0172f; // Will accelerate to v ~= 0.0172 in dt = 1
    fOverM = orientation_ * magnitude;
  }
  else 
  {
    fOverM.setZero();
  }
}

float StarShip::computeTorque() const
{
  float torque = 0.0f;

  // TODO, lägg in variabel kraft?
  if( leftRotationThrusters_ )
  {
    torque += momentOfInertia_ * float( M_PI ) / 50.0f; // Will accelerate to w = pi/50 in dt = 1
  }
  if( rightRotationThrusters_ )
  {
    torque -= momentOfInertia_ * float( M_PI ) / 50.0f;
  }

  return torque;
}

sf::FloatRect StarShip::getBoundingRectDisplayCS() const
{
  auto rect = sprite_.getGlobalBounds();
  rect.width *= 0.5f; // Exclude front end to be able to launch projectile without collision. TODO, make something better!
  return getWorldTransformDisplayCS().transformRect( rect );
}

void StarShip::fire()
{
  isFiring_ = true;
}

unsigned int StarShip::getCategory() const
{
  return Category::Player;
}

void StarShip::drawCurrent_( sf::RenderTarget& target, sf::RenderStates states ) const
{
  target.draw( sprite_, states );
}

void StarShip::updateCurrent_( sf::Time dt, CommandQueue& commands )
{
  checkProjectileLaunch_( dt, commands );
  
  transformable_.setPosition( csHandler_.convertToDisplayCS( position_ ) );
  transformable_.setRotation( csHandler_.computeAngleInDisplayCS( orientation_ ) );

  {
    const float scale = 2.0f * PTable[Projectile::CannonBall].radius * csHandler_.world2DisplayLength() / cannonBall_->getTextureRect().width;
    sf::Vector2f scaleFactor( scale / transformable_.getScale().x, scale / transformable_.getScale().y ); // Scaling should not be done twice (since it is a child of the StarShip)
    scaleFactor.x *= (fireRate_ - fireCountdown_.asSeconds()) / fireRate_;
    scaleFactor.y *= (fireRate_ - fireCountdown_.asSeconds()) / fireRate_;
    cannonBall_->setScaleDisplayCS( scaleFactor );
  }

  healthDisplay_->setString( to_string( getHitPoints() ) + " HP" );
  healthDisplay_->setPositionDisplayCS( 0.f, 500.f ); // TODO, gör att texten snurrar runt skeppet. Fixa!
  healthDisplay_->setRotationDisplayCS( -transformable_.getRotation() );
}

void StarShip::printCurrent_( string& string ) const
{
  string.append( "StarShip: category = " + to_string( getCategory() ) );
  printDefaultData_( string );
  printPOData_( string );
  string.append( ", at = " + to_string( aftThrusters_ ) +
                 ", lt = " + to_string( leftRotationThrusters_ ) +
                 ", rt = " + to_string( rightRotationThrusters_ ) );
  string.append( "\n" );
}

void StarShip::checkProjectileLaunch_( sf::Time dt, CommandQueue& commands )
{
  if( isFiring_ && fireCountdown_ <= sf::Time::Zero )
  {
    commands.push( fireCommand_ );
    fireCountdown_ += sf::seconds( fireRate_ );
    isFiring_ = false;
  }
  else if( fireCountdown_ > sf::Time::Zero )
  {
    fireCountdown_ -= dt;
    isFiring_ = false;
  }
}
  
void StarShip::launchCannonBall_( SceneNode& node,
                                  Projectile::Type type,
                                  vector< PhysicalObject* >& physicalObjects,
                                  const TextureHolder& textureHolder )
{
  std::unique_ptr< Projectile > projectile( new Projectile( type, textureHolder, csHandler_ ) );
  projectile->setPosition( getPosition() + Table[type_].length / 2.0f * getOrientation() );
  projectile->setVelocity( getVelocity() + cannonBallLaunchSpeed_ * getOrientation() ); // TODO, not yet relativistic!

  physicalObjects.push_back( projectile.get() );
  
  // Conservation of momentum
  const float massBefore = mass_ + cannonBallMass_; // Will be excecuted after mass has been reduced!
  const float massAfter = mass_;
  Eigen::Vector2f velocityAfter = (getVelocity() * massBefore - projectile->getVelocity() * projectile->getMass()) / massAfter;
  setVelocity( velocityAfter );

  node.attachChild( std::move( projectile ) );
}
