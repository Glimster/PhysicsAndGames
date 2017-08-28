#include "stdafx.h"
#include "Player.h"

#include "Category.h"
#include "CommandQueue.h"

#include "StarShip.h"

using namespace std;

Player::Player():
  keyBinding_(),
  actionBinding_()
{
  keyBinding_[sf::Keyboard::Up] = AftThrusters;
  keyBinding_[sf::Keyboard::Left] = LeftThrusters;
  keyBinding_[sf::Keyboard::Right] = RightThrusters;
  keyBinding_[sf::Keyboard::Space] = Fire;

  initializeActions_();
  
  for( auto it = actionBinding_.begin(); it != actionBinding_.end(); ++it )
  {
    it->second.category = Category::Player;
  }
}

void Player::handleEvent( const sf::Event& event, CommandQueue& commands )
{
  switch( event.type )
  {
    case sf::Event::KeyPressed:
    {
      // Check if pressed key appears in key binding, trigger command if so
      auto found = keyBinding_.find( event.key.code );
      if( found != keyBinding_.end() && !isRealtimeAction_( found->second ) )
      {
        Command command = actionBinding_[found->second];
        command.action = bind( actionBinding_[found->second].actionP, placeholders::_1, placeholders::_2, true );
        commands.push( command );
      }
      break;
    }
    case sf::Event::KeyReleased:
    {
      auto found = keyBinding_.find( event.key.code );
      if( found != keyBinding_.end() && !isRealtimeAction_( found->second ) )
      {
        Command command = actionBinding_[found->second];
        command.action = bind( actionBinding_[found->second].actionP, placeholders::_1, placeholders::_2, false );
        commands.push( command );
      }
      break;
    }
  }
}

void Player::handleRealtimeInput( CommandQueue& commands )
{
  // Traverse all assigned keys and check if they are pressed
  for each ( auto& pair in keyBinding_ )
  {
    // If key is pressed, lookup action and trigger corresponding command
    if( sf::Keyboard::isKeyPressed( pair.first ) && isRealtimeAction_( pair.second ) )
      commands.push( actionBinding_[pair.second] );
  }
}

void Player::initializeActions_()
{
  // Thruster actions
  actionBinding_[AftThrusters].actionP = derivedActionP< StarShip >( []( StarShip& ship, sf::Time dt, bool isPressed ) { ship.setAftThrusters( isPressed ); } );
  actionBinding_[LeftThrusters].actionP = derivedActionP< StarShip >( []( StarShip& ship, sf::Time dt, bool isPressed ) { ship.setLeftRotationThrusters( isPressed ); } );
  actionBinding_[RightThrusters].actionP = derivedActionP< StarShip >( []( StarShip& ship, sf::Time dt, bool isPressed ) { ship.setRightRotationThrusters( isPressed ); } );

  actionBinding_[Fire].action = derivedAction< StarShip >( std::bind( &StarShip::fire, placeholders::_1 ) );
}

bool Player::isRealtimeAction_( Action action )
{
  switch( action )
  {
    case Fire:
      return true;

    default:
      return false;
  }
}