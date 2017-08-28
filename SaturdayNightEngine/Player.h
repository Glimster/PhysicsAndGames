#pragma once
#include "Command.h"

class CommandQueue;

class Player
{
public:
  enum Action
  {
    AftThrusters,
    LeftThrusters,
    RightThrusters,
    Fire,
    ActionCount
  };

public:
  Player();

  void handleEvent( const sf::Event& event, CommandQueue& commands );
  void handleRealtimeInput( CommandQueue& commands );

private:

  void initializeActions_();
  static bool isRealtimeAction_( Action action );

private:
  std::map< sf::Keyboard::Key, Action > keyBinding_; // TODO, expand so that they can be assigned
  std::map< Action, Command > actionBinding_;
};