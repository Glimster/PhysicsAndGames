#pragma once

class SceneNode;

struct Command
{
public:
  Command();
  
  // What to do
  std::function< void( SceneNode&, sf::Time ) > action;

  // Temporary action with KeyPressed/Released parameter
  // TODO, blev kanske inte så bra lösning...
  std::function< void( SceneNode&, sf::Time, bool ) > actionP;
  
  // Who should recieve the command
  unsigned int category;
};

// Adapter for passing any derived (from SceneNode) object into a command.action
template< typename GameObject, typename Function >
std::function< void( SceneNode&, sf::Time ) > derivedAction( Function fn )
{
  return [=]( SceneNode& node, sf::Time dt )
  {
    assert( dynamic_cast< GameObject* >( &node ) != nullptr );

    fn( static_cast< GameObject& >( node ), dt );
  };
}

// TODO, blev kanske inte så bra lösning...
template< typename GameObject, typename Function >
std::function< void( SceneNode&, sf::Time, bool ) > derivedActionP( Function fn )
{
  return [=]( SceneNode& node, sf::Time dt, bool isPressed )
  {
    assert( dynamic_cast< GameObject* >(&node) != nullptr );

    fn( static_cast< GameObject& >(node), dt, isPressed );
  };
}