#pragma once

namespace Textures
{
  enum class ID { // Planets
                  Mercury, 
                  Venus, 
                  Earth, 
                  Mars, 
                  Jupiter, 
                  Saturn, 
                  Uranus, 
                  Neptune,
                  // Stars
                  Sun,
                  // Ships
                  StarShip,
                  // Backgrounds
                  Space,
                  // Projectiles
                  Rock };
}

namespace Fonts
{
  enum ID
  {
    Main,
  };
}

// TODO, varför behövs dessa? Förstår icke!
// Forward declarations
namespace sf
{
  class Texture;
  class Font;
}

template < typename Resource, typename Identifier >
class ResourceHolder;

typedef ResourceHolder< sf::Texture, Textures::ID > TextureHolder;
typedef ResourceHolder< sf::Font, Fonts::ID > FontHolder;
