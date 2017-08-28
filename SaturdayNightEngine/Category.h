#pragma once

// Categories determines who should recieve commands
namespace Category
{
  // Bitwise identification
  // Enables combination of types using the bitwise OR operator |
  enum Type
  {
    None = 0,
    Scene = 1 << 0,
    Player = 1 << 1,
    Projectile = 1 << 2,
    Planet = 1 << 3,
    Ball = 1 << 4,
  };
}