#pragma once
#include "ResourceIdentifiers.h"

struct StarShipData
{
  float mass;
  float length;
  int hitPoints;
  Textures::ID texture;
  float cannonBallLaunchSpeed;
};

struct ProjectileData
{
  float mass;
  float radius;
  int hitPoints;
  int damage;
  Textures::ID texture;
};

const std::vector< StarShipData > setupStarShipData();
const std::vector< ProjectileData > setupProjectileData();
