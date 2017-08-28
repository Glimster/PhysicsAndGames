#include "stdafx.h"
#include "GameDataTables.h"

#include "StarShip.h"
#include "Projectile.h"

using namespace std;

const vector< StarShipData > setupStarShipData()
{
  vector< StarShipData > data( StarShip::TypeCount );
  data[StarShip::Standard].hitPoints = 100;
  data[StarShip::Standard].mass = 1e-6f; // Solar masses
  data[StarShip::Standard].length = 0.17f; // a.u.
  data[StarShip::Standard].texture = Textures::ID::StarShip;
  data[StarShip::Standard].cannonBallLaunchSpeed = 0.05f; // a.u. / solar day

  return data;
}

const vector< ProjectileData > setupProjectileData()
{
  vector< ProjectileData > data( Projectile::TypeCount );
  data[Projectile::CannonBall].hitPoints = 10;
  data[Projectile::CannonBall].mass = 1e-7f; // Solar masses
  data[Projectile::CannonBall].radius = 0.02f; // a.u.
  data[Projectile::CannonBall].texture = Textures::ID::Rock;

  return data;
}