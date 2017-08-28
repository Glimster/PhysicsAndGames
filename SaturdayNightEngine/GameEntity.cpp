#include "stdafx.h"
#include "GameEntity.h"


GameEntity::GameEntity( float mass, int hitPoints, const CoordinateSystemHandler& csHandler ):
  PhysicalObject( mass ),
  SceneNode(),
  csHandler_( csHandler ),
  hitPoints_( hitPoints )
{}

GameEntity::~GameEntity()
{}
