#include "stdafx.h"
#include "Massive.h"


Massive::Massive( float mass ):
  mass_( mass ),
  momentOfInertia_( 0.0f ),
  angularMomentum_( 0.0f )
{}

Massive::~Massive()
{}

void Massive::updateAngularMomentum( float dL )
{
  angularMomentum_ += dL;
}
