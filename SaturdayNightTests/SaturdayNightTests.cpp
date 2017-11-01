// SaturdayNightTests.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <gtest/gtest.h>

#include "GUITests.h"
#include "PhysicsTests.h"

using namespace std;

int _tmain( int argc, _TCHAR* argv[] )
{
  //PhysicsTests< int >::compileMe();

  // TODO, ska det vara separata projekt?
#if 0
  GUITests test;
  //test.init( GUITests::Setup::Space );
  //test.runSpaceSimulation();
  test.init( GUITests::Setup::RollingBalls );
  test.runCollisionSimulation();
#endif
#if 1

  testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
#endif
#if 0
  // For profiling
  PhysicsTests::test();
#endif
}