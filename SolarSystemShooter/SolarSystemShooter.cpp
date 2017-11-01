// SolarSystemShooter.cpp : Defines the entry point for the console application.

#include "stdafx.h"
#include "Game.h"

using namespace std;

int _tmain( int argc, _TCHAR* argv[] )
{
  try {
    Game game;
    game.run();
  }
  catch( exception& e )
  {
    cout << "\nEXCEPTION: " << e.what() << endl;
  }
}