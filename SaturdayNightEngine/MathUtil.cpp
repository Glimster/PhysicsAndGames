#include "stdafx.h"
#include "MathUtil.h"

using namespace std;

const float MathUtil::oneOverSix_ = 1.0f / 6.0f;

float MathUtil::fociCenterDistanceEllipse( float a, float b )
{
  return sqrt( a * a - b * b );
}

float MathUtil::areaEllipse( float a, float b )
{
  return float( M_PI ) * a * b;
}
float MathUtil::circumferenceEllipse( float a, float b )
{
  return float( M_PI ) * ( 3 * ( a + b ) - sqrt( ( 3 * a + b ) * ( a + 3 * b ) ) );
}