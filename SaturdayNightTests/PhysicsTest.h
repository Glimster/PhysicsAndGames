#pragma once
#include <array>
#include "PhysicalData.h"
#include <gtest/gtest.h>
#include "MotionManager.h"
//class PhysicalObject;

class PhysicsTest : public ::testing::TestWithParam< MotionManager::Integration >
{
public:

  //PhysicsTest();

private:

};