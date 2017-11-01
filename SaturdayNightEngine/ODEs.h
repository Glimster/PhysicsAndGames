#pragma once
#include <vector>

#include "PhysicalConstants.h"

// ----- Functors

// van der Pol's equation (oscillator)
// https://en.wikipedia.org/wiki/Van_der_Pol_oscillator
struct VanDerPol {
  VanDerPol( double eps ) :
    epsilon_( eps )
  {}

  void operator() ( const double x,
                    const std::vector< double >& y,
                    std::vector< double >& dydx )
  {
    dydx[0] = y[1];
    dydx[1] = ((1.0 - y[0] * y[0]) * y[1] - y[0]) / epsilon_;
  }

private:
  const double epsilon_;
};

struct FixedGravitationalField {
  FixedGravitationalField( const double g ) :
    g_( g )
  {}

  void operator() ( const double x,
                    const std::vector< double >& y,
                    std::vector< double >& dydx )
  {
    dydx[0] = y[2]; // dr_x/dt = v_x
    dydx[1] = y[3]; // dr_y/dt = v_y
    dydx[2] = 0.0;  // dv_x/dt = F / m = 0.0
    dydx[3] = -g_;  // dv_y/dt = F / m = -g
  }

private:
  const double g_;
};

struct HookesLaw {
  HookesLaw( const double k, const double m ) :
    k_( k ),
    m_( m )
  {}

  void operator() ( const double x,
                    const std::vector< double >& y,
                    std::vector< double >& dydx )
  {
    dydx[0] = y[1];
    dydx[1] = -k_ * y[0] / m_;
  }

private:
  const double k_;
  const double m_;
};

struct GravitationalNBody {
  GravitationalNBody( const std::vector< double >& masses ):
    masses_( masses ),
    kSquared_( double( Phys::PhysicalConstants::kSquared ) )
  {
    softeningLength2_ = Phys::PhysicalConstants::softeningLength * Phys::PhysicalConstants::softeningLength;
  }

  inline void operator() ( const double x,
                           const std::vector< double >& y,
                           std::vector< double >& dydx )
  {
    // TODO: hur är det med cache-missar och masses_vektorn!? Den ligger väl inte i anknytning
    // till y i minnet!?!?!? Är det ok om två olika vektorer som används ligger sekventiellt i minnet
    // var för sig, men inte tillsammans?
    // Och hur är det med konstanter, t.ex. kSquared?
    for( size_t i = 0; i < y.size(); i += 4 )
    {
      double FOverMx = 0.0, FOverMy = 0.0;
      for( size_t j = 0; j < y.size(); j += 4 )
      {
        if( i == j )
          continue;

        const double rx = y[i] - y[j];
        const double ry = y[i + 1] - y[j + 1];

        // TODO
        const double r2 = rx * rx + ry * ry; // + softeningLength2_;
        const double rNorm = sqrt( r2 );
        const double rxNorm = rx / rNorm;
        const double ryNorm = ry / rNorm;
        const double force = kSquared_ * masses_[j / size_t( 4 )] / r2;
        FOverMx += -force * rxNorm;
        FOverMy += -force * ryNorm;
      }

      dydx[i] = y[i + 2];
      dydx[i + 1] = y[i + 3];
      dydx[i + 2] = FOverMx;
      dydx[i + 3] = FOverMy;
    }
  }

private:
  const std::vector< double >& masses_;
  double softeningLength2_;

  const double kSquared_;
};