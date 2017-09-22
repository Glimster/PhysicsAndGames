// More or less a copy of the Numerical Recipes ODE solver

#pragma once
#include< vector >

class NR_OdeOutput
{
public:
  // No output saved
  NR_OdeOutput();
               
  // Values saved in nSaved number of points. 
  // nSaved == -1 => values saved only at integration limits x1, x2
  NR_OdeOutput( const int nSaved );

  void init( const int neqn,
             const double x1,
             const double x2 );

  void resize();

  template< typename Stepper >
  void saveDense( Stepper& s, const double xout, const double h )
  {
    if( count_ == kmax_ )
      resize();
    for( int i = 0; i < nvar_; ++i )
      ysave_( i, count_ ) = s.dense_out( i, xout, h );
    xsave_[count_++] = xout;
  }
  
  void save( const double x, const std::vector< double >& y );

  template< typename Stepper >
  void out( const int nstp, const double x, const std::vector< double >& y, Stepper& s, const double h )
  {
    if( !dense_ )
      throw runtime_error( "Dense output not set in NR_OdeOutput" );
    if( nstp == -1 )
    {
      save( x, y );
      xout_ += dxout_;
    }
    else
    {
      while( (x - xout_) * (x2_ - x1_) > 0.0 )
      {
        saveDense( s, xout_, h );
        xout_ += dxout_;
      }
    }
  }

public:
  int kmax_;
  int nvar_;
  int nSaved_;
  bool dense_;
  int count_;
  double x1_, x2_, xout_, dxout_;
  std::vector< double > xsave_;
  Eigen::MatrixXd ysave_;

};

