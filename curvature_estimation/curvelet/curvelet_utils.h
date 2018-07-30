// This is curvelet_utils.h
#ifndef curvelet_utils_h
#define curvelet_utils_h

/***************************************************************************************
//file: edgemap.h
//brief: utility functions for the curvelet grouping algorithm
//author: Xiaoyan Li
***************************************************************************************/

#include "edgemap.h"
#include <cmath>

const double k_th = 1e-4;

double sq_distance(const point_2d& A, const point_2d& B)
{
  return (A.first-B.first)*(A.first-B.first)+(A.second-B.second)*(A.second-B.second);
}

double sq_norm(const point_2d& A)
{
    return A.first*A.first+A.second*A.second;
}

//: Convert an angle to [0, 2Pi) range
double To2Pi (double angle)
{
  double a;
  if (angle>=2*M_PI)
    a = fmod(angle,2*M_PI);
  else if (angle < 0)
    a = (2*M_PI+fmod(angle,2*M_PI));
  else 
    a= angle;

  // added by Nhon: these two lines of code is to fix the bug when
  // angle = -1.1721201390607859e-016
  // then after all the computation, we get
  // a = 6.2831853071795862 == 2*vnl_math::pi !!!!!!!
  // the only case this can happen is when a is very close to zero.
  if (!(a>=0 && a<2*M_PI)) {
    a = 0;
  }

  return a;
}

//: Convert an angle to [-Pi, Pi) range
double ToPi (double angle)
{
    double a = angle+M_PI;
    a = To2Pi(a);
    
    return a-M_PI;
}

double vPointPoint (point_2d& start, point_2d& end)
{
  return To2Pi (atan2(end.second - start.second, end.first - start.first) );
}

//: dot product between two angle
double dot (double v1, double v2)
{
  return cos(v1)*cos(v2) + sin(v1)*sin(v2);
}

//: rotate a vector by theta
point_2d rotate(const point_2d &pt, double theta)
{
    return point_2d(cos(theta)*pt.first-sin(theta)*pt.second, sin(theta)*pt.first+cos(theta)*pt.second);
}
#endif
