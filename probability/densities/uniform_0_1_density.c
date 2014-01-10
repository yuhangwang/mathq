////////////////////////////////////////////////////////////////////////////////
// File: uniform_0_1_density.c                                                //
// Routine(s):                                                                //
//    Uniform_0_1_Density                                                     //
////////////////////////////////////////////////////////////////////////////////

#include <float.h>                  // required for DBL_MAX

////////////////////////////////////////////////////////////////////////////////
// double Uniform_0_1_Density( double x )                                     //
//                                                                            //
//  Description:                                                              //
//     The density of the uniform(0,1) distribution is:                       //
//                        f(x) = 1  if 0 < x < 1                              //
//                             = 0  if x < 0 or x > 1                         //
//                        undefined at x = 0 and x = 1.                       //
//                                                                            //
//  Arguments:                                                                //
//     double x   Argument of f(x).                                           //
//                                                                            //
//  Return Values:                                                            //
//     f(x) if x != 0 and x != 1.  If x = 0 or x = 1, then DBL_MAX is         //
//     returned at x = 0 and - DBL_MAX at x = 1.                              //
//                                                                            //
//  Example:                                                                  //
//     double x, pr;                                                          //
//                                                                            //
//     pr = Uniform_0_1_Density(x);                                           //
////////////////////////////////////////////////////////////////////////////////

double Uniform_0_1_Density( double x )
{
   if (x <= 0.0) return (x == 0.0) ? DBL_MAX : 0.0;
   if (x >= 1.0) return (x == 1.0) ? -DBL_MAX : 0.0;
   return  1.0;
}
