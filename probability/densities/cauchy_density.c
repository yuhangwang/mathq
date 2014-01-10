////////////////////////////////////////////////////////////////////////////////
// File: cauchy_density.c                                                     //
// Routine(s):                                                                //
//    Cauchy_Density                                                          //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                   // required for M_PI
                                    
#define PI M_PI

static double const r_pi = 1.0 / PI;

////////////////////////////////////////////////////////////////////////////////
// double Cauchy_Density( double x )                                          //
//                                                                            //
//  Description:                                                              //
//     The density of the Cauchy distribution is                              //
//                          f(x) =  1 / [PI * ( 1 + x^2)].                    //
//                                                                            //
//  Arguments:                                                                //
//     double x   Argument of the density function.  If x <= 0, the result    //
//                is 0 and if x >= 1, the result is 0, otherwise the result   //
//                is the function above except if x = 0 when a < 1 or x = 1   //
//                when b < 1 in which case DBL_MAX is returned.               //
//                                                                            //
//  Return Values:                                                            //
//     A real number between 0 and 1 / PI.                                    //
//                                                                            //
//  Example:                                                                  //
//     double p, x;                                                           //
//                                                                            //
//     p = Cauchy_Density(x);                                                 //
////////////////////////////////////////////////////////////////////////////////

double Cauchy_Density(double x)
{
   return r_pi / (1.0 + x * x);
}

