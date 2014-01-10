////////////////////////////////////////////////////////////////////////////////
// File: cauchy_distribution.c                                                //
// Routine(s):                                                                //
//    Cauchy_Distribution                                                     //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                    // required for atan() and M_PI
                                    
#define PI M_PI

static double const r_pi = 1.0 / PI;

////////////////////////////////////////////////////////////////////////////////
// double Cauchy_Distribution( double x )                                     //
//                                                                            //
//  Description:                                                              //
//     The Cauchy distribution is the integral from -inf to x of the density  //
//                          f(x) =  1 / [PI * ( 1 + x^2)],                    //
//     i.e.                                                                   //
//                   Pr[X < x] = F(x) = 1 / 2 + (1 / PI) arctan(x).           //
//                                                                            //
//  Arguments:                                                                //
//     double x   Argument of the Cauchy distribution Pr[X < x].              //
//                                                                            //
//  Return Values:                                                            //
//     A real number between 0 and 1.                                         //
//                                                                            //
//  Example:                                                                  //
//     double p, x;                                                           //
//                                                                            //
//     p = Cauchy_Distribution(x);                                            //
////////////////////////////////////////////////////////////////////////////////

double Cauchy_Distribution(double x)
{
   return 0.5 + r_pi * atan(x);
}
