////////////////////////////////////////////////////////////////////////////////
// File: t2_density.c                                                         //
// Routine(s):                                                                //
//    t2_Density                                                              //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>              // required for sqrt()

////////////////////////////////////////////////////////////////////////////////
// double t2_Density( double x )                                              //
//                                                                            //
//  Description:                                                              //
//     This function returns the value of the probability density at x where  //
//     X has a t-distribution with two degrees of freedom.                    //
//                                                                            //
//  Arguments:                                                                //
//     double x   Argument of dPr[X < x]/dx where X ~ (1 + x/sqrt(2+x^2)) / 2.//
//                                                                            //
//  Return Values:                                                            //
//     The probability density evaluated at X = x where X is assumed to have  //
//     a t-distribution with two degrees of freedom.                          //
//                                                                            //
//  Example:                                                                  //
//     double x, pr;                                                          //
//                                                                            //
//     pr = t2_Density(x);                                                    //
////////////////////////////////////////////////////////////////////////////////

#define ONE_over_cbrt_DBL_MIN 3.5553731598732436e+102

double t2_Density( double x )
{
   double p;

   if (fabs(x) >= ONE_over_cbrt_DBL_MIN) return 0.0;

   p = 1.0 / (2.0 + x * x);
   return p * sqrt(p);
}
