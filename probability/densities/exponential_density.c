////////////////////////////////////////////////////////////////////////////////
// File: exponential_density.c                                                //
// Routine(s):                                                                //
//    Exponential_Density                                                     //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>              // required for exp()

////////////////////////////////////////////////////////////////////////////////
// double Exponential_Density( double x )                                     //
//                                                                            //
//  Description:                                                              //
//     This function returns the value of the probability density at x where  //
//     x has a standard Exponential distribution.                             //
//                                                                            //
//  Arguments:                                                                //
//     double x   Argument of dPr[X < x]/dx where X ~ Exp(1).                 //
//                                                                            //
//  Return Values:                                                            //
//     The probability density evaluated at X = x where X is assumed to have  //
//     an exponential distribution with mean 1.                               //
//                                                                            //
//  Example:                                                                  //
//     double x, pr;                                                          //
//                                                                            //
//     pr = Exponential_Density(x);                                           //
////////////////////////////////////////////////////////////////////////////////

double Exponential_Density( double x )
{
   return (x < 0.0) ? 0.0 : exp(-x);
}
