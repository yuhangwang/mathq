////////////////////////////////////////////////////////////////////////////////
// File: t2_distribution.c                                                    //
// Routine(s):                                                                //
//    t2_Distribution                                                         //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>              // required for sqrt()

////////////////////////////////////////////////////////////////////////////////
// double t2_Distribution( double x )                                         //
//                                                                            //
//  Description:                                                              //
//     This function returns the probability that a random variable with      //
//     a t-distribution with two degrees of freedom has a value less than "x".//
//                                                                            //
//  Arguments:                                                                //
//     double x   Argument of Pr[X < x] where X ~ (1 + x / sqrt(2+x^2) ) / 2. //
//                                                                            //
//  Return Values:                                                            //
//     The probability of observing a value less than (or equal) to x assuming//
//     a t-distribution with two degrees of freedom.                          //
//                                                                            //
//  Example:                                                                  //
//     double x, pr;                                                          //
//                                                                            //
//    pr = t2_Distribution(x);                                                //
////////////////////////////////////////////////////////////////////////////////

#define sqrt_2_over_DBL_EPSILON 9.4906265624251559e+07
        
double t2_Distribution( double x )
{
   if (x > sqrt_2_over_DBL_EPSILON) return 1.0;
   if (x < -sqrt_2_over_DBL_EPSILON) return 0.0;

   return 0.5 * ( 1.0 + x / sqrt(2.0 + x * x) );
}
