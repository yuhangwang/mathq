////////////////////////////////////////////////////////////////////////////////
// File: gaussian_distribution.c                                              //
// Routine(s):                                                                //
//    Gaussian_Distribution                                                   //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>              // required for erf() and M_SQRT1_2

////////////////////////////////////////////////////////////////////////////////
// double Gaussian_Distribution( double x )                                   //
//                                                                            //
//  Description:                                                              //
//     This function returns the probability that a random variable with      //
//     a standard Normal (Gaussian) distribution has a value less than "x".   //
//                                                                            //
//  Arguments:                                                                //
//     double x   Argument of Pr[X < x] where X ~ N(0,1).                     //
//                                                                            //
//  Return Values:                                                            //
//     The probability of observing a value less than (or equal) to x assuming//
//     a normal (Gaussian) distribution with mean 0 and variance 1.           //
//                                                                            //
//  Example:                                                                  //
//     double x, pr;                                                          //
//                                                                            //
//    pr = Gaussian_Distribution(x);                                          //
////////////////////////////////////////////////////////////////////////////////
        
double Gaussian_Distribution( double x )
{
   return  0.5 * ( 1.0 + erf( M_SQRT1_2 * x ) );
}
