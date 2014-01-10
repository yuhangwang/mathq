////////////////////////////////////////////////////////////////////////////////
// File: exponential_distribution.c                                           //
// Routine(s):                                                                //
//    Exponential_Distribution                                                //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>              // required for exp()

////////////////////////////////////////////////////////////////////////////////
// double Exponential_Distribution( double x )                                //
//                                                                            //
//  Description:                                                              //
//     This function returns the probability that a random variable with      //
//     a standard Exponential distribution has a value less than "x".         //
//                                                                            //
//  Arguments:                                                                //
//     double x   Argument of Pr[X < x] where X ~ Exp(1).                     //
//                                                                            //
//  Return Values:                                                            //
//     The probability of observing a value less than (or equal) to x assuming//
//     an exponential distribution with mean 1.                               //
//                                                                            //
//  Example:                                                                  //
//     double x, pr;                                                          //
//                                                                            //
//    pr = Exponential_Distribution(x);                                       //
////////////////////////////////////////////////////////////////////////////////
        
double Exponential_Distribution( double x )
{
   return (x <= 0.0) ? 0.0 : 1.0 - exp(-x);
}
