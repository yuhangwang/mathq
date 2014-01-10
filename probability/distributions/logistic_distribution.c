////////////////////////////////////////////////////////////////////////////////
// File: logistic_distribution.c                                              //
// Routine(s):                                                                //
//    Logistic_Distribution                                                   //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                                 // required for exp()

////////////////////////////////////////////////////////////////////////////////
// double Logistic_Distribution( double x )                                   //
//                                                                            //
//  Description:                                                              //
//     The Logistic distribution is the integral from -inf to x of the density//
//                       f(x) =  exp(-x) / (1 + exp(-x))^2                    //
//     i.e.                                                                   //
//                      Pr[X < x] = F(x) = 1 / (1 + exp(-x))                  //
//                                                                            //
//  Arguments:                                                                //
//     double x   Argument of the Logistic distribution Pr[X < x].            //
//                                                                            //
//  Return Values:                                                            //
//     A real number between 0 and 1.                                         //
//                                                                            //
//  Example:                                                                  //
//     double p, x;                                                           //
//                                                                            //
//     p = Logistic_Distribution(x);                                          //
////////////////////////////////////////////////////////////////////////////////

double Logistic_Distribution(double x)
{
   return 1.0 / (1.0 + exp(-x));
}
