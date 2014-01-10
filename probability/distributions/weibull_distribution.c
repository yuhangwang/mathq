////////////////////////////////////////////////////////////////////////////////
// File: weibull_distribution.c                                               //
// Routine(s):                                                                //
//    Weibull_Distribution                                                    //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                    // required for pow() and exp()

////////////////////////////////////////////////////////////////////////////////
// double Weibull_Distribution( double x, double a )                          //
//                                                                            //
//  Description:                                                              //
//     The Weibull distribution is the integral from -inf to x of the density //
//                                           0                if x <= 0,      //
//                              f(x) =  a x^(a-1) exp(-x^a)   if x > 0        //
//     i.e.                                                                   //
//                         Pr[X < x] = F(x) = 0               if x <= 0.      //
//                         Pr[X < x] = F(x) = 1 - exp(-x^a)   if x > 0.       //
//                                                                            //
//  Arguments:                                                                //
//     double x   Argument of the Weibull distribution Pr[X < x].             //
//     double a   Shape parameter of the Weibull distribution, a > 0.         //
//                                                                            //
//  Return Values:                                                            //
//     A real number between 0 and 1.                                         //
//                                                                            //
//  Example:                                                                  //
//     double a, p, x;                                                        //
//                                                                            //
//     p = Weibull_Distribution(x,a);                                         //
////////////////////////////////////////////////////////////////////////////////

double Weibull_Distribution(double x, double a)
{
   if (x <= 0.0) return 0.0;
   return 1.0 - exp(-pow(x,a));
}
