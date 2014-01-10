////////////////////////////////////////////////////////////////////////////////
// File: laplace_distribution.c                                               //
// Routine(s):                                                                //
//    Laplace_Distribution                                                    //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                              // required for exp(), fabs()

////////////////////////////////////////////////////////////////////////////////
// double Laplace_Distribution( double x )                                    //
//                                                                            //
//  Description:                                                              //
//     Laplace's distribution is the integral from -inf to x of the density   //
//                           f(x) = (1/2) exp(-|x|)                           //
//     i.e.                                                                   //
//                      Pr[X < x] = F(x) = (1/2) exp(x)      if x < 0 and     //
//                                       = 1 - (1/2) exp(-x) if x > 0.        //
//                                                                            //
//  Arguments:                                                                //
//     double x   Argument of Laplace's distribution Pr[X < x].               //
//                                                                            //
//  Return Values:                                                            //
//     A real number between 0 and 1.                                         //
//                                                                            //
//  Example:                                                                  //
//     double p, x;                                                           //
//                                                                            //
//     p = Laplace_Distribution(x);                                           //
////////////////////////////////////////////////////////////////////////////////

double Laplace_Distribution(double x)
{
   double temp = 0.5 * exp(-fabs(x));

   return (x <= 0.0) ? temp : 1.0 - temp;
}
