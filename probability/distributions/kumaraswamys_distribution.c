////////////////////////////////////////////////////////////////////////////////
// File: kumaraswamys_distribution.c                                          //
// Routine(s):                                                                //
//    Kumaraswamys_Distribution                                               //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                                   // required for pow()

////////////////////////////////////////////////////////////////////////////////
// double Kumaraswamys_Distribution( double x, double a, double b )           //
//                                                                            //
//  Description:                                                              //
//     Kumaraswamys distribution is the integral from -inf to x of the density//
//            f(x) = a b x^(a-1) (1 - x^a)^(b-1) for 0 < x < 1,               //
//                 = 0                           elsewhere                    //
//     i.e.                                                                   //
//            Pr[X < x] = F(x) = 0                 for x <= 0,                //
//                             = 1 - (1 - x^a)^b   for 0 < x < 1,             //
//                             = 1                 for x >= 1.                //
//     The shape parameters, a and b, must be positive.                       //
//     If 0 < a < 1, the probability density function has a discontinuity at  //
//     x = 0 and if 0 < b < 1, it has a discontinuity at x = 1.               //
//                                                                            //
//  Arguments:                                                                //
//     double x   Argument of Kumaraswamys distribution Pr[X < x].            //
//     double a   Shape parameter of Kumarasways distribution, the exponent of//
//                x, a > 0.                                                   //
//     double b   Shape parameter of Kumarasways distribution, the exponent of//
//                (1-x^a), b > 0.                                             //
//                                                                            //
//  Return Values:                                                            //
//     A real number between 0 and 1.                                         //
//                                                                            //
//  Example:                                                                  //
//     double a, b, p, x;                                                     //
//                                                                            //
//     p = Kumaraswamys_Distribution(x, a, b);                                //
////////////////////////////////////////////////////////////////////////////////

double Kumaraswamys_Distribution( double x, double a, double b )
{
   if ( x <= 0.0 ) return 0.0;
   if ( x >= 1.0 ) return 1.0;
   return 1.0 - pow(1.0 - pow(x,a), b);
}
