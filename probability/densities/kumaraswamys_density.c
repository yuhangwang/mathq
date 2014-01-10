////////////////////////////////////////////////////////////////////////////////
// File: kumaraswamys_density.c                                               //
// Routine(s):                                                                //
//    Kumaraswamys_Density                                                    //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                                   // required for pow()
#include <float.h>                                  // required for DBL_MAX

////////////////////////////////////////////////////////////////////////////////
// double Kumaraswamys_Density( double x, double a, double b)                 //
//                                                                            //
//  Description:                                                              //
//     The probability density function of Kumaraswamys distribution is       //
//                   f(x) = a b x^(a-1) (1 - x^a)^(b-1)   for 0 < x < 1,      //
//                        = 0                             elsewhere.          //
//     The shape parameters, a and b, must be positive.                       //
//     If 0 < a < 1, the probability density function has a discontinuity at  //
//     x = 0 and if 0 < b < 1, it has a discontinuity at x = 1.               //
//                                                                            //
//  Arguments:                                                                //
//     double x   Argument of Kumaraswamys density.                           //
//     double a   Shape parameter of Kumarasways distribution, the exponent of//
//                x, a > 0.                                                   //
//     double b   Shape parameter of Kumarasways distribution, the exponent of//
//                (1-x^a), b > 0.                                             //
//                                                                            //
//  Return Values:                                                            //
//     A real number between 0 and 1/2.                                       //
//                                                                            //
//  Example:                                                                  //
//     double a, b, p, x;                                                     //
//                                                                            //
//     p = Kumaraswamys_Density(x, a, b);                                     //
////////////////////////////////////////////////////////////////////////////////

double Kumaraswamys_Density(double x, double a, double b)
{
   double temp;

   if (x <= 0.0 || x >= 1.0) return 0.0;

   temp = pow(x, a-1.0);
   return a * b * temp * pow(1.0 - temp * x, b - 1.0);
}
