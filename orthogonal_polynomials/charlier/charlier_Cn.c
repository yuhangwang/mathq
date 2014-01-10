////////////////////////////////////////////////////////////////////////////////
// File: charlier_Cn.c                                                        //
// Routine(s):                                                                //
//    Charlier_Cn                                                             //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>             // required for fabsl()
#include <float.h>            // required for DBL_MAX

//                        Externally Defined Routines                         //

long double xCharlier_Cn(long double x, long double a, int n);

////////////////////////////////////////////////////////////////////////////////
// double Charlier_Cn(double x, double a, int n)                              //
//                                                                            //
//  Description:                                                              //
//     The Charlier polynomials are orthogonal on the interval [0,inf)        //
//     with Poisson weight function w(x) = Sum [( (exp(-a) a^j / j! ) D(x-j)],//
//     where a > 0 is the mean of the Poisson distribution, D() is the Dirac  //
//     delta function, and the sum is over j = 0, 1, ... .                    //
//     This routine calculates, C[n](x), the n-th Charlier polynomial         //
//     evaluated at x using the function xCharlier_Cn in the file             //
//     xcharlier.c which in turn uses following recursion formula:            //
//               a C[k+1](x) = (k + a - x) C[k](x) - k C[k-1]                 //
//               C[0](x) = 1, C[1](x) = 1 - x/a.                              //
//                                                                            //
//  Arguments:                                                                //
//     double x                                                               //
//        The argument of the Charlier polynomial with parameter a.           //
//     double a                                                               //
//        The parameter of the Charlier polynomial, a > 0.                    //
//     int    n                                                               //
//        The degree of the Charlier polynomial with parameter a.             //
//                                                                            //
//  Return Value:                                                             //
//     Cn(x) if n is a nonnegative integer.  If n is negative, 0 is returned. //
//     If Cn(x) > DBL_MAX, then DBL_MAX is returned and if Cn(x) < -DBL_MAX   //
//     then -DBL_MAX is returned.                                             //
//                                                                            //
//  Example:                                                                  //
//     double Cn;                                                             //
//     double x;                                                              //
//     double a;                                                              //
//     int    n;                                                              //
//                                                                            //
//     (user code to set x, a, and n )                                        //
//                                                                            //
//     Cn = Charlier_Cn(x, a, n);                                             //
////////////////////////////////////////////////////////////////////////////////
double Charlier_Cn(double x, double a, int n)
{
   long double Cn;

   if (n < 0) return 0.0;

   Cn = xCharlier_Cn((long double) x, (long double) a, n);

   if (fabsl(Cn) < DBL_MAX) return (double) Cn;
   return (Cn > 0.0L) ? DBL_MAX : -DBL_MAX;

}
