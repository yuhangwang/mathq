////////////////////////////////////////////////////////////////////////////////
// File: hermite_Hen.c                                                        //
// Routine(s):                                                                //
//    Hermite_Hen                                                             //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>             // required for fabsl()                       
#include <float.h>            // required for DBL_MAX

//                        Externally Defined Routines                         //

extern long double xHermite_Hen(long double x, int n);

////////////////////////////////////////////////////////////////////////////////
// double Hermite_Hen(double x, int n)                                        //
//                                                                            //
//  Description:                                                              //
//     There are several systems of orthogonal polynomials called Hermite     //
//     polynomials.  The version programmed here is the system with weight    //
//     function w(x) = exp(-x^2 / 2), denoted by He_n, and normalized so that //
//     the coefficient of the leading term of He_n is 1.  The inner product is//
//     given by:                                                              //
//             <He_n,He_m> = 0             if n != m,                         //
//             <He_n,He_n> = sqrt(2pi) n!  if n >= 0.                         //
//     This routine calculates He[n](x) using the function xHermite_Hen in the//
//     file xhermite_Hen.c which in turn uses the recursion formula:          //
//            He[k+1](x) = x He[k](x) - k He[k-1](x), k = 1,...,n-1           //
//            He[0](x) = 1, He[1](x) = x                                      //
//     to evaluate He_n at x.                                                 //
//                                                                            //
//  Arguments:                                                                //
//     double x                                                               //
//        The argument of the Hermite polynomial, He_n, of degree n.          //
//     int    n                                                               //
//        The degree of the Hermite polynomial.                               //
//                                                                            //
//  Return Value:                                                             //
//     He_n(x) if n is a nonnegative integer.  If n is negative, 0 is         //
//     returned. If He_n(x) > DBL_MAX, then DBL_MAX is returned and if        //
//     He_n(x) < -DBL_MAX then -DBL_MAX is returned.                          //
//                                                                            //
//  Example:                                                                  //
//     double Hen;                                                            //
//     double x;                                                              //
//     int    n;                                                              //
//                                                                            //
//     (user code to set x and n)                                             //
//                                                                            //
//     Hen = Hermite_Hen(x, n);                                               //
////////////////////////////////////////////////////////////////////////////////
double Hermite_Hen(double x, int n)
{
   long double Hen;

   if (n < 0) return 0.0;
   Hen = xHermite_Hen((long double)x, n);
   if (fabsl(Hen) < DBL_MAX) return (double) Hen;
   return (Hen > 0.0L) ? DBL_MAX : -DBL_MAX;
}
