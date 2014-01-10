////////////////////////////////////////////////////////////////////////////////
// File: hermite_Hn.c                                                         //
// Routine(s):                                                                //
//    Hermite_Hn                                                              //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>             // required for fabsl()                       
#include <float.h>            // required for DBL_MAX

//                        Externally Defined Routines                         //

extern long double xHermite_Hn(long double x, int n);

////////////////////////////////////////////////////////////////////////////////
// double Hermite_Hn(double x, int n)                                         //
//                                                                            //
//  Description:                                                              //
//     There are several systems of orthogonal polynomials called Hermite     //
//     polynomials.  The version programmed here is the system with weight    //
//     function w(x) = exp(-x^2).  Normalized so that the coefficient of the  //
//     leading term of Hn is 2^n.  The inner product is given by:             //
//             <Hn,Hm> = 0    if n != m,                                      //
//             <Hn,Hn> = n! 2^n sqrt(pi) if n >= 0.                           //
//     This routine calculates H[n](x) using the function xHermite_Hn in the  //
//     file xhermite_Hn.c which in turn uses the recursion formula:           //
//            H[k+1](x) = 2x H[k](x) - 2k H[k-1](x), k = 1,...,n-1            //
//            H[0](x) = 1, H[1](x) = 2x                                       //
//     to evaluate Hn at x.                                                   //
//                                                                            //
//  Arguments:                                                                //
//     double x                                                               //
//        The argument of the Hermite polynomial, Hn, of degree n.            //
//     int    n                                                               //
//        The degree of the Hermite polynomial.                               //
//                                                                            //
//  Return Value:                                                             //
//     Hn(x) if n is a nonnegative integer.  If n is negative, 0 is returned. //
//     If Hn(x) > DBL_MAX, then DBL_MAX is returned and if Hn(x) < -DBL_MAX   //
//     then -DBL_MAX is returned.                                             //
//                                                                            //
//  Example:                                                                  //
//     double Hn;                                                             //
//     double x;                                                              //
//     int    n;                                                              //
//                                                                            //
//     (user code to set x and n)                                             //
//                                                                            //
//     Hn = Hermite_Hn(x, n);                                                 //
////////////////////////////////////////////////////////////////////////////////
double Hermite_Hn(double x, int n)
{
   long double Hn;

   if (n < 0) return 0.0;
   Hn = xHermite_Hn((long double)x, n);
   if (fabsl(Hn) < DBL_MAX) return (double) Hn;
   return (Hn > 0.0L) ? DBL_MAX : -DBL_MAX;
}
