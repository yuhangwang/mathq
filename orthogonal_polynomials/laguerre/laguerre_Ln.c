////////////////////////////////////////////////////////////////////////////////
// File: laguerre_Ln.c                                                        //
// Routine(s):                                                                //
//    Laguerre_Ln                                                             //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>             // required for fabsl()
#include <float.h>            // required for DBL_MAX

//                        Externally Defined Routines                         //

extern long double xLaguerre_Ln(long double x, int n);

////////////////////////////////////////////////////////////////////////////////
// double Laguerre_Ln(double x, int n)                                        //
//                                                                            //
//  Description:                                                              //
//     The Laguerre polynomials, Ln,  are orthogonal on the interval [0,inf)  //
//     with weight function w(x) = exp(-x) on [0, inf) and 0 elsewhere.  The  //
//     Laguerre polynomials are normalized so that Ln(0) = 1.                 //
//             <Ln,Lm> = 0    if n != m,                                      //
//             <Ln,Ln> = 1    if n >= 0.                                      //
//     This routine calculates Ln(x) using the function xLaguerre_Ln in the   //
//     file xlaguerre.c which in turn uses the recursion formula:             //
//       (k+1) L[k+1](x) = (2k+1-x) L[k](x) - k L[k-1](x), k = 1,...,n-1      //
//              L[0](x) = 1, L[1](x) = 1-x.                                   //
//                                                                            //
//  Arguments:                                                                //
//     double x                                                               //
//        The argument of the Laguerre polynomial Ln.                         //
//     int    n                                                               //
//        The degree of the Laguerre polynomial Ln.                           //
//                                                                            //
//  Return Value:                                                             //
//     Ln(x) if n is a nonnegative integer.  If n is negative, 0 is returned. //
//     If Ln(x) > DBL_MAX, then DBL_MAX is returned and if Ln(x) < -DBL_MAX   //
//     then -DBL_MAX is returned.                                             //
//                                                                            //
//  Example:                                                                  //
//     double Ln;                                                             //
//     double x;                                                              //
//     int    n;                                                              //
//                                                                            //
//     (user code to set x and n)                                             //
//                                                                            //
//     Ln = Laguerre_Ln(x, n);                                                //
////////////////////////////////////////////////////////////////////////////////
double Laguerre_Ln(double x, int n)
{
   long double Ln;

   if (n < 0) return 0.0;
   Ln = xLaguerre_Ln((long double)x, n);
   if (fabsl(Ln) < DBL_MAX) return (double) Ln;
   return (Ln > 0.0L) ? DBL_MAX : -DBL_MAX;

}
