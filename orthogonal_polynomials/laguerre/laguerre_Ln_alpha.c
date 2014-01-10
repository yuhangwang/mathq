////////////////////////////////////////////////////////////////////////////////
// File: laguerre_Ln_alpha.c                                                  //
// Routine(s):                                                                //
//    Laguerre_Ln_alpha                                                       //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>             // required for fabsl()
#include <float.h>            // required for DBL_MAX

//                        Externally Defined Routines                         //

long double xLaguerre_Ln_alpha(long double x, long double alpha, int n);

////////////////////////////////////////////////////////////////////////////////
// double Laguerre_Ln_alpha(double x, double alpha, int n)                    //
//                                                                            //
//  Description:                                                              //
//     The generalized Laguerre polynomials with parameter alpha > -1,        //
//     L^(alpha)n(x), are orthogonal on the interval [0,inf) with weight      //
//     function w(x) = x^alpha exp(-x) on [0,inf) and 0 elsewhere.            //
//     For convenience in typing, the superscript (alpha) is dropped so that  //
//     L^(alpha)[n](x) will be denoted by L[n](x), i.e. the alpha is to be    //
//     understood.  The generalized Laguerre polynomials are normalized so    //
//     that Ln(0) = gamma(alpha + n + 1) / (gamma(alpha+1) n!).               //
//             <Ln,Lm> = 0                         if n != m,                 //
//             <Ln,Ln> = gamma(alpha + n + 1) / n! if n >= 0.                 //
//     This routine calculates, Ln(x), the generalized Laguerre polynomial    //
//     with parameter alpha of degree n evaluated at x using the function     //
//     xLaguerre_Ln_alpha in the file xlaguerre_Ln_alpha.c which in turn uses //
//     the recursion formula:                                                 //
//       (k+1) L[k+1](x) = (2k+alpha+1-x) L[k](x)                             //
//                                  - (k+alpha) L[k-1](x), k = 1,...,n-1      //
//              L[0](x) = 1, L[1](x) = alpha + 1 - x.                         //
//                                                                            //
//  Arguments:                                                                //
//     double x                                                               //
//        The argument of the generalized Laguerre polynomial with parameter  //
//        alpha.                                                              //
//     double alpha                                                           //
//        The parameter of the generalized Laguerre polynomial.               //
//     int    n                                                               //
//        The degree of the generalized Laguerre polynomial with parameter    //
//        alpha.                                                              //
//                                                                            //
//  Return Value:                                                             //
//     Ln(x) if n is a nonnegative integer.  If n is negative, 0 is returned. //
//     If Ln(x) > DBL_MAX, then DBL_MAX is returned and if Ln(x) < -DBL_MAX   //
//     then -DBL_MAX is returned.                                             //
//                                                                            //
//  Example:                                                                  //
//     double x, alpha, Ln;                                                   //
//     int    n;                                                              //
//                                                                            //
//     (user code to set x, alpha, and n)                                     //
//                                                                            //
//     Ln = Laguerre_Ln_alpha(x, alpha, n);                                   //
////////////////////////////////////////////////////////////////////////////////
double Laguerre_Ln_alpha(double x, double alpha, int n)
{
   long double Ln;

   if (n < 0) return 0.0;
   Ln = xLaguerre_Ln_alpha((long double) x, (long double) alpha, n);
   if (fabsl(Ln) < DBL_MAX) return (double) Ln;
   return (Ln > 0.0L) ? DBL_MAX : -DBL_MAX;
}
