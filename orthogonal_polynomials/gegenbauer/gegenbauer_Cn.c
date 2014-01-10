////////////////////////////////////////////////////////////////////////////////
// File: gegenbauer_Cn.c                                                      //
// Routine(s):                                                                //
//    Gegenbauer_Cn                                                           //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>             // required for fabsl()
#include <float.h>            // required for DBL_MAX

//                        Externally Defined Routines                         //

extern long double xGegenbauer_Cn(long double x, long double alpha, int n);

////////////////////////////////////////////////////////////////////////////////
// double Gegenbauer_Cn(double x, double alpha, int n)                        //
//                                                                            //
//  Description:                                                              //
//     The Gegenbauer polynomials with parameter alpha > -1/2, C^(alpha)n(x)  //
//     are orthogonal on the interval [-1,1] with weight function             //
//     w(x) = (1-x^2)^(alpha-1/2) on [-1,1] and 0 elsewhere.  Further they are//
//     normalized so that C^(alpha)n(1) = gamma(n+2alpha)/(n! gamma(2*alpha)).//
//     Gegenbauer polynomials are also called the ultraspherical polynomials. //
//     For convenience in typing, the superscript (alpha) is dropped so that  //
//     C^(alpha)n(x) will be denoted by Cn(x), i.e. the parameter alpha is    //
//     to be understood.                                                      //
//      <Cn,Cm> = 0                                      if n != m,           //
//      <Cn,Cn> = 2^(1-2alpha) gamma(n+2alpha) pi                             //
//                / (n! (n+alpha)gamma^2(alpha) )        if n >= 0.           //
//     This routine calculates, Cn(x), the Gegenbauer polynomial of degree n  //
//     with parameter alpha evaluated at x  using the function xGegenbauer_Cn //
//     in the file xgegenbauer_Cn.c which in turn uses the recursion formula: //
//     (k+1) C[k+1](x) = 2(k + alpha) x C[k](x) - (k + 2 alpha - 1) C[k-1](x),//
//                                                         k = 1,...,n-1.     //
//     C[0](x) = 1, C[1](x) = 2*alpha*x.                                      //
//                                                                            //
//  Arguments:                                                                //
//     double x                                                               //
//        The argument of the Gegenbauer polynomial with parameter alpha of   //
//        degree n.                                                           //
//     double alpha                                                           //
//        The parameter of the Gegenbauer polynomial.                         //
//     int    n                                                               //
//        The degree of the Gegenbauer polynomial.                            //
//                                                                            //
//  Return Value:                                                             //
//     Cn(x) if n is a nonnegative integer.  If n is negative, 0 is returned. //
//     If Cn(x) > DBL_MAX, then DBL_MAX is returned and if Cn(x) < -DBL_MAX   //
//     then -DBL_MAX is returned (this applies only if |x| > 1).              //
//                                                                            //
//  Example:                                                                  //
//     double x, alpha, Cn;                                                   //
//     int    n;                                                              //
//                                                                            //
//     (user code to set x, alpha, and n)                                     //
//                                                                            //
//     Cn = Gegenbauer_Cn(x, alpha, n);                                       //
////////////////////////////////////////////////////////////////////////////////
double Gegenbauer_Cn(double x, double alpha, int n)
{
   long double Cn;

   if (n < 0) return 0.0;
   Cn = xGegenbauer_Cn((long double)x, (long double) alpha, n);
   if (fabsl(Cn) < DBL_MAX) return (double) Cn;
   return (Cn > 0.0L) ? DBL_MAX : -DBL_MAX;
}
