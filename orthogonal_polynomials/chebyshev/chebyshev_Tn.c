////////////////////////////////////////////////////////////////////////////////
// File: chebyshev_Tn.c                                                       //
// Routine(s):                                                                //
//    Chebyshev_Tn                                                            //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>             // required for fabsl()                           
#include <float.h>            // required for DBL_MAX

//                        Externally Defined Routines                         //

extern long double xChebyshev_Tn(long double x, int n);

////////////////////////////////////////////////////////////////////////////////
// double Chebyshev_Tn(double x, int n)                                       //
//                                                                            //
//  Description:                                                              //
//     The Chebyshev polynomials of the first kind, Tn(x), are orthogonal on  //
//     the interval [-1,1] with weight function w(x) = 1 / sqrt(1-x^2) on the //
//     interval (-1,1) and 0 elsewhere and are normalized so that Tn(1) = 1.  //
//             <Tn,Tm> = 0    if n != m,                                      //
//             <Tn,Tn> = pi/2 if n > 0,                                       //
//             <T0,T0> = pi.                                                  //
//     This routine uses xChebyhev_Tn to calculate Tn(x) which in turn uses   //
//     the following techniques:                                              //
//     If n < N or if |x| > 1, where N is defined in xchebyshev_Tn.c, then    //
//     use the recursion                                                      //
//              T[k+1](x) = 2x T[k](x) - T[k-1](x), k = 1,2,...,n-1           //
//              T[0](x) = 1, T[1](x) = x                                      //
//     otherwise use the explicit formula T[n](x) = cos(n*acos(x)).           //
//     This routine calculates T[n](x), the Chebyshev polynomial of the       //
//     first kind of degree n evaluated at x where x and the return value are //
//     of type double.                                                        //
//                                                                            //
//  Arguments:                                                                //
//     double x                                                               //
//        The argument of the Chebyshev polynomial Tn.                        //
//     int    n                                                               //
//        The degree of the Chebyshev polynomial.                             //
//                                                                            //
//  Return Value:                                                             //
//     Tn(x) if n is a nonnegative integer.  If n is negative, 0 is returned. //
//     If Tn(x) > DBL_MAX, then DBL_MAX is returned and if Tn(x) < -DBL_MAX   //
//     then -DBL_MAX is returned (this applies only if |x| > 1).              //
//                                                                            //
//  Example:                                                                  //
//     double x, Tn;                                                          //
//     int    n;                                                              //
//                                                                            //
//     (user code to set x and n)                                             //
//                                                                            //
//     Tn = Chebyshev_Tn(x, n);                                               //
////////////////////////////////////////////////////////////////////////////////
double Chebyshev_Tn(double x, int n)
{
   long double Tn;

   if (n < 0) return 0.0;
   Tn = xChebyshev_Tn((long double)x, n);
   if (fabsl(Tn) < DBL_MAX) return (double) Tn;
   return (Tn > 0.0L) ? DBL_MAX : -DBL_MAX;
}
