////////////////////////////////////////////////////////////////////////////////
// File: chebyshev_shifted_Un.c                                               //
// Routine(s):                                                                //
//    Chebyshev_Shifted_Un                                                    //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>        // required for fabsl()
#include <float.h>       // required for DBL_MAX

//                        Externally Defined Routines                         //
long double xChebyshev_Shifted_Un(long double x, int n);

////////////////////////////////////////////////////////////////////////////////
// double Chebyshev_Shifted_Un(double x, int n)                               //
//                                                                            //
//  Description:                                                              //
//     The shifted Chebyshev polynomials of the second kind, Un*(x)=Un(2x-1), //
//     are orthogonal on the interval [0,1] with weight function              //
//     w(x) = sqrt(x-x^2) on [0,1] and 0 elsewhere.  Further Un* is normalized//
//     so that U*n(1) = n + 1.                                                //
//             <Un*,Um*> = 0    if n != m,                                    //
//             <Un*,Un*> = pi/8 if n >= 0,                                    //
//     This routine uses xChebyhev_Shift_Un to calculate Un*(x) which in turn //
//     uses the following techniques:                                         //
//     If n < N or if |2x-1| > 1, where N is defined in xchebyshev_shift_Un.c,//
//     then use the recursion                                                 //
//              U*[k+1](x) = 2(2x - 1) U*[k](x) - U*[k-1](x), k = 1,...,n-1   //
//              U*[0](x) = 1, U*[1](x) = 4x - 2.                              //
//     otherwise use the explicit formula                                     //
//        U*[n](x) = sin((n+1)*theta)/sin(theta), where theta = acos(2x-1).   //
//                                                                            //
//  Arguments:                                                                //
//     double x                                                               //
//        The argument of the shifted Chebyshev polynomial Un*.               //
//     int    n                                                               //
//        The degree of the shifted Chebyshev polynomial.                     //
//                                                                            //
//  Return Value:                                                             //
//     Un*(x) if n is a nonnegative integer.  If n is negative, 0 is returned.//
//     If Un*(x) > DBL_MAX, then DBL_MAX is returned and if Un*(x) < -DBL_MAX //
//     then -DBL_MAX is returned (this applies only if x < 0 or x > 1.)       //
//                                                                            //
//  Example:                                                                  //
//     double x, Un;                                                          //
//     int    n;                                                              //
//                                                                            //
//     (user code to set x and n)                                             //
//                                                                            //
//     Un = Chebyshev_shifted_Un(x, n);                                       //
////////////////////////////////////////////////////////////////////////////////
double Chebyshev_Shifted_Un(double x, int n)
{
   long double Un;

   if (n < 0) return 0.0;
   Un = xChebyshev_Shifted_Un((long double)x, n);
   if (fabsl(Un) < DBL_MAX) return (double) Un;
   return (Un > 0.0L) ? DBL_MAX : -DBL_MAX;

}
