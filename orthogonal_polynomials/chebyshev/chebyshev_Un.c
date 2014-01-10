////////////////////////////////////////////////////////////////////////////////
// File: chebyshev_Un.c                                                       //
// Routine(s):                                                                //
//    Chebyshev_Un                                                            //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>             // required for fabsl()
#include <float.h>            // required for DBL_MAX

//                        Externally Defined Routines                         //

extern long double xChebyshev_Un(long double x, int n);

////////////////////////////////////////////////////////////////////////////////
// double Chebyshev_Un(double x, int n)                                       //
//                                                                            //
//  Description:                                                              //
//     The Chebyshev polynomials of the second kind, Un(x), are orthogonal on //
//     the interval [-1,1] with weight function w(x) = sqrt(1-x^2) on the     //
//     interval [-1,1] and 0 elsewhere and are normalized so that Un(1) = n+1.//
//             <Un,Um> = 0    if n != m,                                      //
//             <Un,Un> = pi/2 if n >= 0.                                      //
//     This routine calculates Un(x) using the following techniques:          //
//     If n < N or if |x| > 1, where N is defined in xchebyshev_Un.c, then    //
//     use the recursion                                                      //
//              U[k+1](x) = 2x U[k](x) - U[k-1](x), k = 1,2,...,n-1           //
//              U[0](x) = 1, U[1](x) = 2x.                                    //
//     otherwise use U[n](x) = sin((n+1)*theta)/sin(theta), where             //
//     theta = acos(x).                                                       //
//                                                                            //
//  Arguments:                                                                //
//     double x                                                               //
//        The argument of the Chebyshev polynomial Un.                        //
//     int    n                                                               //
//        The degree of the Chebyshev polynomial.                             //
//                                                                            //
//  Return Value:                                                             //
//     Un(x) if n is a nonnegative integer.  If n is negative, 0 is returned. //
//     If Un(x) > DBL_MAX, then DBL_MAX is returned and if Un(x) < -DBL_MAX   //
//     then -DBL_MAX is returned (this applies only if |x| > 1).              //
//                                                                            //
//  Example:                                                                  //
//     double x, Un;                                                          //
//     int    n;                                                              //
//                                                                            //
//     (user code to set x and n)                                             //
//                                                                            //
//     Un = Chebyshev_Un(x, n);                                               //
////////////////////////////////////////////////////////////////////////////////
double Chebyshev_Un(double x, int n)
{
   long double Un;

   if (n < 0) return 0.0;
   Un = xChebyshev_Un((long double)x, n);
   if (fabsl(Un) < DBL_MAX) return (double) Un;
   return (Un > 0.0L) ? DBL_MAX : -DBL_MAX;
}
