////////////////////////////////////////////////////////////////////////////////
// File: xchebyshev_Un.c                                                      //
// Routine(s):                                                                //
//    xChebyshev_Un                                                           //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>             // required for fabsl(), acosl(), sinl() and
                              //   cosl()
#define N 8                   // cross-over between recursion and using       
                              //   Un = sin((n+1)theta)/sin(theta),
                              //   where x = cos(theta)

////////////////////////////////////////////////////////////////////////////////
// long double xChebyshev_Un(long double x, int n)                            //
//                                                                            //
//  Description:                                                              //
//     The Chebyshev polynomials of the second kind, Un(x), are orthogonal on //
//     the interval [-1,1] with weight function w(x) = sqrt(1-x^2) on the     //
//     interval [-1,1] and 0 elsewhere and are normalized so that Un(1) = n+1.//
//             <Un,Um> = 0    if n != m,                                      //
//             <Un,Un> = pi/2 if n >= 0.                                      //
//     This routine calculates Un(x) using the following techniques:          //
//     If n < N or if |x| > 1, where N is defined below then use the recursion//
//              U[k+1](x) = 2x U[k](x) - U[k-1](x), k = 1,2,...,n-1           //
//              U[0](x) = 1, U[1](x) = 2x.                                    //
//     otherwise use U[n](x) = sin((n+1)*theta)/sin(theta), where             //
//     theta = acos(x).                                                       //
//                                                                            //
//  Arguments:                                                                //
//     long double x                                                          //
//        The argument of the Chebyshev polynomial Un.                        //
//     int    n                                                               //
//        The degree of the Chebyshev polynomial.                             //
//                                                                            //
//  Return Value:                                                             //
//     Un(x) if n is a nonnegative integer.  If n is negative, 0 is returned. //
//                                                                            //
//  Example:                                                                  //
//     long double Un;                                                        //
//     long double x;                                                         //
//     int    n;                                                              //
//                                                                            //
//     (user code to set x and n)                                             //
//                                                                            //
//     Un = xChebyshev_Un(x, n);                                              //
////////////////////////////////////////////////////////////////////////////////
long double xChebyshev_Un(long double x, int n)
{
   long double two_x = x + x;
   long double theta;
   long double sin_theta;
   long double U0, U1, Un;
   int k;

   if (n < 0) return 0.0L;
   if ( fabsl(x) == 1.0L ) { 
      if (x > 0.0L) return (long double)(n + 1);
      else if (n % 2 == 0) return (long double)(n + 1);
      else return -(long double)(n + 1);
   }
    
   if ( (n > N) && (fabsl(x) < 1.0L) ) {
      theta = acosl(x);
      sin_theta = sinl(theta);
      if (sin_theta != 0.0L)
         return sinl((long double)(n + 1) * theta) / sin_theta;
      else 
         return (long double) (n + 1) * cosl((long double)(n + 1) * theta) / x;
   }
      
                    // Initialize the recursion process. //
 
   if (n == 0) return 1.0L;
   if (n == 1) return two_x;
   U0 = 1.0L;
   U1 = two_x;

                             // Calculate Un(x) //

   for (k = 2; k <= n; k++) {
      Un = two_x * U1 - U0;
      U0 = U1;
      U1 = Un;
   }

   return Un;
}
