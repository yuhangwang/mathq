////////////////////////////////////////////////////////////////////////////////
// File: xchebyshev_shifted_Un.c                                              //
// Routine(s):                                                                //
//    xChebyshev_Shifted_Un                                                   //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>        // required for fabsl(), acosl(), cosl(), & sinl()
#define N 6              // cross-over between recursion and explicit form

////////////////////////////////////////////////////////////////////////////////
// long double xChebyshev_Shifted_Un(long double x, int n)                    //
//                                                                            //
//  Description:                                                              //
//     The shifted Chebyshev polynomials of the second kind, Un*(x)=Un(2x-1), //
//     are orthogonal on the interval [0,1] with weight function              //
//     w(x) = sqrt(x-x^2) on [0,1] and 0 elsewhere.  Further Un* is normalized//
//     so that U*n(1) = n + 1.                                                //
//             <Un*,Um*> = 0    if n != m,                                    //
//             <Un*,Un*> = pi/8 if n >= 0,                                    //
//     This routine calculates U*[n](x) using the following techniques:       //
//     If n < N or if |2x-1| > 1, where N is defined below the use recursion  //
//              U*[k+1](x) = 2(2x - 1) U*[k](x) - U*[k-1](x), k = 1,...,n-1   //
//              U*[0](x) = 1, U*[1](x) = 4x - 2.                              //
//     otherwise use the explicit formula                                     //
//        U*[n](x) = sin((n+1)*theta)/sin(theta), where theta = acos(2x-1).   //
//                                                                            //
//  Arguments:                                                                //
//     long double x                                                          //
//        The argument of the shifted Chebyshev polynomial Un*.               //
//     int    n                                                               //
//        The degree of the shifted Chebyshev polynomial.                     //
//                                                                            //
//  Return Value:                                                             //
//     Un*(x) if n is a nonnegative integer.  If n is negative, 0 is returned.//
//                                                                            //
//  Example:                                                                  //
//     long double x, Un;                                                     //
//     int    n;                                                              //
//                                                                            //
//     (user code to set x and n)                                             //
//                                                                            //
//     Un = xChebyshev_shifted_Un(x, n);                                      //
////////////////////////////////////////////////////////////////////////////////
long double xChebyshev_Shifted_Un(long double x, int n)
{
   long double two_x_m1 = x + x - 1.0L;
   long double four_x_m2 = two_x_m1 + two_x_m1;
   long double theta;
   long double sin_theta;
   long double U0, U1, Un;
   long double xn1 = (long double)(n + 1);
   int k;

   if (n < 0) return 0.0L;
   if ( x == 1.0L ) return xn1;
   if ( x == 0.0L ) return (n % 2 == 0) ? xn1 : -xn1;
    
   if ( (n > N) && (fabsl(two_x_m1) < 1.0L) ) {
      theta = acosl(two_x_m1);
      sin_theta = sinl(theta);
      if (sin_theta != 0.0L)
         return sinl(xn1 * theta) / sin_theta;
      else 
         return xn1 * cosl(xn1 * theta) / two_x_m1;
   }
      
                    // Initialize the recursion process. //
 
   if (n == 0) return 1.0L;
   if (n == 1) return four_x_m2;
   U0 = 1.0L;
   U1 = four_x_m2;

                             // Calculate U*nx) //

   for (k = 2; k <= n; k++) {
      Un = four_x_m2 * U1 - U0;
      U0 = U1;
      U1 = Un;
   }

   return Un;
}
