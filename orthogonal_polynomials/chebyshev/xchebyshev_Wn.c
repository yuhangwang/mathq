////////////////////////////////////////////////////////////////////////////////
// File: xchebyshev_Wn.c                                                      //
// Routine(s):                                                                //
//    xChebyshev_Wn                                                           //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>             // required for fabsl(), acosl(), sinl() and 
                              //   cosl()
#define N 8                   // cross-over between recursion and using       
                              //   Wn = sin((n+1/2)theta)/sin(theta/2),
                              //   where x = cos(theta)
////////////////////////////////////////////////////////////////////////////////
// long double xChebyshev_Wn(long double x, int n)                            //
//                                                                            //
//  Description:                                                              //
//     The Chebyshev polynomials of the fourth kind, Wn(x), are orthogonal on //
//     the interval [-1,1] with weight function w(x) = sqrt((1-x)/(1+x)) on   //
//     [-1,1] and 0 elsewhere and are normalized so that Wn(1) = 2n+1.        //
//             <Wn,Wm> = 0    if n != m,                                      //
//             <Wn,Wn> = pi   if n >= 0.                                      //
//     This routine calculates W[n](x) using the following techniques:        //
//     If n < N or if |x| > 1, where N is defined above then use the recursion//
//              W[k+1](x) = 2x W[k](x) - W[k-1](x), k = 1,2,...,n-1           //
//              W[0](x) = 1, W[1](x) = 2x+1.                                  //
//     otherwise use W[n](x) = sin((n+1/2)*theta)/sin(theta/2), where         //
//     theta = acos(x).                                                       //
//                                                                            //
//  Arguments:                                                                //
//     long double x                                                          //
//        The argument of the Chebyshev polynomial Wn.                        //
//     int    n                                                               //
//        The degree of the Chebyshev polynomial.                             //
//                                                                            //
//  Return Value:                                                             //
//     Wn(x) if n is a nonnegative integer.  If n is negative, 0 is returned. //
//                                                                            //
//  Example:                                                                  //
//     long double x, Wn;                                                     //
//     int    n;                                                              //
//                                                                            //
//     (user code to set x and n)                                             //
//                                                                            //
//     Wn = xChebyshev_Wn(x, n);                                              //
////////////////////////////////////////////////////////////////////////////////
long double xChebyshev_Wn(long double x, int n)
{
   long double two_x = x + x;
   long double theta;
   long double cos_theta2;
   long double sin_theta2;
   long double W0, W1, Wn;
   int k;

   if (n < 0) return 0.0L;
   if ( fabsl(x) == 1.0L ) { 
      if (x > 0.0L) return (long double)(n+n+1);
      else
         return (n % 2 == 0) ? 1.0L : -1.0L;
   }
    
   if ( (n > N) && (fabsl(x) < 1.0L) ) {
      theta = acosl(x);
      cos_theta2 = cosl(theta/2.0L);
      sin_theta2 = sinl(theta/2.0L);
      if (cos_theta2 != 1.0L)
         return sinl(((long double)(n) + 0.5L) * theta) / sin_theta2;
      else 
         if ( x > 0.0L) return (long double) (n+n+1);
         else return (n % 2 == 0) ? 1.0L : -1.0L;
   }
      
                    // Initialize the recursion process. //
 
   if (n == 0) return 1.0L;
   if (n == 1) return two_x + 1.0L;
   W0 = 1.0L;
   W1 = two_x + 1.0L;

                             // Calculate Wn(x) //

   for (k = 2; k <= n; k++) {
      Wn = two_x * W1 - W0;
      W0 = W1;
      W1 = Wn;
   }

   return Wn;
}
