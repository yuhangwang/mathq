////////////////////////////////////////////////////////////////////////////////
// File: xchebyshev_shifted_Wn.c                                              //
// Routine(s):                                                                //
//    xChebyshev_Shifted_Wn                                                   //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>        // required for fabsl(), acosl(), sinl(), & cosl()
#define N 4              // cross-over between recursion and explicit form
////////////////////////////////////////////////////////////////////////////////
// long double xChebyshev_Shifted_Wn(long double x, int n)                    //
//                                                                            //
//  Description:                                                              //
//     The shifted Chebyshev polynomials of the fourth kind, Wn*(x)=Wn(2x-1), //
//     are orthogonal on the interval [0,1] with weight function              //
//     w(x) = sqrt((1-x)/x) on [0,1] and 0 elsewhere.  Further Wn*(x) is      //
//     normalized so that W*n(1) = 2n+1.                                      //
//             <Wn*,Wm*> = 0    if n != m,                                    //
//             <Wn*,Wn*> = pi/2 if n >= 0,                                    //
//     This routine calculates Wn*(x) using the following techniques:         //
//     If n < N or if x > 1 or x < 0, where N is defined above then recurse   //
//              W*[k+1](x) = 2(2x - 1) W*[k](x) - W*[k-1](x), k = 1,...,n-1   //
//              W*[0](x) = 1, W*[1](x) = 4x - 1.                              //
//     otherwise use the explicit formula                                     //
//      W*[n](x) = sin((n+1/2)*theta)/sin(theta/2), where theta = acos(2x-1). //
//                                                                            //
//  Arguments:                                                                //
//     long double x                                                          //
//        The argument of the shifted Chebyshev polynomial Wn*.               //
//     int    n                                                               //
//        The degree of the shifted Chebyshev polynomial.                     //
//                                                                            //
//  Return Value:                                                             //
//     Wn*(x) if n is a nonnegative integer.  If n is negative, 0 is returned.//
//                                                                            //
//  Example:                                                                  //
//     long double x, Wn;                                                     //
//     int    n;                                                              //
//                                                                            //
//     (user code to set x and n)                                             //
//                                                                            //
//     Wn = xChebyshev_Shifted_Wn(x, n);                                      //
////////////////////////////////////////////////////////////////////////////////
long double xChebyshev_Shifted_Wn(long double x, int n)
{
   long double two_x_m1 = x + x - 1.0L;
   long double four_x_m2 = two_x_m1 + two_x_m1;
   long double theta;
   long double sin_theta2;
   long double cos_theta2;
   long double W0, W1, Wn;
   int k;

   if (n < 0) return 0.0L;
   if ( x == 1.0L ) return (long double)(n+n+1);
   if ( x == 0.0L ) return (n % 2 == 0) ? 1.0L : -1.0L;
    
   if ( (n > N) && (fabsl(two_x_m1) < 1.0L) ) {
      theta = acosl(two_x_m1);
      sin_theta2 = sinl(theta/2.0L);
      cos_theta2 = cosl(theta/2.0L);
      if (cos_theta2 != 1.0L)
         return sinl(((long double)(n) + 0.5L) * theta) / sin_theta2;
      else 
         return (n % 2 == 0) ? 1.0L : -1.0L;
   }
      
                    // Initialize the recursion process. //
 
   if (n == 0) return 1.0L;
   if (n == 1) return four_x_m2 + 1.0L;
   W0 = 1.0L;
   W1 = four_x_m2 + 1.0L;

                             // Calculate V*nx) //

   for (k = 2; k <= n; k++) {
      Wn = four_x_m2 * W1 - W0;
      W0 = W1;
      W1 = Wn;
   }

   return Wn;
}
