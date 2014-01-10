////////////////////////////////////////////////////////////////////////////////
// File: xchebyshev_shifted_Vn.c                                              //
// Routine(s):                                                                //
//    xChebyshev_Shifted_Vn                                                   //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>         // required for fabsl(), acosl(), cosl(), & sinl()
#define N 6               // cross-over between recursion and explicit form
////////////////////////////////////////////////////////////////////////////////
// long double xChebyshev_Shifted_Vn(long double x, int n)                    //
//                                                                            //
//  Description:                                                              //
//     The shifted Chebyshev polynomials of the third kind, Vn*(x) = Vn(2x-1),//
//     are orthogonal on the interval [0,1] with weight function              //
//     w(x) = sqrt(x/(1-x)) on [0,1] and 0 elsewhere.  Further Vn* is         //
//     normalized so that Vn*(1) = 1.                                         //
//             <Vn*,Vm*> = 0    if n != m,                                    //
//             <Vn*,Vn*> = pi/2 if n >= 0.                                    //
//     This routine calculates V*[n](x), the shifted Chebyshev polynomials    //
//     of the third kind evaluated at x.                                      //
//     This routine calculates Vn*(x) using the following methods:            //
//     If n < N or if |2x-1| > 1, where N is defined above then use recursion //
//              V*[k+1](x) = 2(2x - 1) V*[k](x) - V*[k-1](x), k = 1,...,n-1   //
//              V*[0](x) = 1, V*[1](x) = 4x - 3.                              //
//     otherwise use V*[n](x) = cos((n+1/2)*theta)/cos(theta/2),              //
//     where theta = acos(2x-1).                                              //
//                                                                            //
//  Arguments:                                                                //
//     long double x                                                          //
//        The argument of the shifted Chebyshev polynomial Vn*.               //
//     int    n                                                               //
//        The degree of the shifted Chebyshev polynomial.                     //
//                                                                            //
//  Return Value:                                                             //
//     Vn*(x) if n is a nonnegative integer.  If n is negative, 0 is returned.//
//                                                                            //
//  Example:                                                                  //
//     long double x, Vn;                                                     //
//     int    n;                                                              //
//                                                                            //
//     (user code to set x and n)                                             //
//                                                                            //
//     Vn = xChebyshev_Shifted_Vn(x, n);                                      //
////////////////////////////////////////////////////////////////////////////////
long double xChebyshev_Shifted_Vn(long double x, int n)
{
   long double two_x_m1 = x + x - 1.0L;
   long double four_x_m2 = two_x_m1 + two_x_m1;
   long double theta;
   long double sin_theta2;
   long double cos_theta2;
   long double V0, V1, Vn;
   int k;

   if (n < 0) return 0.0L;
   if ( x == 1.0L ) return 1.0L;
   if ( x == 0.0L )
      return (n % 2 == 0) ? (long double)(n+n+1) : -(long double)(n+n+1);
    
   if ( (n > N) && (fabsl(two_x_m1) < 1.0L) ) {
      theta = acosl(two_x_m1);
      sin_theta2 = sinl(theta/2.0L);
      cos_theta2 = cosl(theta/2.0L);
      if (sin_theta2 != 1.0L)
         return cosl(((long double)(n) + 0.5L) * theta) / cos_theta2;
      else 
         return (n % 2 == 0) ? (long double) (n+n+1) : -(long double) (n+n+1);
   }
      
                    // Initialize the recursion process. //
 
   if (n == 0) return 1.0L;
   if (n == 1) return four_x_m2 - 1.0L;
   V0 = 1.0L;
   V1 = four_x_m2 - 1.0L;

                             // Calculate V*nx) //

   for (k = 2; k <= n; k++) {
      Vn = four_x_m2 * V1 - V0;
      V0 = V1;
      V1 = Vn;
   }

   return Vn;
}
