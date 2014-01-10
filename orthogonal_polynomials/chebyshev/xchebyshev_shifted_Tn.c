////////////////////////////////////////////////////////////////////////////////
// File: xchebyshev_shifted_Tn.c                                              //
// Routine(s):                                                                //
//    xChebyshev_Shifted_Tn                                                   //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>             // required for fabsl(), acosl(), and cosl()
#define N 6                   // cross-over between recursion and cos(acos(x))

////////////////////////////////////////////////////////////////////////////////
// long double xChebyshev_Shifted_Tn(long double x, int n)                    //
//                                                                            //
//  Description:                                                              //
//     The shifted Chebyshev polynomials of the first kind, Tn*(x) = Tn(2x-1),//
//     are orthogonal on the interval [0,1] with weight function              //
//     w(x) = 1 / sqrt(x-x^2) on [0,1] and 0 elsewhere.  Tn* is normalized    //
//     so that T*n(1) = 1.                                                    //
//             <Tn*,Tm*> = 0    if n != m,                                    //
//             <Tn*,Tn*> = pi/2 if n > 0,                                     //
//             <T0*,T0*> = pi.                                                //
//     This routine calculates Tn*(x) using the following techniques:         //
//     If n < N or if |2x-1| > 1, where N is defined above then use recursion //
//              T*[k+1](x) = 2(2x - 1) T*[k](x) - T*[k-1](x), k = 1,...,n-1   //
//              T*[0](x) = 1, T*[1](x) = 2x - 1.                              //
//     otherwise use the explicit formula T*[n](x) = cos(n*acos(2x-1))        //
//                                                                            //
//  Arguments:                                                                //
//     long double x                                                          //
//        The argument of the shifted Chebyshev polynomial Tn*.               //
//     int    n                                                               //
//        The degree of the shifted Chebyshev polynomial.                     //
//                                                                            //
//  Return Value:                                                             //
//     Tn*(x) if n is a nonnegative integer.  If n is negative, 0 is returned.//
//                                                                            //
//  Example:                                                                  //
//     long double x, Tn;                                                     //
//     int    n;                                                              //
//                                                                            //
//     (user code to set x and n)                                             //
//                                                                            //
//     Tn = xChebyshev_Shifted_Tn(x, n);                                      //
////////////////////////////////////////////////////////////////////////////////
long double xChebyshev_Shifted_Tn(long double x, int n)
{
   long double two_x_m1 = x + x - 1.0L;
   long double four_x_m2 = two_x_m1 + two_x_m1;
   long double theta;
   long double T0, T1, Tn;
   int k;

   if (n < 0) return 0.0L;
   if ( x == 1.0L ) return 1.0L;
   if ( x == 0.0L ) return (n % 2 == 0) ? 1.0L : -1.0L;
    
   if ( (n > N) && (fabsl(two_x_m1) < 1.0L) ) {
      theta = acosl(two_x_m1);
      return cosl(n * theta);
   }
      
                    // Initialize the recursion process. //
 
   if (n == 0) return 1.0L;
   if (n == 1) return two_x_m1;
   T0 = 1.0L;
   T1 = two_x_m1;

                             // Calculate T*nx) //

   for (k = 2; k <= n; k++) {
      Tn = four_x_m2 * T1 - T0;
      T0 = T1;
      T1 = Tn;
   }

   return Tn;
}
