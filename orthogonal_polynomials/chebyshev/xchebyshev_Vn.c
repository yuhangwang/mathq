////////////////////////////////////////////////////////////////////////////////
// File: xchebyshev_Vn.c                                                      //
// Routine(s):                                                                //
//    xChebyshev_Vn                                                           //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>         // required for fabsl(), acosl(), sinl(), cosl() 
#define N 8               // cross-over between recursion and using       
                              //   Vn = cos((n+1/2)theta)/cos(theta/2),
                              //   where x = cos(theta)

////////////////////////////////////////////////////////////////////////////////
// long double xChebyshev_Vn(long double x, int n)                            //
//                                                                            //
//  Description:                                                              //
//     The Chebyshev polynomials of the third kind, Vn(x) are orthogonal on   //
//     the interval [-1,1] with weight function w(x) = sqrt((1+x)/(1-x)) on   //
//     interval [-1,1] and 0 elsewhere and are normalized so that Vn(1) = 1.  //
//             <Vn,Vm> = 0    if n != m,                                      //
//             <Vn,Vn> = pi   if n >= 0.                                      //
//     This routine calculates V[n](x) using the following techniques:        //
//     If n < N or if |x| > 1, where N is defined below then use the recursion//
//              V[k+1](x) = 2x V[k](x) - V[k-1](x), k = 1,2,...,n-1           //
//              V[0](x) = 1, V[1](x) = 2x-1.                                  //
//     otherwise use V[n](x) = cos((n+1/2)*theta)/cos(theta/2), where         //
//     theta = acos(x).                                                       //
//                                                                            //
//  Arguments:                                                                //
//     long double x                                                          //
//        The argument of the Chebyshev polynomial Vn.                        //
//     int    n                                                               //
//        The degree of the Chebyshev polynomial.                             //
//                                                                            //
//  Return Value:                                                             //
//     Vn(x) if n is a nonnegative integer.  If n is negative, 0 is returned. //
//                                                                            //
//  Example:                                                                  //
//     long double x, Vn;                                                     //
//     int    n;                                                              //
//                                                                            //
//     (user code to set x and n)                                             //
//                                                                            //
//     Vn = xChebyshev_Vn(x, n);                                              //
////////////////////////////////////////////////////////////////////////////////
long double xChebyshev_Vn(long double x, int n)
{
   long double two_x = x + x;
   long double theta;
   long double cos_theta2;
   long double sin_theta2;
   long double V0, V1, Vn;
   int k;

   if (n < 0) return 0.0L;
   if ( fabsl(x) == 1.0L ) { 
      if (x > 0.0L) return 1.0L;
      else
         return (n % 2 == 0) ? (long double) (n+n+1) : -(long double) (n+n+1);
   }
    
   if ( (n > N) && (fabsl(x) < 1.0L) ) {
      theta = acosl(x);
      cos_theta2 = cosl(theta/2.0L);
      sin_theta2 = sinl(theta/2.0L);
      if (sin_theta2 != 1.0L)
         return cosl(((long double)(n) + 0.5L) * theta) / cos_theta2;
      else 
         return (n % 2 == 0) ? (long double) (n+n+1) : -(long double) (n+n+1);
   }
      
                    // Initialize the recursion process. //
 
   if (n == 0) return 1.0L;
   if (n == 1) return two_x - 1.0L;
   V0 = 1.0L;
   V1 = two_x - 1.0L;

                             // Calculate Vn(x) //

   for (k = 2; k <= n; k++) {
      Vn = two_x * V1 - V0;
      V0 = V1;
      V1 = Vn;
   }

   return Vn;
}
