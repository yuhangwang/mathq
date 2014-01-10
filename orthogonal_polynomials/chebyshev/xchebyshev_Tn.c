////////////////////////////////////////////////////////////////////////////////
// File: xchebyshev_Tn.c                                                      //
// Routine(s):                                                                //
//    xChebyshev_Tn                                                           //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>             // required for fabsl(), acosl(), and cosl()
#define N 6                   // cross-over between recursion and cos(acos(x))

//                        Internally Defined Routines                         //

long double xChebyshev_Tn(long double x, int n);

////////////////////////////////////////////////////////////////////////////////
// long double xChebyshev_Tn(long double x, int n)                            //
//                                                                            //
//  Description:                                                              //
//     The Chebyshev polynomials of the first kind, Tn(x), are orthogonal on  //
//     the interval [-1,1] with weight function w(x) = 1 / sqrt(1-x^2) on the //
//     interval (-1,1) and 0 elsewhere and are normalized so that Tn(1) = 1.  //
//             <Tn,Tm> = 0    if n != m,                                      //
//             <Tn,Tn> = pi/2 if n > 0,                                       //
//             <T0,T0> = pi.                                                  //
//     This routine calculate Tn(x) using the following techniques:           //
//     If n < N or if |x| > 1, where N is defined below then use the recursion//
//              T[k+1](x) = 2x T[k](x) - T[k-1](x), k = 1,2,...,n-1           //
//              T[0](x) = 1, T[1](x) = x                                      //
//     otherwise use the explicit formula T[n](x) = cos(n*acos(x)).           //
//     This routine calculates T[n](x), the Chebyshev polynomial of the       //
//     first kind of degree n evaluated at x where x and the return value are //
//     of type long double.                                                   //
//                                                                            //
//  Arguments:                                                                //
//     long double x                                                          //
//        The argument of the Chebyshev polynomial Tn.                        //
//     int    n                                                               //
//        The degree of the Chebyshev polynomial.                             //
//                                                                            //
//  Return Value:                                                             //
//     Tn(x) if n is a nonnegative integer.  If n is negative, 0 is returned. //
//                                                                            //
//  Example:                                                                  //
//     long double x, Tn;                                                     //
//     int    n;                                                              //
//                                                                            //
//     (user code to set x and n)                                             //
//                                                                            //
//     Tn = xChebyshev_Tn(x, n);                                              //
////////////////////////////////////////////////////////////////////////////////
long double xChebyshev_Tn(long double x, int n)
{
   long double two_x = x + x;
   long double theta;
   long double T0, T1, Tn;
   int k;

   if (n < 0) return 0.0L;
   if ( fabsl(x) == 1.0L ) { 
      if (x > 0.0L) return 1.0L;
      else if (n % 2 == 0) return 1.0L;
      else return -1.0L;
   }
    
   if ( (n > N) && (fabsl(x) < 1.0L) ) {
      theta = acosl(x);
      return cosl(n * theta);
   }
      
                    // Initialize the recursion process. //
 
   if (n == 0) return 1.0L;
   if (n == 1) return x;
   T0 = 1.0L;
   T1 = x;
                             // Calculate Tn(x) //

   for (k = 2; k <= n; k++) {
      Tn = two_x * T1 - T0;
      T0 = T1;
      T1 = Tn;
   }

   return Tn;
}
