////////////////////////////////////////////////////////////////////////////////
// File: chebyshev_tn.c                                                       //
// Routine(s):                                                                //
//    Chebyshev_tn                                                            //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>             // required for fabsl()
#include <float.h>            // required for DBL_MAX

//                        Externally Defined Routines                         //

long double xChebyshev_tn(long double x, int N, int n);

////////////////////////////////////////////////////////////////////////////////
// double Chebyshev_tn_(double x, int N, int n)                               //
//                                                                            //
//  Description:                                                              //
//     The discrete Chebyshev polynomials are orthogonal on the interval      //
//     [0,N-1] with discrete uniform weight function w(x) = Sum [D(x-j)],     //
//     where D() is the Dirac delta function and the sum is over              //
//     j = 0, 1, ... N - 1.                                                   //
//     I.e. the support points of the weight function are {0,1,...,N-1}.      //
//     This routine calculates, t[n](x), the n-th discrete Chebyshev          //
//     polynomial evaluated at x using the function xChebyshev_tn in the file //
//     xchebyshev_tn.c which in turn uses the following recursion formula:    //
//      (k+1) t[k+1](x) = (2k+1)(2x-(N-1)) t[k](x) - k(N^2 - k^2)t[k-1](x)    //
//      t[0](x) = 1, t[1](x) = 2x-(N-1).                                      //
//     Note that only t[n](x), n = 0,...,N-1 are orthogonal wrt the weight    //
//                                                                            //
//  Arguments:                                                                //
//     double x                                                               //
//        The argument of the discrete Chebyshev polynomials with weight      //
//        function support points {0,1,2,...,N-1}.                            //
//     int    N                                                               //
//        The number of equally spaced support points of the weight function. //
//     int    n                                                               //
//        The degree of the discrete Chebyshev polynomial, 0 <= n <= N-1.     //
//                                                                            //
//  Return Value:                                                             //
//     tn(x) if n is a nonnegative integer less than N.  If n is negative     //
//     then 0 is returned.                                                    //
//     If tn(x) > DBL_MAX, then DBL_MAX is returned and if tn(x) < -DBL_MAX   //
//     then -DBL_MAX is returned.                                             //
//                                                                            //
//  Example:                                                                  //
//     double x,tn;                                                           //
//     int    n,N;                                                            //
//                                                                            //
//     (user code to set x, N, and n )                                        //
//                                                                            //
//     tn = Chebyshev_tn(x, N, n);                                            //
////////////////////////////////////////////////////////////////////////////////
double Chebyshev_tn(double x, int N, int n)
{
   long double tn;

   if (n < 0) return 0.0;
   tn = xChebyshev_tn((long double) x, N, n);
   if (fabsl(tn) < DBL_MAX) return (double) tn;
   return (tn > 0.0L) ? DBL_MAX : -DBL_MAX;

}
