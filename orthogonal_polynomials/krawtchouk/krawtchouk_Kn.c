////////////////////////////////////////////////////////////////////////////////
// File: krawtchouk_Kn.c                                                      //
// Routine(s):                                                                //
//    Krawtchouk_Kn                                                           //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>             // required for fabsl()
#include <float.h>            // required for DBL_MAX

//                        Externally Defined Routines                         //

long double xKrawtchouk_Kn(long double x, long double p,int N, int n);

////////////////////////////////////////////////////////////////////////////////
// double Krawtchouk_Kn(double x, double p, int N, int n)                     //
//                                                                            //
//  Description:                                                              //
//     The Krawtchouk polynomials are orthogonal on the interval [0,N] with   //
//     binomial distributed weight function w(x) = Sum Bin(N,p,j) [D(x-j)],   //
//     where D() is the Dirac delta function, Bin(N,p,j) is the binomial term //
//     Bin(N,p,j) = C(N,j) p^j (1-p)^(N-j) where C(N,j) is the combination of //
//     N objects taken j at a time, and the sum is over j = 0, 1, ... N.      //
//     I.e. the support points of the weight function are {0,1,...,N} with    //
//     values Bin(N,p,j) at x = j.                                            //
//     This routine calculates, K[n](x), the n-th Krawtchouk polynomial       //
//     evaluated at x using the following recursion formula:                  //
//      (k+1) K[k+1](x) = [x - (k+p(N-2k))] K[k](x) - (N-n+1)p(1-p) K[k-1](x) //
//      K[0](x) = 1, K[1](x) = x-pN.                                          //
//     Note that only K[n](x), n = 0,...,N are orthogonal wrt the weight      //
//     function w(x) given above.                                             //
//                                                                            //
//  Arguments:                                                                //
//     double x                                                               //
//        The argument of the Krawtchouk polynomials with weight function     //
//        function support points {0,1,2,...,N} and values Bin(N,p,x).        //
//     double p                                                               //
//        The parameter p of the binomial weight function, 0 < p < 1.         //
//     int    N                                                               //
//        The parameter N of the binomial weight function.  The number of     //
//        equally spaced support points of the weight function is N + 1.      //
//     int    n                                                               //
//        The degree of the Krawtchouk polynomial, 0 <= n <= N.               //
//                                                                            //
//  Return Value:                                                             //
//     Kn(x) if n is a nonnegative integer otherwise 0 is returned.           //
//     If Kn(x) > DBL_MAX, then DBL_MAX is returned and if Kn(x) < -DBL_MAX   //
//     then -DBL_MAX is returned.                                             //
//                                                                            //
//  Example:                                                                  //
//     double x,p,Kn;                                                         //
//     int    n,N;                                                            //
//                                                                            //
//     (user code to set x, p, N, and n )                                     //
//                                                                            //
//     Kn = Krawtchouk_Kn(x, p, N, n);                                        //
////////////////////////////////////////////////////////////////////////////////
double Krawtchouk_Kn(double x, double p, int N, int n)
{
   long double Kn;

   Kn = xKrawtchouk_Kn((long double) x, (long double) p, N, n);

   if (fabsl(Kn) < DBL_MAX) return (double) Kn;
   return (Kn > 0.0L) ? DBL_MAX : -DBL_MAX;
}
