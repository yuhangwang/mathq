////////////////////////////////////////////////////////////////////////////////
// File: xchebyshev_tn.c                                                      //
// Routine(s):                                                                //
//    xChebyshev_tn                                                           //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// long double xChebyshev_tn(long double x, int N, int n)                     //
//                                                                            //
//  Description:                                                              //
//     The discrete Chebyshev polynomials are orthogonal on the interval      //
//     [0,N-1] with discrete uniform weight function w(x) = Sum [D(x-j)],     //
//     where D() is the Dirac delta function and the sum is over              //
//     j = 0, 1, ... N - 1.                                                   //
//     I.e. the support points of the weight function are {0,1,...,N-1}.      //
//     This routine calculates, t[n](x), the n-th discrete Chebyshev          //
//     polynomial evaluated at x using the following recursion formula:       //
//      (k+1) t[k+1](x) = (2k+1)(2x-(N-1)) t[k](x) - k(N^2 - k^2)t[k-1](x)    //
//      t[0](x) = 1, t[1](x) = 2x-(N-1).                                      //
//     Note that only t[n](x), n = 0,...,N-1 are orthogonal wrt the weight    //
//     function w(x) given above.                                             //
//                                                                            //
//  Arguments:                                                                //
//     long double x                                                          //
//        The argument of the discrete Chebyshev polynomial.                  //
//     int    N                                                               //
//        The number of equally spaced support points of the weight function. //
//     int    n                                                               //
//        The degree of the discrete Chebyshev polynomial.                    //
//                                                                            //
//  Return Value:                                                             //
//     tn(x) if n is a nonnegative integer less than N.  If n is negative or  //
//     if n >= N then 0 is returned.                                          //
//                                                                            //
//  Example:                                                                  //
//     long double x, tn;                                                     //
//     int         n,N;                                                       //
//                                                                            //
//     (user code to set x, N, and n)                                         //
//                                                                            //
//     tn = xChebyshev_tn( x, N, n);                                          //
////////////////////////////////////////////////////////////////////////////////
long double xChebyshev_tn(long double x, int N, int n)
{
   long double a = 2.0L * x - (long double)(N - 1);
   long double t0, t1, tn;
   int N2 = N * N;
   int k;

   if (n < 0) return 0.0L;
   if (n >= N) return 0.0L;
      
                    // Initialize the recursion process. //
 
   t0 = 1.0L;
   if (n == 0) return t0;
   t1 = a;
   if (n == 1) return t1;

                             // Calculate tn(x) //

   for (k = 1; k < n; k++) {
      tn = ( ( (long double)(k + k + 1) * a ) * t1
             - (long double)k * (long double) (N2 - k * k)  * t0 )
                                                      / (long double)(k + 1);
      t0 = t1;
      t1 = tn;
   }

   return tn;
}
