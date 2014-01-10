////////////////////////////////////////////////////////////////////////////////
// File: xkrawtchouk_Kn.c                                                     //
// Routine(s):                                                                //
//    xKrawtchouk_Kn                                                          //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// long double xKrawtchouk_Kn(long double x, long double p, int N, int n)     //
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
//     long double x                                                          //
//        The argument of the Krawtchouk polynomials with weight function     //
//        function support points {0,1,2,...,N} and values Bin(N,p,x).        //
//     long double p                                                          //
//        The parameter p of the binomial weight function, 0 < p < 1.         //
//     int    N                                                               //
//        The parameter N of the binomial weight function.  The number of     //
//        equally spaced support points of the weight function is N + 1.      //
//     int    n                                                               //
//        The degree of the discrete Krawtchouk polynomial.                   //
//                                                                            //
//  Return Value:                                                             //
//     Kn(x) if n is a nonnegative integer otherwise 0 is returned.           //
//                                                                            //
//  Example:                                                                  //
//     long double x,p,Kn                                                     //
//     int    n,N;                                                            //
//                                                                            //
//     (user code to set x, p, N, and n )                                     //
//                                                                            //
//     Kn = xKrawtchouk_Kn(x, p, N, n);                                       //
////////////////////////////////////////////////////////////////////////////////
long double xKrawtchouk_Kn(long double x, long double p, int N, int n)
{
   long double k0, k1, kn;
   long double a, b;
   long double pq = p * (1.0L - p);
   int k;
                    // Initialize the recursion process. //

   if (n < 0) return 0.0L;
   if (n > N) return 0.0L;
   k0 = 1.0L;
   if (n == 0) return k0;
   k1 = x - (long double)N * p;
   if (n == 1) return k1;

                             // Calculate Kn(x) //
   for (k = 1; k < n; k++) {
      a = x - ( (long double) k + p * (long double)( N - k - k));
      b = (long double)(N - k + 1) * pq;
      kn = ( a * k1 - b * k0 ) / (long double)(k + 1);
      k0 = k1;
      k1 = kn;
   }

   return kn;
}
