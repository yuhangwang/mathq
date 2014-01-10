////////////////////////////////////////////////////////////////////////////////
// File: xchebyshev_tn_sequence.c                                             //
// Routine(s):                                                                //
//    xChebyshev_tn_Sequence                                                  //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// void xChebyshev_tn_Sequence(long double t[], long double x, int N,         //
//                                                                int max_n)  //
//                                                                            //
//  Description:                                                              //
//     The discrete Chebyshev  polynomials are orthogonal on the interval     //
//     [0,N-1] with weight function w(x) = Sum [D(x-j)], where D() is the     //
//     Dirac delta function, and the sum is over j = 0, 1, ... .              //
//     I.e. the support points of the weight function are {0,1,...,N-1}.      //
//     This routine calculates, t[n](x), the n-th discrete Chebyshev          //
//     polynomial evaluated at x for n = 0, ,,,, max_n.                       //
//     These procedures use the recursion formula:                            //
//      (k+1) t[k+1](x) = (2k+1)(2x-(N-1)) t[k](x) - k(N^2 - k^2)t[k-1](x)    //
//      t[0](x) = 1, t[1](x) = 2x-(N-1).                                      //
//                                                                            //
//  Arguments:                                                                //
//     long double t[]                                                        //
//        On output, t[n] is the value of the discrete Chebyshev polynomial,  //
//        tn(), evaluated at x where 0 <= n <= max_n and the support points   //
//        of the weight function are {0,1,2,...,N-1}.  The calling routine    //
//        must have defined t as long double t[K] where K >= max_n + 1.       //
//     long double x                                                          //
//        The argument of the discrete Chebyshev polynomials with weight      //
//        function support points {0,1,2,...,N-1}.                            //
//     int    N                                                               //
//        The number of equally spaced support points of the weight function. //
//     int    max_n                                                           //
//        The maximum order of the sequence, 0 <= max_n <= N-1.               //
//                                                                            //
//  Return Value:                                                             //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     long double tn[K];                                                     //
//     long double x;                                                         //
//     int    N;                                                              //
//     int    max_deg = K - 1;                                                //
//                                                                            //
//     (user code to set x and N)                                             //
//                                                                            //
//     xChebyshev_tn_Sequence(tn, x, N, max_deg);                             //
////////////////////////////////////////////////////////////////////////////////
void xChebyshev_tn_Sequence(long double t[], long double x, int N, int max_n)
{
   long double a = 2.0L * (long double)x - (long double)(N-1);
   int N2 = N * N;
   int k;

   if (max_n < 0) return;
   if (max_n >= N-1) max_n = N - 1;

                    // Initialize the recursion process. //
 
   t[0] = 1.0L;
   if (max_n == 0) return;

   t[1] = a;
   if (max_n == 1) return;

                   // Calculate tn(x) for n = 2,...,max_n //

   for (k = 1; k < max_n; k++) {
      t[k+1] = (((long double)(k + k + 1) * a) * t[k] 
                    - ((long double) k * (long double) (N2 - k * k)) * t[k-1]);
      t[k+1] /= (long double)(k + 1);
   }
}
