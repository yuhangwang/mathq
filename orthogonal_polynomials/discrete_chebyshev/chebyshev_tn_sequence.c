////////////////////////////////////////////////////////////////////////////////
// File: chebyshev_tn_sequence.c                                              //
// Routine(s):                                                                //
//    Chebyshev_tn_Sequence                                                   //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// void Chebyshev_tn_Sequence(double t[], double x, int N, int max_n)         //
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
//     double t[]                                                             //
//        On output, t[n] is the value of the discrete Chebyshev polynomial,  //
//        tn(), evaluated at x where 0 <= n <= max_n and the support points   //
//        of the weight function are {0,1,2,...,N-1}.  The calling routine    //
//        must have defined t as double t[K] where K >= max_n + 1.            //
//     double x                                                               //
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
//     double tn[K];                                                          //
//     double x;                                                              //
//     int    N;                                                              //
//     int    max_deg = K - 1;                                                //
//                                                                            //
//     (user code to set x and N)                                             //
//                                                                            //
//     Chebyshev_tn_Sequence(tn, x, N, max_deg);                              //
////////////////////////////////////////////////////////////////////////////////
void Chebyshev_tn_Sequence(double t[], double x, int N, int max_n)
{
   long double a = 2.0L * (long double)x - (long double)(N-1);
   long double tn, tn1, tn2;
   int N2 = N * N;
   int k;

   if (max_n < 0) return;
   if (max_n >= N-1) max_n = N - 1;

                    // Initialize the recursion process. //
 
   tn2 = 1.0L;
   t[0] = (double) tn2;
   if (max_n == 0) return;

   tn1 = a;
   t[1] = (double) tn1;
   if (max_n == 1) return;

                   // Calculate tn(x) for n = 2,...,max_n //

   for (k = 1; k < max_n; k++) {
      tn = (((long double)(k + k + 1) * a) * tn1 
                      - ((long double) k * (long double) (N2 - k * k)) * tn2);
      tn /= (long double)(k + 1);
      t[k+1] = (double)tn;
      tn2 = tn1;
      tn1 = tn;
   }
}
