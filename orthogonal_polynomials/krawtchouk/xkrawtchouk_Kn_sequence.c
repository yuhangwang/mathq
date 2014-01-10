////////////////////////////////////////////////////////////////////////////////
// File: xkrawtchouk_Kn_sequence.c                                            //
// Routine(s):                                                                //
//    xKrawtchouk_Kn_Sequence                                                 //
////////////////////////////////////////////////////////////////////////////////
#define Krawtchouk_Step(a, kn1, b, kn2) ( a * kn1 - b * kn2 )
////////////////////////////////////////////////////////////////////////////////
// void xKrawtchouk_Kn(long double K[], long double x, long double p, int N,  //
//                                                                 int max_n) //
//                                                                            //
//  Description:                                                              //
//     The Krawtchouk polynomials are orthogonal on the interval [0,N] with   //
//     binomial distributed weight function w(x) = Sum Bin(N,p,j) [D(x-j)],   //
//     where D() is the Dirac delta function, Bin(N,p,j) is the binomial term //
//     Bin(N,p,j) = C(N,j) p^j (1-p)^(N-j) where C(N,j) is the combination of //
//     N objects taken j at a time, and the sum is over j = 0, 1, ... N.      //
//     I.e. the support points of the weight function are {0,1,...,N} with    //
//     values Bin(N,p,j) at x = j.                                            //
//     This routine calculatex, K[n](x), the n-th Krawtchouk polynomial       //
//     evaluated at x for n = 0, ,,,, max_n.                                  //
//     This procedure uses the recursion formula:                             //
//      (k+1) K[k+1](x) = [x - (k+p(N-2k))] K[k](x) - (N-n+1)p(1-p) K[k-1](x) //
//      K[0](x) = 1, K[1](x) = x-pN,    k = 1,...,max_n-1.                    //
//                                                                            //
//  Arguments:                                                                //
//     long double K[]                                                        //
//        On output, K[n] is the value of the Krawtchouk polynomial Kn(),     //
//        evaluated at x where 0 <= n <= max_n and the support points         //
//        of the weight function are {0,1,2,...,N}.  The calling routine      //
//        must have defined K as double K[L] where L >= max_n + 1.            //
//     long double x                                                          //
//        The argument of the Krawtchouk polynomials Kn with weight function  //
//        function support points {0,1,2,...,N} and values Bin(N,p,x)         //
//        for n = 0, ..., max_n.                                              //
//     long double p                                                          //
//        A parameter of the Krawtchouk polynomials Kn, n = 0, ..., max_n     //
//        where 0 < p < 1.                                                    //
//     int    N                                                               //
//        N + 1 is the number of equally spaced support points of the weight  //
//        function.                                                           //
//     int    max_n                                                           //
//        The maximum order of the sequence, 0 <= max_n <= N.                 //
//                                                                            //
//  Return Value:                                                             //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     double Kn[L];                                                          //
//     double x;                                                              //
//     double p;                                                              //
//     int    N;                                                              //
//     int    max_deg = L - 1;                                                //
//                                                                            //
//     (user code to set x, p, and N )                                        //
//                                                                            //
//     xKrawtchouk_Kn_Sequence(Kn, x, p, N, max_deg);                         //
////////////////////////////////////////////////////////////////////////////////
void xKrawtchouk_Kn_Sequence(long double K[], long double x, long double p, 
                                                              int N, int max_n)
{
   long double a;
   long double b;
   long double pq = (1.0L - (long double) p) * (long double) p;
   int k;

   if (max_n < 0) return;
   if (max_n > N) max_n = N;

                    // Initialize the recursion process. //
 
   K[0] = 1.0L;
   if (max_n == 0) return;

   K[1] = x - (long double)N * p;
   if (max_n == 1) return;

                   // Calculate Kn(x) for n = 2,...,max_n //

   for (k = 1; k < max_n; k++) {
      a = x - ( (long double) k + p * (long double)( N - k - k));
      b = (long double)(N - k + 1) * pq;
      K[k+1] = ( Krawtchouk_Step(a, K[k], b, K[k-1]) ) / (long double)(k + 1);
   }
}
