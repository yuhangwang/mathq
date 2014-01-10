////////////////////////////////////////////////////////////////////////////////
// File: krawtchouk_Kn_sequence.c                                             //
// Routine(s):                                                                //
//    Krawtchouk_Kn_Sequence                                                  //
////////////////////////////////////////////////////////////////////////////////
#define Krawtchouk_Step(a, kn1, b, kn2) ( a * kn1 - b * kn2 )
////////////////////////////////////////////////////////////////////////////////
// void Krawtchouk_Kn(double K[], double x, double p, int N, int max_n)       //
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
//     double K[]                                                             //
//        On output, K[n] is the value of the Krawtchouk polynomial Kn(),     //
//        evaluated at x where 0 <= n <= max_n and the support points         //
//        of the weight function are {0,1,2,...,N}.  The calling routine      //
//        must have defined K as double K[L] where L >= max_n + 1.            //
//     double x                                                               //
//        The argument of the Krawtchouk polynomials Kn with weight function  //
//        function support points {0,1,2,...,N} and values Bin(N,p,x)         //
//        for n = 0, ..., max_n.                                              //
//     double p                                                               //
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
//     Krawtchouk_Kn_sequence(Kn, x, p, N, max_deg);                          //
////////////////////////////////////////////////////////////////////////////////
void Krawtchouk_Kn_Sequence(double K[], double x, double p, int N, int max_n)
{
   long double a;
   long double b;
   long double pq = (1.0L - (long double) p) * (long double) p;
   long double kn, kn1, kn2;
   int k;

   if (max_n < 0) return;
   if (max_n > N) max_n = N;

                    // Initialize the recursion process. //
 
   kn2 = 1.0L;
   K[0] = (double) kn2;
   if (max_n == 0) return;

   kn1 = (long double) x - (long double)N * (long double)p;
   K[1] = (double) kn1;
   if (max_n == 1) return;

                   // Calculate Kn(x) for n = 2,...,max_n //

   for (k = 1; k < max_n; k++) {
      a = (long double) x - ( (long double) k
                              + (long double) p * (long double)( N - k - k));
      b = (long double)(N - k + 1) * pq;
      kn = Krawtchouk_Step(a, kn1, b, kn2);
      kn /= (long double)(k + 1);
      K[k+1] = (double)kn;
      kn2 = kn1;
      kn1 = kn;
   }
}
