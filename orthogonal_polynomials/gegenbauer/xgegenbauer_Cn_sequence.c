////////////////////////////////////////////////////////////////////////////////
// File: xgegenbauer_Cn_sequence.c                                            //
// Routine(s):                                                                //
//    xGegenbauer_Cn_Sequence                                                 //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// void xGegembauer_Cn_Sequence(long double C[], long double x,               //
//                                             long double alpha, int max_n)  //
//                                                                            //
//  Description:                                                              //
//     The Gegenbauer polynomials with parameter alpha > -1/2, C^(alpha)n(x)  //
//     are orthogonal on the interval [-1,1] with weight function             //
//     w(x) = (1-x^2)^(alpha-1/2) on [-1,1] and 0 elsewhere.  Further they are//
//     normalized so that C^(alpha)n(1) = gamma(n+2alpha)/(n! gamma(2*alpha)).//
//     Gegenbauer polynomials are also called the ultraspherical polynomials. //
//     For convenience in typing, the superscript (alpha) is dropped so that  //
//     C^(alpha)n(x) will be denoted by Cn(x), i.e. the parameter alpha is    //
//     to be understood.                                                      //
//      <Cn,Cm> = 0                                      if n != m,           //
//      <Cn,Cn> = 2^(1-2alpha) gamma(n+2alpha) pi                             //
//                / (n! (n+alpha)gamma^2(alpha) )        if n >= 0.           //
//     This routine calculates, C^(alpha)[n](x), the Gegenbauer polynomials   //
//     with parameter alpha > -1/2 evaluated at x of degree n for             //
//     n = 0,...,max_n.  The Gegenbauer polynomials are also called the       //
//     ultraspherical polynomials.                                            //
//     This routine calculates, Cn(x), the Gegenbauer polynomial of degree n  //
//     with parameter alpha evaluated at x  using recursion formula:          //
//     (k+1) C[k+1](x) = 2(k + alpha) x C[k](x) - (k + 2 alpha - 1) C[k-1](x),//
//                                                         k = 1,...,max_n-1  //
//     C[0](x) = 1, C[1](x) = 2*alpha*x.                                      //
//                                                                            //
//  Arguments:                                                                //
//     long double C[]                                                        //
//        On output, C[n] is the value of the Gegenbauer polynomial, Cn(),    //
//        evaluated at x where 0 <= n <= max_n.  The calling routine must     //
//        have defined C as long double C[N] where N >= max_n + 1.            //
//     double x                                                               //
//        The argument of the Gegenbauer polynomials with parameter alpha of  //
//        degree n, 0 <= n <= max_n.                                          //
//     double alpha                                                           //
//        The parameter of the Gegenbauer polynomials Cn, 0 <= n <= max_n.    //
//     int    max_n                                                           //
//        The maximum order of the sequence.                                  //
//                                                                            //
//  Return Value:                                                             //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     long double Cn[N];                                                     //
//     long double x;                                                         //
//     long double alpha;                                                     //
//     int    max_deg = N - 1;                                                //
//                                                                            //
//     (user code to set x and alpha)                                         //
//                                                                            //
//     xGegenbauer_Cn_Sequence(Cn, x, alpha, max_deg);                        //
////////////////////////////////////////////////////////////////////////////////
void xGegenbauer_Cn_Sequence(long double C[], long double x, long double alpha,
                                                                     int max_n)
{
   long double two_alpha = 2.0L * alpha;
   long double kp1;
   long double beta;
   long double gamma;
   int k;
                    // Initialize the recursion process. //

   if (max_n < 0) return; 
   C[0] = 1.0L;
   if (max_n == 0) return;
   C[1] = two_alpha * x;
   if (max_n == 1) return;

                   // Calculate Cn(x) for n = 2,...,max_n //

   for (k = 1; k < max_n; k++) {
      kp1 = (long double)(k + 1);
      beta = ( (long double)(k + k) + two_alpha ) * x;
      gamma = (long double)(k - 1) + two_alpha;
      C[k+1] = ( beta * C[k] - gamma * C[k-1] ) / kp1;
   }
}
