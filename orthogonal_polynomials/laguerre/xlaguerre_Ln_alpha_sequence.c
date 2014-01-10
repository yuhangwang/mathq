////////////////////////////////////////////////////////////////////////////////
// File: xlaguerre_Ln_alpha_sequence.c                                        //
// Routine(s):                                                                //
//    xLaguerre_Ln_alpha_Sequence                                             //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// void xLaguerre_Ln_alpha_Sequence(long double L[], long double x,           //
//                                              long double alpha, int max_n) //
//                                                                            //
//  Description:                                                              //
//     The generalized Laguerre polynomials with parameter alpha > -1,        //
//     L^(alpha)n(x), are orthogonal on the interval [0,inf) with weight      //
//     function w(x) = x^alpha exp(-x) on [0,inf) and 0 elsewhere.            //
//     For convenience in typing, the superscript (alpha) is dropped so that  //
//     L^(alpha)[n](x) will be denoted by L[n](x), i.e. the alpha is to be    //
//     understood.  The generalized Laguerre polynomials are normalized so    //
//     that Ln(0) = gamma(alpha + n + 1) / (gamma(alpha+1) n!).               //
//             <Ln,Lm> = 0                         if n != m,                 //
//             <Ln,Ln> = gamma(alpha + n + 1) / n! if n >= 0.                 //
//     This routine calculates, L[n](x), Laguerre polynomials evaluated at    //
//     x for n = 0, ,,,, max_n.                                               //
//     These routines calculate, Ln(x), the generalized Laguerre polynomial   //
//     with parameter alpha of degree n evaluated at x using the recursion    //
//     formula:                                                               //
//       (k+1) L[k+1](x) = (2k+alpha+1-x) L[k](x)                             //
//                                  - (k+alpha) L[k-1](x), k = 1,...,max_n-1  //
//              L[0](x) = 1, L[1](x) = alpha + 1 - x.                         //
//                                                                            //
//  Arguments:                                                                //
//     long double L[]                                                        //
//        On output, L[n] is the value of the generalized Laguerre polynomial //
//        with parameter alpha, Ln(), evaluated at x where 0 <= n <= max_n.   //
//        The calling routine must have defined L as long double L[N] where   //
//        N >= max_n + 1.                                                     //
//     long double x                                                          //
//        The argument of the generalized Laguerre polynomials with parameter //
//        alpha, Ln, for n = 0,...,max_n.                                     //
//     long double alpha                                                      //
//        The parameter of the generalized Laguerre polynomials, alpha > -1.  //
//     int    max_n                                                           //
//        The maximum order of the sequence.                                  //
//                                                                            //
//  Return Value:                                                             //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     long double Ln[N];                                                     //
//     long double x;                                                         //
//     long double alpha;                                                     //
//     int    max_deg = N - 1;                                                //
//                                                                            //
//     (user code to set x and alpha)                                         //
//                                                                            //
//     xLaguerre_Ln_alpha_Sequence(Ln, x, alpha, max_deg);                    //
////////////////////////////////////////////////////////////////////////////////
void xLaguerre_Ln_alpha_Sequence(long double L[], long double x,
                                                  long double alpha, int max_n)
{
   long double alpha_one_mx = alpha + 1.0L - x;
   long double beta, gamma;
   int k;

   if (max_n < 0) return;

                    // Initialize the recursion process. //
 
   L[0] = 1.0L;
   if (max_n == 0) return;
   L[1] = alpha_one_mx;
   if (max_n == 1) return;

                   // Calculate Ln(x) for n = 2,...,max_n //

   for (k = 1; k < max_n; k++) {
      beta = (long double)(k + k) + alpha_one_mx;
      gamma = (long double) k + alpha;
      L[k+1] = (beta * L[k] - gamma * L[k-1]) / (long double)(k+1);
   }
}
