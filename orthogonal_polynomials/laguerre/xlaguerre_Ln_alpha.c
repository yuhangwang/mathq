////////////////////////////////////////////////////////////////////////////////
// File: xlaguerre_Ln_alpha.c                                                 //
// Routine(s):                                                                //
//    xLaguerre_Ln_alpha                                                      //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// long double xLaguerre_Ln_alpha(long double x, long double alpha, int n)    //
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
//     This routine calculates, Ln(x), the generalized Laguerre polynomial    //
//     with parameter alpha of degree n evaluated at x using the recursion    //
//     formula:                                                               //
//       (k+1) L[k+1](x) = (2k+alpha+1-x) L[k](x)                             //
//                                  - (k+alpha) L[k-1](x), k = 1,...,n-1      //
//              L[0](x) = 1, L[1](x) = alpha + 1 - x.                         //
//                                                                            //
//  Arguments:                                                                //
//     long double x                                                          //
//        The argument of the generalized Laguerre polynomial with parameter  //
//        alpha.                                                              //
//     long double alpha                                                      //
//        The parameter of the generalized Laguerre polynomial.               //
//     int    n                                                               //
//        The degree of the generalized Laguerre polynomial with parameter    //
//        alpha.                                                              //
//                                                                            //
//  Return Value:                                                             //
//     Ln(x) if n is a nonnegative integer.  If n is negative, 0 is returned. //
//                                                                            //
//  Example:                                                                  //
//     long double x, alpha, Ln;                                              //
//     int    n;                                                              //
//                                                                            //
//     (user code to set x, alpha, and n)                                     //
//                                                                            //
//     Ln = xLaguerre_Ln_alpha( x, alpha, n);                                 //
////////////////////////////////////////////////////////////////////////////////
long double xLaguerre_Ln_alpha(long double x, long double alpha, int n)
{
   long double alpha_one_mx = alpha + 1.0L - x;
   long double L0, L1, Ln;
   int k;

   if (n < 0) return 0.0L;
      
                    // Initialize the recursion process. //
 
   if (n == 0) return 1.0L;
   if (n == 1) return alpha_one_mx;
   L0 = 1.0L;
   L1 = alpha_one_mx;

                             // Calculate Ln(x) //

   for (k = 1; k < n; k++) {
      Ln = ( ( (long double)(k+k) + alpha_one_mx ) * L1
                     - (alpha +(long double)k) * L0 ) / (long double)(k + 1);
      L0 = L1;
      L1 = Ln;
   }

   return Ln;
}
