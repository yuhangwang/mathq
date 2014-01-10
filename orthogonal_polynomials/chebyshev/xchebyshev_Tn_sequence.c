////////////////////////////////////////////////////////////////////////////////
// File: xchebyshev_Tn_sequence.c                                             //
// Routine(s):                                                                //
//    xChebyshev_Tn_Sequence                                                  //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// void xChebyshev_Tn_Sequence(long double T[], long double x, int max_n)     //
//                                                                            //
//  Description:                                                              //
//     The Chebyshev polynomials of the first kind, Tn(x), are orthogonal on  //
//     the interval [-1,1] with weight function w(x) = 1 / sqrt(1-x^2) on the //
//     interval [-1,1] and 0 elsewhere and are normalized so that Tn(1) = 1.  //
//             <Tn,Tm> = 0    if n != m,                                      //
//             <Tn,Tn> = pi/2 if n > 0,                                       //
//             <T0,T0> = pi.                                                  //
//     This routine calculates, T[n](x), the Chebyshev polynomials of the     //
//     first kind evaluated at x for n = 0, ..., max_n using the recursion    //
//     This procedure uses the recursion:                                     //
//              T[n+1](x) = 2x T[n](x) - T[n-1](x), n = 1,...,max_n - 1       //
//              T[0](x) = 1, T[1](x) = x.                                     //
//                                                                            //
//  Arguments:                                                                //
//     long double T[]                                                        //
//        On output, T[n] is the value of the Chebyshev polynomial, Tn(),     //
//        evaluated at x where 0 <= n <= max_n.  The calling routine must have//
//        defined T as long double T[N] where N >= max_n + 1.                 //
//     double x                                                               //
//        The argument of the Chebyshev polynomials Tn, n = 0,...,max_n.      //
//     int    max_n                                                           //
//        The maximum order of the sequence.                                  //
//                                                                            //
//  Return Value:                                                             //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     long double x, Tn[N];                                                  //
//     int    max_deg = N - 1;                                                //
//                                                                            //
//     (user code to set x)                                                   //
//                                                                            //
//     xChebyshev_Tn_Sequence(Tn, x, max_deg);                                //
////////////////////////////////////////////////////////////////////////////////
void xChebyshev_Tn_Sequence(long double T[], long double x, int max_n)
{
   long double two_x = x + x;
   int k;
                    // Initialize the recursion process. //
 
   if (max_n < 0) return;
   T[0] = 1.0L;
   if (max_n == 0) return;
   T[1] = x;
   if (max_n == 1) return;

                   // Calculate Tn(x) for n = 2,...,max_n //

   for (k = 2; k <= max_n; k++) 
      T[k] = two_x * T[k-1] - T[k-2];
}
