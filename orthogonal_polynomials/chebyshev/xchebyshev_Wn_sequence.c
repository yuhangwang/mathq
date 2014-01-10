////////////////////////////////////////////////////////////////////////////////
// File: xchebyshev_Wn_sequence.c                                             //
// Routine(s):                                                                //
//    xChebyshev_Wn_Sequence                                                  //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// void xChebyshev_Wn_Sequence(long double W[], long double x, int max_n)     //
//                                                                            //
//  Description:                                                              //
//     The Chebyshev polynomials of the fourth kind, Wn(x), are orthogonal on //
//     the interval [-1,1] with weight function w(x) = sqrt((1-x)/(1+x)) on   //
//     [-1,1] and 0 elsewhere and are normalized so that Wn(1) = 2n+1.        //
//             <Wn,Wm> = 0    if n != m,                                      //
//             <Wn,Wn> = pi   if n >= 0.                                      //
//     This routine calculates, W[n](x), the Chebyshev polynomials of the     //
//     fourth kind evaluated at x for n = 0, ..., max_n using the recursion   //
//              W[n+1](x) = 2x W[n](x) - W[n-1](x), n = 1,...,max_n - 1       //
//              W[0](x) = 1, W[1](x) = 2x+1.                                  //
//                                                                            //
//  Arguments:                                                                //
//     long double W[]                                                        //
//        On output, W[n] is the value of the Chebyshev polynomial, Wn(),     //
//        evaluated at x where 0 <= n <= max_n.  The calling routine must have//
//        defined W as long double W[N] where N >= max_n + 1.                 //
//     double x                                                               //
//        The argument of the Chebyshev polynomials Wn, n = 0,..,max_n.       //
//     int    max_n                                                           //
//        The maximum order of the sequence.                                  //
//                                                                            //
//  Return Value:                                                             //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     long double x, Wn[N];                                                  //
//     int    max_deg = N - 1;                                                //
//                                                                            //
//     (user code to set x)                                                   //
//                                                                            //
//     xChebyshev_Wn_Sequence(Wn, x, max_deg);                                //
////////////////////////////////////////////////////////////////////////////////
void xChebyshev_Wn_Sequence(long double W[], long double x, int max_n)
{
   long double two_x = x + x;
   int k;
                    // Initialize the recursion process. //

   if (max_n < 0) return; 
   W[0] = 1.0L;
   if (max_n == 0) return;
   W[1] = two_x + 1.0L;
   if (max_n == 1) return;

                   // Calculate Wn(x) for n = 2,...,max_n //

   for (k = 2; k <= max_n; k++) 
      W[k] = two_x * W[k-1] - W[k-2];

}
