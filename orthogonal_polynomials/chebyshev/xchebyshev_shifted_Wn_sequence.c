////////////////////////////////////////////////////////////////////////////////
// File: xchebyshev_shifted_Wn_sequence.c                                     //
// Routine(s):                                                                //
//    xChebyshev_Shifted_Wn_Sequence                                          //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// void xChebyshev_Shifted_Wn_Sequence(long double W[], long double x,        //
//                                                                int max_n)  //
//                                                                            //
//  Description:                                                              //
//     The shifted Chebyshev polynomials of the fourth kind are orthogonal on //
//     the interval [0,1] with weight function w(x) = sqrt((1-x)/x).          //
//     Normalized so that W*n(1) = 2n+1.                                      //
//             <Wn*,Wm*> = 0    if n != m,                                    //
//             <Wn*,Wn*> = pi/2 if n >= 0.                                    //
//     These routines calculate W*[n](x), the shifted Chebyshev polynomials of//
//     the fourth kind evaluated at x for n = 0, ..., max_n using the         //
//     recursion                                                              //
//         W*[n+1](x) = 2(2x - 1) W*[n](x) - W*[n-1](x), n = 1,...,max_n - 1  //
//              W*[0](x) = 1, W*[1](x) = 4x - 1.                              //
//                                                                            //
//  Arguments:                                                                //
//     long double W[]                                                        //
//        On output, W[n] is the value of the shifted Chebyshev polynomial,   //
//        W*n(), evaluated at x, where 0 <= n <= max_n.  The calling routine  //
//        must have defined W as long double W[N] where N >= max_n + 1.       //
//     double x                                                               //
//        The argument of the shifted Chebyshev polynomials Wn*, n = 0,...,   //
//        max_n.                                                              //
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
//     xChebyshev_Shifted_Wn_Sequence(Wn, x, max_deg);                        //
////////////////////////////////////////////////////////////////////////////////
void xChebyshev_Shifted_Wn_Sequence(long double W[], long double x, int max_n)
{
   long double two_x_m1 = x + x - 1.0L;
   long double four_x_m2 = two_x_m1 + two_x_m1;
   int k;
                    // Initialize the recursion process. //
 
   if (max_n < 0) return;
   W[0] = 1.0L;
   if (max_n == 0) return;
   W[1] = four_x_m2 + 1.0L;
   if (max_n == 1) return;

                   // Calculate Wn(x) for n = 2,...,max_n //

   for (k = 2; k <= max_n; k++) 
      W[k] = four_x_m2 * W[k-1] - W[k-2];
}
