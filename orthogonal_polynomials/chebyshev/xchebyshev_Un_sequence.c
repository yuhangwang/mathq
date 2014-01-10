////////////////////////////////////////////////////////////////////////////////
// File: xchebyshev_Un_sequence.c                                             //
// Routine(s):                                                                //
//    xChebyshev_Un_Sequence                                                  //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// void xChebyshev_Un_Sequence(long double U[], long double x, int max_n)     //
//                                                                            //
//  Description:                                                              //
//     The Chebyshev polynomials of the second kind, Un(x), are orthogonal on //
//     the interval [-1,1] with weight function w(x) = sqrt(1-x^2) on the     //
//     interval [-1,1] and 0 elsewhere and are normalized so that Un(1) = n+1.//
//             <Un,Um> = 0    if n != m,                                      //
//             <Un,Un> = pi/2 if n >= 0.                                      //
//     This routine calculates, U[n](x), the Chebyshev polynomials of the     //
//     second kind evaluated at x for n = 0, ..., max_n using the recursion   //
//              U[n+1](x) = 2x U[n](x) - U[n-1](x), n = 1,...,max_n - 1       //
//              U[0](x) = 1, U[1](x) = 2x.                                    //
//                                                                            //
//  Arguments:                                                                //
//     long double U[]                                                        //
//        On output, U[n] is the value of the Chebyshev polynomial, Un(),     //
//        evaluated at x where 0 <= n <= max_n.  The calling routine must have//
//        defined U as long double U[N] where N >= max_n + 1.                 //
//     double x                                                               //
//        The argument of the Chebyshev polynomials Un, n = 0,...,max_n.      //
//     int    max_n                                                           //
//        The maximum order of the sequence.                                  //
//                                                                            //
//  Return Value:                                                             //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     long double x, Un[N];                                                  //
//     int    max_deg = N - 1;                                                //
//                                                                            //
//     (user code to set x)                                                   //
//                                                                            //
//     xChebyshev_Un_Sequence(Un, x, max_deg);                                //
////////////////////////////////////////////////////////////////////////////////
void xChebyshev_Un_Sequence(long double U[], long double x, int max_n)
{
   long double two_x = x + x;
   int k;

                    // Initialize the recursion process. //
   
   if (max_n < 0) return;
   U[0] = 1.0L;
   if (max_n == 0) return;
   U[1] = two_x;
   if (max_n == 1) return;

                   // Calculate Tn(x) for n = 2,...,max_n //

   for (k = 2; k <= max_n; k++) 
      U[k] = two_x * U[k-1] - U[k-2];
}
