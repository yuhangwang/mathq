////////////////////////////////////////////////////////////////////////////////
// File: xchebyshev_shifted_Un_sequence.c                                     //
// Routine(s):                                                                //
//    xChebyshev_Shifted_Un_Sequence                                          //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// void xChebyshev_Shifted_Un_Sequence(long double U[], long double x,        //
//                                                                int max_n)  //
//                                                                            //
//  Description:                                                              //
//     The shifted Chebyshev polynomials of the second kind, Un*(x)=Un(2x-1), //
//     are orthogonal on the interval [0,1] with weight function              //
//     w(x) = sqrt(x-x^2) on [0,1] and 0 elsewhere.  Further Un* is normalized//
//     so that U*n(1) = n + 1.                                                //
//             <Un*,Um*> = 0    if n != m,                                    //
//             <Un*,Un*> = pi/8 if n >= 0,                                    //
//     This routine calculates U*[n](x), the shifted Chebyshev polynomials    //
//     of the second kind evaluated at x for n = 0, ..., max_n using the      //
//     recursion:                                                             //
//        U*[k+1](x) = 2(2x - 1) U*[k](x) - U*[k-1](x), k = 1,...,max_n-1     //
//        U*[0](x) = 1, U*[1](x) = 4x - 2.                                    //
//                                                                            //
//  Arguments:                                                                //
//     long double U[]                                                        //
//        On output, U[n] is the value of the shifted Chebyshev polynomial,   //
//        U*n(), evaluated at x where 0 <= n <= max_n.  The calling routine   //
//        must have defined U as double U[N] where N >= max_n + 1.            //
//     double x                                                               //
//        The argument of the shifted Chebyshev polynomials Un*, n = 0,...    //
//        max_n.                                                              //
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
//     xChebyshev_Shifted_Un_Sequence(Un, x, max_deg);                        //
////////////////////////////////////////////////////////////////////////////////
void xChebyshev_Shifted_Un_Sequence(long double U[], long double x, int max_n)
{
   long double two_x_m1 = x + x - 1.0L;
   long double four_x_m2 = two_x_m1 + two_x_m1;
   int k;
                    // Initialize the recursion process. //
 
   if (max_n < 0) return;
   U[0] = 1.0L;
   if (max_n == 0) return;
   U[1] = four_x_m2;
   if (max_n == 1) return;

                   // Calculate Un(x) for n = 2,...,max_n //

   for (k = 2; k <= max_n; k++) 
      U[k] = four_x_m2 * U[k-1] - U[k-2];
}
