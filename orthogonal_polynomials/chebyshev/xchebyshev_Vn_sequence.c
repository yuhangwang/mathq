////////////////////////////////////////////////////////////////////////////////
// File: xchebyshev_Vn_sequence.c                                             //
// Routine(s):                                                                //
//    xChebyshev_Vn_Sequence                                                  //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// void xChebyshev_Vn_Sequence(long double V[], long double x, int max_n)     //
//                                                                            //
//  Description:                                                              //
//     The Chebyshev polynomials of the third kind, Vn(x) are orthogonal on   //
//     the interval [-1,1] with weight function w(x) = sqrt((1+x)/(1-x)) on   //
//     interval [-1,1] and 0 elsewhere and are normalized so that Vn(1) = 1.  //
//             <Vn,Vm> = 0    if n != m,                                      //
//             <Vn,Vn> = pi   if n >= 0.                                      //
//     This routine calculates V[n](x), the Chebyshev polynomials of the      //
//     third kind evaluated at x for n = 0, ..., max_n using the recursion    //
//              V[n+1](x) = 2x V[n](x) - V[n-1](x), n = 1,...,max_n - 1       //
//              V[0](x) = 1, V[1](x) = 2x-1.                                  //
//                                                                            //
//  Arguments:                                                                //
//     long double V[]                                                        //
//        On output, V[n] is the value of the Chebyshev polynomial, Vn(),     //
//        evaluated at x where 0 <= n <= max_n.  The calling routine must have//
//        defined V as long double V[N] where N >= max_n + 1.                 //
//     double x                                                               //
//        The argument of the Chebyshev polynomials Vn, n = 0,...,max_n.      //
//     int    max_n                                                           //
//        The maximum order of the sequence.                                  //
//                                                                            //
//  Return Value:                                                             //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     long double x, Vn[N];                                                  //
//     int    max_deg = N - 1;                                                //
//                                                                            //
//     (user code to set x)                                                   //
//                                                                            //
//     xChebyshev_Vn_Sequence(Vn, x, max_deg);                                //
////////////////////////////////////////////////////////////////////////////////
void xChebyshev_Vn_Sequence(long double V[], long double x, int max_n)
{
   long double two_x = x + x;
   int k;
                    // Initialize the recursion process. //

   if (max_n < 0) return; 
   V[0] = 1.0L;
   if (max_n == 0) return;
   V[1] = two_x - 1.0L;
   if (max_n == 1) return;

                   // Calculate Vn(x) for n = 2,...,max_n //

   for (k = 2; k <= max_n; k++) 
      V[k] = two_x * V[k-1] - V[k-2];
}
