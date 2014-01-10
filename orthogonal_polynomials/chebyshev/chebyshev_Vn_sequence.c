////////////////////////////////////////////////////////////////////////////////
// File: chebyshev_Vn_sequence.c                                              //
// Routine(s):                                                                //
//    Chebyshev_Vn_Sequence                                                   //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// void Chebyshev_Vn_Sequence(double V[], double x, int max_n)                //
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
//     double V[]                                                             //
//        On output, V[n] is the value of the Chebyshev polynomial, Vn(),     //
//        evaluated at x where 0 <= n <= max_n.  The calling routine must have//
//        defined V as double V[N] where N >= max_n + 1.                      //
//     double x                                                               //
//        The argument of the Chebyshev polynomials Vn, n = 0,...,max_n.      //
//     int    max_n                                                           //
//        The maximum order of the sequence.                                  //
//                                                                            //
//  Return Value:                                                             //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     double x, Vn[N];                                                       //
//     int    max_deg = N - 1;                                                //
//                                                                            //
//     (user code to set x)                                                   //
//                                                                            //
//     Chebyshev_Vn_Sequence(Vn, x, max_deg);                                 //
////////////////////////////////////////////////////////////////////////////////
void Chebyshev_Vn_Sequence(double V[], double x, int max_n)
{
   long double vn, vn1, vn2;
   long double two_x = (long double) x + (long double) x;
   int k;
                    // Initialize the recursion process. //
 
   if (max_n < 0) return; 
   vn2 = 1.0L;
   V[0] = (double) vn2;
   if (max_n == 0) return;
   vn1 = two_x - 1.0L;
   V[1] = (double) vn1;
   if (max_n == 1) return;

                   // Calculate Vn(x) for n = 2,...,max_n //

   for (k = 2; k <= max_n; k++, vn2 = vn1, vn1 = vn) {
      vn = two_x * vn1 - vn2;
      V[k] = (double)vn;
   }
}
