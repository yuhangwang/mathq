////////////////////////////////////////////////////////////////////////////////
// File: chebyshev_shifted_Vn_sequence.c                                      //
// Routine(s):                                                                //
//    Chebyshev_Shifted_Vn_Sequence                                           //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// void Chebyshev_Shifted_Vn_Sequence(double V[], double x, int max_n)        //
//                                                                            //
//  Description:                                                              //
//     The shifted Chebyshev polynomials of the third kind, Vn*(x) = Vn(2x-1),//
//     are orthogonal on the interval [0,1] with weight function              //
//     w(x) = sqrt(x/(1-x)) on [0,1] and 0 elsewhere.  Further Vn* is         //
//     normalized so that Vn*(1) = 1.                                         //
//             <Vn*,Vm*> = 0    if n != m,                                    //
//             <Vn*,Vn*> = pi/2 if n >= 0.                                    //
//     This routine calculates V*[n](x), the shifted Chebyshev polynomials of //
//     the third kind evaluated at x for n = 0, ..., max_n using the recursion//
//         V*[n+1](x) = 2(2x - 1) V*[n](x) - V*[n-1](x), n = 1,...,max_n - 1  //
//              V*[0](x) = 1, V*[1](x) = 4x - 3.                              //
//                                                                            //
//  Arguments:                                                                //
//     double V[]                                                             //
//        On output, V[n] is the value of the shifted Chebyshev polynomial,   //
//        V*n(), evaluated at x where 0 <= n <= max_n.  The calling routine   //
//        must have defined V as double V[N] where N >= max_n + 1.            //
//     double x                                                               //
//        The argument of the shifted Chebyshev polynomials Vn*, n = 0,...,   //
//        max_n.                                                              //
//     int    max_n                                                           //
//        The maximum order of the sequence.                                  //
//                                                                            //
//  Return Value:                                                             //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     double Vn[N];                                                          //
//     double x;                                                              //
//     int    max_deg = N - 1;                                                //
//                                                                            //
//     (user code to set x)                                                   //
//                                                                            //
//     Chebyshev_Shifted_Vn_Sequence(Vn, x, max_deg);                         //
////////////////////////////////////////////////////////////////////////////////
void Chebyshev_Shifted_Vn_Sequence(double V[], double x, int max_n)
{
   long double vn, vn1, vn2;
   long double two_x_m1 = (long double) x + (long double) x - 1.0L;
   long double four_x_m2 = two_x_m1 + two_x_m1;
   int k;
                    // Initialize the recursion process. //
 
   if (max_n < 0) return;
   vn2 = 1.0L;
   V[0] = (double) vn2;
   if (max_n == 0) return;
   vn1 = (long double) four_x_m2 - 1.0L;
   V[1] = (double) vn1;
   if (max_n == 1) return;

                   // Calculate Vn(x) for n = 2,...,max_n //

   for (k = 2; k <= max_n; k++, vn2 = vn1, vn1 = vn) {
      vn = four_x_m2 * vn1 - vn2;
      V[k] = (double)vn;
   }
}
