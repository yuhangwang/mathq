////////////////////////////////////////////////////////////////////////////////
// File: chebyshev_Tn_sequence.c                                              //
// Routine(s):                                                                //
//    Chebyshev_Tn_Sequence                                                   //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// void Chebyshev_Tn_Sequence(double T[], double x, int max_n)                //
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
//     double T[]                                                             //
//        On output, T[n] is the value of the Chebyshev polynomial, Tn(),     //
//        evaluated at x, where 0 <= n <= max_n.  The calling routine must    //
//        have defined T as double T[N] where N >= max_n + 1.                 //
//     double x                                                               //
//        The argument of the Chebyshev polynomials Tn, n = 0,...,max_n.      //
//     int    max_n                                                           //
//        The maximum order of the sequence.                                  //
//                                                                            //
//  Return Value:                                                             //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     double x, Tn[N];                                                       //
//     int    max_deg = N - 1;                                                //
//                                                                            //
//     (user code to set x)                                                   //
//                                                                            //
//     Chebyshev_Tn_Sequence(Tn, x, max_deg);                                 //
////////////////////////////////////////////////////////////////////////////////
void Chebyshev_Tn_Sequence(double T[], double x, int max_n)
{
   long double tn, tn1, tn2;
   long double two_x = (long double) x + (long double) x;
   int k;
                    // Initialize the recursion process. //
   
   if (max_n < 0) return; 
   tn2 = 1.0L;
   tn1 = (long double) x;
   T[0] = (double) tn2;
   if (max_n == 0) return;
   T[1] = (double) tn1;
   if (max_n == 1) return;

                   // Calculate Tn(x) for n = 2,...,max_n //

   for (k = 2; k <= max_n; k++, tn2 = tn1, tn1 = tn) {
      tn = two_x * tn1 - tn2;
      T[k] = (double)tn;
   }
}
